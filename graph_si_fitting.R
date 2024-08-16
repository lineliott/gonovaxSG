library(ggplot2)
source("script_calibration_support.R")
source("script_calibration_plot.R")
stat_boxplot_custom <- function(mapping = NULL, data = NULL,
                                geom = "boxplot", position = "dodge",
                                ...,
                                qs = c(0.025, 0.025, 0.5, 0.975, 0.975),
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatBoxplotCustom,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      qs = qs,
      ...
    )
  )
}

StatBoxplotCustom <- ggproto("StatBoxplotCustom", Stat,
                             required_aes = c("x", "y"),
                             non_missing_aes = "weight",
                             setup_params = function(data, params) {
                               params$width <- ggplot2:::"%||%"(
                                 params$width, (resolution(data$x) * 0.75)
                               )
                               
                               if (is.double(data$x) && !ggplot2:::has_groups(data) && any(data$x != data$x[1L])) {
                                 warning(
                                   "Continuous x aesthetic -- did you forget aes(group=...)?",
                                   call. = FALSE
                                 )
                               }
                               
                               params
                             },
                             
                             compute_group = function(data, scales, width = NULL, na.rm = FALSE, qs = c(.05, .25, 0.5, 0.75, 0.95)) {
                               
                               if (!is.null(data$weight)) {
                                 mod <- quantreg::rq(y ~ 1, weights = weight, data = data, tau = qs)
                                 stats <- as.numeric(stats::coef(mod))
                               } else {
                                 stats <- as.numeric(stats::quantile(data$y, qs))
                                 stats[3] <- mean(data$y) # replace median by mean
                               }
                               names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
                               iqr <- diff(stats[c(2, 4)])
                               
                               outliers <- (data$y < stats[1]) | (data$y > stats[5])
                               
                               if (length(unique(data$x)) > 1)
                                 width <- diff(range(data$x)) * 0.9
                               
                               
                               df <- as.data.frame(as.list(stats))
                               df$outliers <- list(data$y[outliers])
                               
                               if (is.null(data$weight)) {
                                 n <- sum(!is.na(data$y))
                               } else {
                                 # Sum up weights for non-NA positions of y and weight
                                 n <- sum(data$weight[!is.na(data$y) & !is.na(data$weight)])
                               }
                               
                               df$notchupper <- df$middle + 1.58 * iqr / sqrt(n)
                               df$notchlower <- df$middle - 1.58 * iqr / sqrt(n)
                               
                               df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
                               df$width <- width
                               df$relvarwidth <- sqrt(n)
                               df
                             }
)





stream <- c("sgmsm","sgmale","ukmsm")
projection_combined <- list()
observation_combined <- list()
for (i in stream) {
  if (i == "sgmsm") {
    sg_data <- read_csv("MSM.csv")
    mcmc_sample <- readRDS("calibration_sgmsm/calibration/outputs/mcmc_sample_tuned.rds")
    title <- "Main scenario"
  } else if (i == "sgmale"){
    sg_data <- read_csv("male.csv")
    mcmc_sample <- readRDS("calibration_sgmale/calibration/outputs/mcmc_sample_tuned.rds")
    title <- "Upper bound"
  } else if (i == "ukmsm"){
    sg_data <- read_csv("male.csv")
    sg_data$diagnosed <- sg_data$diagnosed*0.7
    mcmc_sample <- readRDS("calibration_ukmsm/calibration/outputs/mcmc_sample_tuned.rds")
    title <- "Lower bound"
  }
  
  data <- sg_data
  n <- nrow(mcmc_sample$pars)
  y <- lapply(seq_len(n), function(i) {
    gonovax::run(tt = c(0, data$gonovax_year),
                 gono_params = mcmc_sample$transform(mcmc_sample$pars[i, ]),
                 demographic_params = demographic_params_sg,
                 transform = FALSE)
  })
  projection <- list(extract_trajectories(y, "tot_treated"))
  projection <- reshape2::melt(projection)
  projection$L1 <- title
  colnames(projection) <- c("year", "sample", "value", "scenario")
  year_mapping <- 2004:2018
  projection$year <- rep(year_mapping, length.out = nrow(projection) / length(unique(projection$sample)))
  projection$type <- "Projection"
  # Add the projection to the combined list
  projection_combined[[i]] <- projection
  
  # observation 
  observation <- sg_data[,2:3]
  observation$scenario <- title
  observation$type <- "Observation"
  observation_combined[[i]] <- observation
  }
# Combine all projections into one data frame
projection_combined <- do.call(rbind, projection_combined)
observation_combined <- do.call(rbind,observation_combined)

# Specify the strategy order
order <- c("Main scenario", "Upper bound", "Lower bound")
projection_combined$scenario <- factor(projection_combined$scenario, levels = order)
observation_combined$scenario <- factor(observation_combined$scenario, levels = order)
  
p <- ggplot() +
    # Boxplot with custom statistics
    stat_boxplot_custom(data = projection_combined, 
                        aes(x = factor(year), 
                            y = value, 
                            color = factor(type), 
                            fill = factor(type)), 
                        qs = c(0.025, 0.025, 0.5, 0.975, 0.975),
                        outlier.shape = NA,
                        alpha = 0.7) +
    # stat_summary(fun = mean, 
    #              geom ="point", 
    #              shape = 18, 
    #              size = 10) +
    # Points for the observations
    geom_point(data = observation_combined, 
               aes(x = factor(year), 
                   y = diagnosed, 
                   color = factor(type), 
                   fill = factor(type)), 
                   stroke = 1.5,
                   shape = 23, 
                   size = 6) +
    ylab("Annual incidence of diagnosed") +
    xlab("Year")+
    scale_color_manual(name = "Type",
                       label = c("Projection", "Observation"),
                       values = c(Projection = "darkred",
                                  Observation = "cadetblue1")) +
    scale_fill_manual(name = "Type",
                      label = c("Projection", "Observation"),
                      values = c(Projection = "darksalmon",
                                 Observation = "deepskyblue4")) +
    scale_x_discrete(breaks = c(2004, 2009, 2013, 2018)) + 
    scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 1000)) +
    coord_capped_cart(bottom='both', left='both') +
    facet_rep_wrap(~scenario, ncol = 6, scales = "free_x") + 
    theme(
      axis.title.y = element_text(size = 30),
      axis.title.x = element_text(size = 30),
      axis.text.x = element_text(size = 30, hjust = 0.5),
      axis.text.y = element_text(size = 30, hjust = 0.5),
      plot.title = element_blank(),
      panel.background = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 30, hjust = 0.5),
      panel.spacing = unit(50, "pt"),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 30),
      panel.border = element_blank(),
      axis.line = element_line()
    )

ggsave(paste0("SI/fitting.pdf"), plot = p, width = 24, height = 8)
