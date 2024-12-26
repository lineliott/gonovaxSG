# these functions are developed, adapted/cited from gonovax package in support of plotting calibration result
library(ggplot2)
library(lemon)
library(gridExtra)


write_png <- function(filename, code, asp = 5 / 8, res = 200, pointsize = 8,
                      width = 17, ...) {
  png(filename, width = width, units = "cm", height = width * asp, res = res,
      pointsize = pointsize, ...)
  options(scipen = 5)
  par(bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), oma = rep(1, 4))
  on.exit(dev.off())
  force(code)
}

plot_axis <- function(at, samples) {
  plot(0, 0, type = "n", ylim = c(0, max(at)), xlim = c(0, 1), axes = FALSE,
       ylab = "", xlab = "")
  axis(2, pos = 1, at, names(samples), las = 1, tick = FALSE, font = 3)
}

plot_bars <- function(samples, what, at, add = FALSE, xlab = NULL, xmax = NULL,
                      col = "grey20", parameters_prior) {
  x <- lapply(samples, "[", , what)
  if (!add) {
    xmax <- xmax %||% max(unlist(x), na.rm = TRUE)
    xlab <- xlab %||% what
    plot(0, 0, type = "n", ylim = c(0, max(at)), xlim = c(0, xmax), yaxt = "n",
         ylab = "", xlab = xlab)
  }
  
  type <- parameters_prior[what, "type"]
  
  if (type %in% c("beta", "gamma", "lnorm")) {
    qs <- parameters_prior[what, sprintf("%s_%s", type, c("mean", "q2.5", "q97.5"))]
    abline(v = qs, lty = c(1, 2, 2), col = "grey20")
  }
  
  mapply(ci_bar, x, at, width = 0.1, col = col, horiz = TRUE)
}

plot_bars_t <- function(samples, what, at, xlab = NULL, cols, xmax = NULL) {
  plot_bars(samples, paste0(what, 2), at + 0.1, FALSE, xlab, xmax, cols[3])
  plot_bars(samples, what, at, TRUE, col = cols[2])
  plot_bars(samples, paste0(what, 1), at - 0.1, TRUE, col = cols[1])
}

compare_models <- function(samples, parameters_prior) {
  at <- seq_along(samples)
  cols <- hcl.colors(3)
  plot_axis(at, samples)
  
  MoreArgs <- list(samples = samples, at = at, col = cols[2],
                   parameters_prior = parameters_prior)
  mapply(plot_bars,
         what = c("prev_Aslpc", "prev_Ashpc","beta_2004", "beta2020", "epsilon"),
         xmax = c(2, 10, 1, 1, 1),
         MoreArgs = MoreArgs)
  # plot_bars_t(reformatted_samples, "eta", at, xlab = "eta", xmax = 0.5, cols)
  
  plot_axis(at, samples)
  mapply(plot_bars,
         what = c( "psipc", paste0("1/",c("sigma", "nu", "mu", "rho"))),
         xmax = c( 100, 10, 365 * 2, 10, 10),
         MoreArgs = MoreArgs)
  
  plot_axis(at, samples)
  mapply(plot_bars,
         what = c("eta_l", "eta_h", "phi_eta_l", "phi_eta_h", "omega"),
         xmax = c(1, 1, 2, 2, 1),
         MoreArgs = MoreArgs)
  
}

plot_posterior <- function(mcmc_sample, what, xlab = NULL, col = "#D7CADEFF") {
  
  xlab <- xlab %||% what
  x <- mcmc_sample$pars[, what]
  qs <- quantile(x, c(0.5, 0.025, 0.975))
  hist(x, freq = FALSE, col = col, xlab = xlab, main = "", breaks = 20)
  
  mtext(sprintf("%.2g (%.2g, %.2g)", qs[1], qs[2], qs[3]),
        side = 3, line = 1, cex = 0.8, font = 3)
}

plot_corr <- function(sample, x, y, n = 100, palette = "Mint") {
  lp <- sample$probabilities[, "log_posterior"]
  pars <- data.frame(sample$pars[order(lp, decreasing = FALSE), ],
                     col = ceiling(seq(1e-10, n, length.out = length(ll))))
  cols <- hcl.colors(n, palette)
  plot(pars[, x], pars[, y], pch = 20, col = cols[pars$col], xlab = x, ylab = y)
}

plot_pairs <- function(sample) {
  
  pars <- colnames(sample$pars)
  n <- length(pars)
  corr <- cor(sample$pars)
  
  par(mfcol = c(n, n))
  for (x in seq_len(n)) {
    for (y in seq_len(n)) {
      if (x == y) {
        plot_posterior(sample, pars[x])
      } else if (x < y) {
        plot_corr(sample, pars[x], pars[y])
      } else {
        plot.new()
        text(0.5, 0.5, sprintf("r = %.2g", corr[pars[x], pars[y]]), cex = 1.2)
      }
    }
  } 
}

plot_traces <- function(mcmc_results, mcmc_metrics, skip_spot = FALSE) {
  
  n_chains <- max(mcmc_results$chain)
  col <- hcl.colors(n_chains, "Geyser", rev = TRUE)
  
  par_names <- ylab <- colnames(mcmc_results$pars)
  ylab[1:2] <- paste(ylab[1:2], "(%)")
  fac <- ifelse(ylab == par_names, 1, 100)
  
  Map(plot_trace, what = par_names, ylab = ylab, fac = fac,
      MoreArgs = list(mcmc_results = mcmc_results,
                      mcmc_metrics = mcmc_metrics,
                      col = col, title = TRUE))
  plot.new()
  legend("topleft", fill = col, legend = sprintf("chain%i", seq_len(n_chains)),
         bty = "n", ncol = 2)
  ar <- mcmc_metrics$acceptance_rate[[1]] * 100
  mtext(sprintf("AR = %.2fpc", ar),
        side = 3, line = 1, cex = 0.8, font = 3)

  if (skip_spot) {
    plot.new()
  }

  Map(plot_trace, what = colnames(mcmc_results$probabilities),
      MoreArgs = list(mcmc_results = mcmc_results,
                      mcmc_metrics = mcmc_metrics,
                      col = col))
  
  
}

plot_trace <- function(mcmc_results, mcmc_metrics, what, ylab = NULL, col,
                       title = FALSE, fac = 1) {
  
  results <- cbind(mcmc_results$pars, mcmc_results$probabilities)
  y <- matrix(results[, what], ncol = max(mcmc_results$chain))
  ylab <- ylab %||% what
  
  matplot(y, type = "l", col = col, xlab = "Iteration", ylab = ylab, lty = 1)
  
  if (title) {
    qs <- quantile(y, c(0.5, 0.025, 0.975)) * fac
    ess <- mcmc_metrics$ess[what]
    gr <- mcmc_metrics$gr$psrf[what, 2]
    
    mtext(sprintf("%s = %.3g (%.3g, %.3g);\nESS = %.f, GR = %.2f",
                  ylab, qs[1], qs[2], qs[3], ess, gr),
          side = 3, line = 0.5, cex = 0.7, font = 3)
  }
}

plot_posteriors <- function(mcmc_sample, priors) {
  Map(plot_posterior, what = colnames(mcmc_sample$pars),
      MoreArgs = list(mcmc_sample = mcmc_sample, priors = priors))
}

plot_posterior <- function(mcmc_sample, priors, what,xlab = NULL) {
  
  xlab <- xlab %||% what
  x <- mcmc_sample$pars[, what]
  qs <- quantile(x, c(0.5, 0.025, 0.975))
  hist(x, freq = FALSE, col = "#D7CADEFF", xlab = xlab, main = "", breaks = 20) 
       #xlim = c(min(parameters_info[parameters_info$name==what,]$min), 
       #        max(x)))
  
  prior <- function(x) {
    f <- priors[[what]] %||% function(x) pmin(x, 1e-9)
    exp(f(x))
  }
  
  curve(prior(x), min(x), max(x), 100, TRUE, col = grey(0.2), lwd = 1.5)
  
  mtext(sprintf("%.3g (%.3g, %.3g)", qs[1], qs[2], qs[3]),
        side = 3, line = 1, cex = 0.8, font = 3)
  
}


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



boxplot_fits_eliott <- function(data,mcmc_sample) {
  
  data <- sg_data
  n <- nrow(mcmc_sample$pars)
  y <- lapply(seq_len(n), function(i) {
        gonovax::run(tt = c(0, data$gonovax_year),
        gono_params = mcmc_sample$transform(mcmc_sample$pars[i, ]),
        demographic_params = demographic_params_sg,
        transform = FALSE)
  })
  
  
  #diag_a <- extract_trajectories(y, "cum_diag_a[1,1]") +
  #  extract_trajectories(y, "cum_diag_a[2,1]")
  # diag_s <- extract_trajectories(y, "cum_diag_s[1,1]") +
  #   extract_trajectories(y, "cum_diag_s[2,1]")
  # data$p_symp <- data$n_symptomatic / data$n_reported * 100
  
  fits <- list(diagnosed = extract_trajectories(y, "tot_treated")
               # attended  = extract_trajectories(y, "tot_attended"), 
               # p_symp = diag_s / (diag_a + diag_s) * 100
               )
  
  # Convert 'fits' to a data frame suitable for ggplot2
  fits_df <- reshape2::melt(fits$diagnosed)
  colnames(fits_df) <- c("year", "sample", "value")
  fits_df$year <- 2004:2018

  #boxplot_fit(fits_df, "diagnosed", data)

  fits$year <- factor(fits$year, levels = 2004:2018)
  # Create a dummy variable for legend assignment
  fits_df$legend <- "Projection"
  sg_data$legend <- "Observation"
  
  p<-ggplot() +
    # Boxplot with custom statistics
    stat_boxplot_custom(data = fits_df, 
                        aes(x = factor(year), 
                            y = value, 
                            color = legend, 
                            fill = legend), 
                        qs = c(0.025, 0.025, 0.5, 0.975, 0.975),
                        outlier.shape = NA,
                        alpha = 0.7) +
    # Points for the data
    geom_point(data = sg_data, 
               aes(x = factor(year), y = diagnosed, 
                   color = legend, 
                   fill = legend), 
               shape = 18, size = 10) +
    ylab("Annual incidence of diagnosed") +
    xlab("Year")+
    scale_color_manual(name = "Type",
                       values = c("Projection" = "darkred",
                                  "Observation" = "deepskyblue4")) +
    scale_fill_manual(name = "Type",
                      values = c("Projection" = "darksalmon",
                                 "Observation" = "cadetblue1")) +
    scale_x_discrete(breaks = c(2004, 2009, 2013, 2018)) + 
    scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 1000)) +
    coord_capped_cart(bottom = 'both', left = 'both') +
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
  
  
    return(p)
  
# plot_fit(fits, "diagnosed", data, ylab = "Total diagnoses",
#          ylim = c(0, 1e4))
#  plot_fit(fits, "attended", data, ylab = "Total tests performed",
#           ylim = c(0, 3e5))
#  plot_fit(fits, "p_symp", data, ylab = "Diagnoses with symptoms (%)",
#           ylim = c(0, 100))
}


extract_trajectories <- function(y, what) {
  apply(sapply(y, "[", , what), 2, diff)
}
# plot_fit <- function(fits, what, data, ...) {
#   matplot(data$year, fits[[what]], type = "l", xlab = "",
#           lty = 1, col = grey(0.8), xlim = c(data$year[1], 2020), ...)
#   points(data$year, data[[what]], bg = "darkred", pch = 23)
# }



traceplot_fits_eliott <- function(data, mcmc_sample) {
  
  data <- sg_data
  n <- nrow(mcmc_sample$pars)
  y <- lapply(seq_len(n), function(i) {
    gonovax::run(tt = c(0, data$gonovax_year),
        gono_params = mcmc_sample$transform(mcmc_sample$pars[i, ]),
        demographic_params = demographic_params_sg,
        transform = FALSE)
  })
  
  fits <- list(diagnosed = extract_trajectories(y, "tot_treated"))
  
  # Convert 'fits' to a data frame suitable for ggplot2
  fits_df <- reshape2::melt(fits$diagnosed)
  colnames(fits_df) <- c("year", "sample", "value")
  fits_df$year <- 2004:2018
  
  # Calculate quantiles for each year
  fits_df_quantiles <- aggregate(value ~ year, data = fits_df, function(x) {
    c(Q1 = quantile(x, probs = 0.25), Q2 = median(x), Q3 = quantile(x, probs = 0.75))
  })
  fits_df_quantiles <- do.call(data.frame, fits_df_quantiles)
  colnames(fits_df_quantiles) <- c("year", "Q1", "Q2","Q3")
  
  # Merge quantiles with fits_df
  fits_df <- merge(fits_df, fits_df_quantiles, by = "year")
  
  traceplot_fit(fits_df, "diagnosed", data)
  
}

traceplot_fit <- function(fits_df, what, data) {
  p <- ggplot(fits_df, aes(x = factor(year), y = Q2)) +
    geom_ribbon(aes(ymin = Q1, ymax = Q3, group = sample), fill='grey', alpha=0.5) +
    geom_line(aes(group = sample), alpha = 0.5, linewidth=1) +
    geom_point(data = data, aes(x = factor(year), y = .data[[what]]), color = "darkred", size=3) +
    scale_x_discrete(breaks = seq(from = 2004, to = 2018, by = 4)) +
    scale_y_continuous(limits = c(NA, 5000)) +  # Set the upper limit to 6000
    labs(x = "Year", y = "Total Diagnosed") +
    theme_minimal()+
    theme(axis.title = element_text(#face = "bold", 
                                    size = 20),
          axis.text.x = element_text(#face = "bold", 
                                     size = 20),
          axis.text.y = element_text(#face = "italic", 
                                     size = 20))
  return(p)
}

# Define the combined function for plotting MCMC traces and metrics
plot_mcmc_traces_ggplot <- function(mcmc_results, mcmc_metrics, path) {
  n_chains <- max(mcmc_results$chain)
  col <- hcl.colors(n_chains, "Geyser", rev = TRUE)
  
  par_names <- colnames(mcmc_results$pars)
  ylab <- par_names
  ylab[1:4] <- paste(ylab[1:4], "(%)")
  fac <- ifelse(ylab == par_names, 1, 100)
  
  plots <- list()
  
  # Generate trace plots for parameters
  for (i in seq_along(par_names)) {
    qs <- quantile(mcmc_results$pars[, i], c(0.5, 0.025, 0.975)) * fac[i]
    mean <- mean(mcmc_results$pars[, i]) * fac[i]
    ess <- mcmc_metrics$ess[i]
    gr <- mcmc_metrics$gr$psrf[i, 2]
    title_text <- sprintf("%s = %.3g (%.3g, %.3g);\nESS = %.f, GR = %.2f",
                          ylab[i], mean, qs[2], qs[3], ess, gr)
    
    # Create a sequence for the iteration and chain assignment
    n_iterations <- length(mcmc_results$pars[, i]) / n_chains
    iteration_seq <- rep(seq(n_iterations), n_chains)
    chain_seq <- rep(seq(n_chains), each = n_iterations)
    
    # Create the data frame for plotting and filter off first 1000 iterations
    plot_data <- data.frame(
      iteration = iteration_seq,
      value = mcmc_results$pars[, i],
      chain = chain_seq
    ) %>%
      group_by(chain) %>%
      filter(iteration > 1000) %>%
      ungroup()
    
    plot <- ggplot(data = plot_data, 
                   aes(x = iteration/100000, y = value, color = factor(chain))) +
      geom_line() +
      labs(x = "Iteration per 100,000", y = ylab[i], color = "Chain") +
      ggtitle(title_text) +
      scale_color_manual(values = col) +
      theme(axis.title = element_text(size = 30),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            plot.title = element_text(size = 30, hjust = 0.5),
            legend.text = element_text(size = 30),  # Increase legend text size
            legend.title = element_text(size = 30),  # Increase legend title size
            panel.background = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(size = 30), # pseudo x axis title
            panel.spacing = unit(50, "pt"),
            panel.border = element_blank(),
            axis.line = element_line(),
            plot.margin = margin(t = 10, r = 10, b = 20, l = 10))
    
    plots[[length(plots) + 1]] <- plot
  }
  
  # Generate trace plots for probabilities
  for (j in colnames(mcmc_results$probabilities)) {
    
    # Create a sequence for the iteration and chain assignment
    n_iterations <- length(mcmc_results$probabilities[, j]) / n_chains
    iteration_seq <- rep(seq(n_iterations), n_chains)
    chain_seq <- rep(seq(n_chains), each = n_iterations)
    
    # Create the data frame for plotting and filter off first 1000 iterations
    plot_data <- data.frame(
      iteration = iteration_seq,
      value = mcmc_results$probabilities[, j],
      chain = chain_seq
    ) %>%
      group_by(chain) %>%
      filter(iteration > 1000) %>%
      ungroup()
    
    
    ar <- mcmc_metrics$acceptance_rate[[1]] * 100
    title_text <- sprintf("%s\nAcceptance Ratio = %.2f%%", j, ar) 
    plot <- ggplot(data = plot_data,
                   aes(x = iteration/100000, y = value, color = factor(chain))) +
      geom_line() +
      labs(x = "Iteration per 100,000", y = j, color = "Chain") +
      scale_color_manual(values = col) +
      ggtitle(title_text)+
      theme(axis.title = element_text(size = 30),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            plot.title = element_text(size = 30, hjust = 0.5),
            legend.text = element_text(size = 30),  # Increase legend text size
            legend.title = element_text(size = 30),  # Increase legend title size
            panel.background = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(size = 30), # pseudo x axis title
            panel.spacing = unit(50, "pt"),
            panel.border = element_blank(),
            axis.line = element_line(),
            plot.margin = margin(t = 10, r = 10, b = 20, l = 10))
    
    plots[[length(plots) + 1]] <- plot
  }
  
  # Add a legend plot
  legend_plot <- ggplot(data = data.frame(chain = factor(seq_len(n_chains))), aes(x = 1, fill = chain)) +
    geom_bar() +
    scale_fill_manual(values = col) +
    theme_void() +
    theme(legend.position = "none") +
    labs(fill = "Chain")
  
  # Function to modify the plot based on its position
  adjust_theme <- function(plot, position, ncol, nrow) {
    row <- ceiling(position / ncol)
    col <- position %% ncol
    col <- ifelse(col == 0, ncol, col)
    
    plot <- plot +
      theme(
        legend.position = if (col == ncol) "right" else "none",
        axis.text.x = if (row == nrow) element_text(size = 30) else element_blank(),
        axis.title.x = if (row == nrow) element_text(size = 30) else element_blank()
      )
    
    return(plot)
  }
  
  # Adjust the theme for each plot based on its position in the grid
  ncol <- 4
  nrow <- ceiling(length(plots) / ncol)
  
  for (i in seq_along(plots)) {
    plots[[i]] <- adjust_theme(plots[[i]], i, ncol, nrow)
  }
  
  # Save the plots to a PDF file with a 3x5 layout
  pdf(paste0(path,".pdf"), width = 48, height = 32)
  grid.arrange(grobs = plots, ncol = 4)
  dev.off()
  
  # Save the plots to a PNG file with 320 DPI
  png(paste0(path,".png"), width = 48, height = 32, units = "in", res = 320)
  grid.arrange(grobs = plots, ncol = 4)
  dev.off()
}

