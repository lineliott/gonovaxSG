# This script produces figure 3 (heatmap: sensitivity analysis)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lemon)
library(patchwork)
library(showtext)
showtext_auto()
path <- "calibration/graph/multivax/bivariates_10/"
dir.create(path, FALSE, TRUE)


# generate statistical summary
stats_summary <- function(vector) {
  if (anyNA(vector)){
    return(c("NaN produced","NaN produced","NaN produced","NaN produced","NaN produced","NaN produced"))
  }else{
    ci <- quantile(vector, c(0.025, 0.5, 0.975))
    ci["mean"] <- mean(vector)  # add on mean
    ci["min"] <- min(vector) # min
    ci["max"] <- max(vector) # max
    ci <- round(ci, 4)  # Round to 2 decimals
    return(ci) # 2.5;50;97.5;mean;min;max
  }
}

run_grid_control <- function(short_run, title, revax) {
  vbe_uptake <- 0 # no vbe background
  
  if (short_run) {
    eff <- c(0.31)
    dur <- c(4)
  } else {
    eff <- seq(0.1, 1, 0.1)
    dur <- seq(1, 10, 1)
  }
  if (revax) {
    model <- run_onevax_xvwv
  } else {
    stop("no-revaccination not supported")
  }
  
  if (title == "VbE") {
    vbe_uptake <- 0.873 # vbe strategy
  }
  
  list(eff = eff, 
       dur = dur,
       model = model, vbe = vbe_uptake)
}

strategies <- c("VbE", "VoD(H)", "VoD","VoA",  "VoA(H)", "VoD+VoA(H)")
for (i in strategies) {
  short_run <- FALSE
  uptake_total<- 0.33
  graph <- "heatmap"
  elapsed <- 10
  vaccine_algorithm <- "xvwv"
  stream <- "msm"
  fix_par_t <- TRUE
  strategy <- if (i == "VbE") NULL else i
  # run result
  source("script_projection.R")
  # process
  
  # Create empty vectors to store results
  eff <- numeric(length(intervention[["results"]]))
  dur <- numeric(length(intervention[["results"]]))
  averted_cases <- numeric()
  averted_cases_per_dose <- numeric()
  averted_cases_percentage <- numeric()
  mean_averted_case <- numeric(length(intervention[["results"]]))
  
  mean_perdose <- numeric(length(intervention[["results"]]))
  mean_percentage <- numeric(length(intervention[["results"]]))
  front_perdose <- numeric(length(intervention[["results"]]))
  front_percentage <- numeric(length(intervention[["results"]]))
  tail_perdose <- numeric(length(intervention[["results"]]))
  tail_percentage <- numeric(length(intervention[["results"]]))
  
  
  # Iterate through each list element
  for (i in seq_along(intervention[["results"]])) {
    result_key <- names(intervention[["results"]][i])
    # Split names into eff and dur
    eff_dur <- strsplit(result_key, "_dur")
    eff[i] <- sapply(eff_dur, function(x) as.numeric(sub("eff", "", x[1])))
    dur[i] <- sapply(eff_dur, function(x) as.numeric(sub("_dur", "", x[2])))
    
    # Retrieve wanted metrics by samples in 10 year expectation
    cases_averted_per_dose <- intervention[["results"]][[i]][["cases_averted_per_dose"]][elapsed,]
    inc_cum_doses <- intervention[["results"]][[i]][["inc_cum_doses"]][elapsed,]
    averted_case <- cases_averted_per_dose * inc_cum_doses
    averted_cases_percentage <- ((novax_2044_flows[["cum_treated"]][elapsed,]-intervention[["results"]][[i]][["cum_treated"]][elapsed,])/novax_2044_flows[["cum_treated"]][elapsed,])
    if (any(is.na(averted_case))) {
      averted_case <- inc_cum_doses
    }
    
    if (title == "VbE") {
      singapore_male_entrance <- demographic_params_sg$male_enr
      vaccinated <- singapore_male_entrance*elapsed*control$vbe
      averted_cases <- -intervention[["results"]][[i]][["inc_cum_treated"]][elapsed,]
      cases_averted_per_dose <- averted_cases/(vaccinated*2) 
    }
    # averted_cases <- c(averted_cases, averted_case)
    # averted_cases_per_dose <- c(averted_cases_per_dose, cases_averted_per_dose)
    # averted_cases_percentage <- c(averted_cases_percentage, ((novax_2044_flows[["cum_treated"]][10,]-intervention[["results"]][[i]][["cum_treated"]][10,])/novax_2044_flows[["cum_treated"]][10,])*100)
    
    # Calculate stats
    perdose <- stats_summary(cases_averted_per_dose)
    percentage <- stats_summary(averted_cases_percentage)
    
    # mean_averted_case[i] <- mean(averted_case)
    #front_perdose[i] <- perdose [[1]]
    mean_perdose[i] <- perdose [[4]]
    #tail_perdose[i] <- perdose [[3]]
    
    #front_percentage[i] <- percentage[[1]]
    mean_percentage[i] <- percentage [[4]]
    #tail_percentage[i] <- percentage [[3]]
  }
  
  
  # Create the data frame with the modified object naming convention
  results_summary <- data.frame(eff = eff, 
                                dur = dur, 
                               # uptake = uptake_total,
                                mean_perdose = mean_perdose,
                                mean_percentage = mean_percentage,
                                strategy = title)
  
  
  # save data frame according to the naming convention
  object_name <- paste(path, title, "_multivax_summary.rds", sep = "")
  saveRDS(results_summary, object_name)
  print(paste0(object_name, " analysis finished"))
  }
###########################################################################
VoA_multivax_summary <- readRDS(paste0(path,"VoA_multivax_summary.rds"))
`VoA(H)_multivax_summary` <- readRDS(paste0(path,"VoA(H)_multivax_summary.rds"))
`VaR_multivax_summary` <- readRDS(paste0(path,"VaR_multivax_summary.rds"))
VoD_multivax_summary <- readRDS(paste0(path,"VoD_multivax_summary.rds"))
`VoD(H)_multivax_summary` <- readRDS(paste0(path,"VoD(H)_multivax_summary.rds"))
#novax_multivax_summary <- readRDS(paste0(path,"Heatmap_novax_multivax_summary.rds"))
VbE_multivax_summary <- readRDS(paste0(path,"VbE_multivax_summary.rds"))

combined_df <- bind_rows(VoA_multivax_summary,
                         `VoA(H)_multivax_summary`,
                         `VaR_multivax_summary`, 
                         VoD_multivax_summary, 
                         `VoD(H)_multivax_summary`,
                        # novax_multivax_summary, 
                         VbE_multivax_summary, .id = "id")

banchmark <- subset(combined_df, combined_df$eff==0.3 & combined_df$dur==4)

percentage_mark <- 0.4026
perdose_mark <- 0.24

perdose_superior <- subset(combined_df,combined_df$mean_perdose > perdose_mark)
percentage_superior <- subset(combined_df, combined_df$mean_percentage > percentage_mark)




# Specify the strategy order
strategy_order <- c("VbE", "VoD(H)","VoD", "VoA", "VoA(H)", "VaR")
combined_df$strategy <- factor(combined_df$strategy, levels = strategy_order)
percentage_superior$strategy <- factor(percentage_superior$strategy, levels = strategy_order)
perdose_superior$strategy <- factor(perdose_superior$strategy, levels = strategy_order)



plot1 <- ggplot(combined_df, aes(x = eff, y = dur, fill = mean_percentage*100)) +
  geom_tile() + # must before add on point
  geom_point(data = percentage_superior, aes(x = eff, y = dur, shape = paste0("≥ ", round(percentage_mark*100, 2),"%")), size = 10,  color = "darkgoldenrod1") +
  facet_rep_wrap(~strategy, ncol = 6, scales = "fixed", 
                 strip.position = "top") + 
  scale_x_continuous(limit = c(0, 1.05), breaks = seq (0.1,1,0.3) )+
  # scale_fill_gradient(name = "A. P.\n", 
  #                     low = "white", high = "deepskyblue4") + 
  
  scale_fill_gradient(name = "Averted cases\nby percentage\n",
                      low = "white", high = "deepskyblue4",
                      breaks = c(min(combined_df$mean_percentage*100), 
                                 percentage_mark*100, 
                                 max(combined_df$mean_percentage*100)),
                      labels = c(paste0(round(min(combined_df$mean_percentage) * 100, 2),"%"),
                                 paste0(round(percentage_mark*100, 2),"%"),
                                 paste0(round(max(combined_df$mean_percentage) * 100, 2),"%"))) + 
  
  scale_shape_manual(name = " ",
                     values = setNames(18, paste0("≥ ", round(percentage_mark * 100, 2),"%"))) +  # Set shape to 18
  guides(
    fill = guide_colorbar(order = 1, barwidth=2, barheight = 20),  # Tile legend first
    color = guide_legend(order = 2),  # Segment legend second
    shape = guide_legend(order = 3)  # Point legend third
  ) +
  labs(#title = "mean Averted Cases by Percentage by Strategy",
    title = NULL,
    x = NULL,
    y = "Duration of protection (years)",
    fill = "CP") +
  coord_capped_cart(bottom='both', left='both') +
  theme(axis.title = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 30),
        plot.title = element_text(size = 30, 
                                  # face = text_font_weight, 
                                  hjust = 0),
        #strip.text = element_blank(),  # Remove subtitles
        legend.text = element_text(size = 30),  # Increase legend text size
        legend.title = element_text(size = 30),  # Increase legend title size
        # Hide panel borders and remove grid lines
        panel.background = element_blank(),
        # Facet wrap adjustment
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 30),
        panel.spacing = unit(50, "pt"),
        panel.border=element_blank(), axis.line=element_line(),
        plot.margin = margin(t = 10, r = 10, b = 20, l = 10))


plot2 <- ggplot(combined_df, aes(x = eff, y = dur, fill = mean_perdose)) +
  geom_tile() + # must before add on point
  geom_point(data = perdose_superior, aes(x = eff, y = dur, shape =  paste0("≥ ", round(perdose_mark, 2))), size = 10,  color = "mediumturquoise") +
  facet_rep_wrap(~strategy, ncol = 6, scales = "fixed",
                 # psuedo x axis label
                 labeller=as_labeller(c("VbE"="Efficacy (%)",
                                        "VoD(H)"="Efficacy (%)",
                                        "VoD"="Efficacy (%)",
                                        "VoA"="Efficacy (%)",
                                        "VoA(H)"="Efficacy (%)",
                                        "VaR"="Efficacy (%)")),
                 strip.position = "bottom") + 
  scale_x_continuous(limit = c(0, 1.05), breaks = seq (0.1,1,0.3) )+
  scale_fill_gradient(name = "Averted cases\nper dose\n",
                      low = "white", high = "indianred4",
                      breaks = c(min(combined_df$mean_perdose), 
                                 perdose_mark, 
                                 max(combined_df$mean_perdose)),
                      labels = c(paste0(min(combined_df$mean_perdose)),
                                 paste0(round(perdose_mark, 2)),
                                 paste0(round(max(combined_df$mean_perdose), 2))))+ 
  scale_shape_manual(name = " ",
                     values = setNames(18, paste0("≥ ", round(perdose_mark, 2)))) +  # Set shape to 18
 
  guides(
    fill = guide_colorbar(order = 1, barwidth=2, barheight = 20),  # Tile legend first
    color = guide_legend(order = 2),  # Segment legend second
    shape = guide_legend(order = 3)  # Point legend third
  ) +
  labs(#title = "mean Averted Cases by Strategy",
    title = NULL,
    x = NULL,
    y = "Duration of protection (years)",
    fill = "CD") +
  coord_capped_cart(bottom='both', left='both') +
  theme(axis.title = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        plot.title = element_text(size = 30, 
                                  # face = text_font_weight, 
                                  hjust = 0),
        legend.text = element_text(size = 30),  # Increase legend text size
        legend.title = element_text(size = 30),  # Increase legend title size
        # Hide panel borders and remove grid lines
        panel.background = element_blank(),
        # Facet wrap adjustment
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 30), # pseudo x axis title
        panel.spacing = unit(50, "pt"),
        panel.border=element_blank(), axis.line=element_line(),
        plot.margin = margin(t = 10, r = 10, b = 20, l = 10))



# Combine plots into a single big plot
big_plot <- plot1 + plot2 + plot_layout(ncol=1)


# Save the plot
#ggsave(file = paste0(path,"combined_plot.png"), big_plot, width = 48, height = 16, dpi = 320)
# Save the combined plot to a PDF file

ggsave(file = paste0(path,"combined_plot.pdf"), big_plot, width = 48, height = 16)
