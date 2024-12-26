# This script produces figure S2 (mcmc convergency check)
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

options(scipen = 999)

# Function to plot MCMC traces and metrics
plot_mcmc_traces <- function(mcmc_results, title, plots) {
  
  n_chains <- max(mcmc_results$chain)
  col <- hcl.colors(n_chains, "Geyser", rev = TRUE)
  
  par_names <- colnames(mcmc_results$pars)
  ylab <- par_names

  # Generate trace plots for parameters
  for (i in seq_along(par_names)) {

    ess <- mcmc_results$metrics$ess[i]
    gr <- mcmc_results$metrics$gr$psrf[i, 2]
    # title_text <- sprintf("%.2f (%.2f, %.2f)\nESS = %.f, GR = %.2f",
    #                       mean, qs[2], qs[3], ess, gr)
    title_text <- sprintf("ESS = %.f, GR = %.2f",
                           ess, gr)
    
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
                   aes(x = iteration / 100000, y = value, color = factor(chain))) +
      geom_line() +
      labs(x = "Iteration per 100,000", y = ylab[i], color = title) +  # Use 'title' instead of 'Chain'
      ggtitle(title_text) +
      scale_color_manual(values = col) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +  # Specify 3 breaks on the y-axis
      theme(axis.title = element_text(size = 30),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 30),
            plot.title = element_text(size = 30, hjust = 0.5),
            legend.text = element_text(size = 30),
            legend.title = element_text(size = 30),  # Adjust legend title size
            panel.background = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(size = 30),
            panel.spacing = unit(50, "pt"),
            panel.border = element_blank(),
            axis.line = element_line(),
            plot.margin = margin(t = 10, r = 10, b = 20, l = 10))
    
    plots[[length(plots) + 1]] <- plot
  }
  
  # Generate trace plot for log_posterior
  j <- "log_posterior"
  
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
  
  ar <- mcmc_results$metrics$acceptance_rate[[1]] * 100
  title_text <- sprintf("%s\nAcceptance Ratio = %.2f%%", j, ar) 
  plot <- ggplot(data = plot_data,
                 aes(x = iteration / 100000, y = value, color = factor(chain))) +
    geom_line() +
    labs(x = "Iteration per 100,000", y = j, color = title) +  # Use 'title' instead of 'Chain'
    scale_color_manual(values = col) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +  # Specify 3 breaks on the y-axis
    ggtitle(title_text) +
    theme(axis.title = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          plot.title = element_text(size = 30, hjust = 0.5),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30),  # Adjust legend title size
          panel.background = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 30),
          panel.spacing = unit(50, "pt"),
          panel.border = element_blank(),
          axis.line = element_line(),
          plot.margin = margin(t = 10, r = 10, b = 20, l = 10))
  
  plots[[length(plots) + 1]] <- plot
  
  # Return list of plots
  return(plots)
}

# Combine the plotting logic into a function
combine_plots <- function(plots) {
  # Adjust the theme for each plot based on its position in the grid
  ncol <- 7
  nrow <- ceiling(length(plots) / ncol)
  
  for (i in seq_along(plots)) {
    plots[[i]] <- plots[[i]] +
      theme(
        legend.position = if (i %% ncol == 0) "right" else "none",
        axis.text.x = if (ceiling(i / ncol) == nrow) element_text(size = 30) else element_blank(),
        axis.title.x = if (ceiling(i / ncol) == nrow) element_text(size = 30) else element_blank()
      )
  }
  return(plots)
}

# Stream loop
stream <- c("sgmsm", "sgmale", "ukmsm")

for (i in seq_along(stream)) {
  
  if (stream[i] == "sgmsm") {
    mcmc_results <- readRDS("calibration_sgmsm/calibration/outputs/mcmc_results_tuned.rds")
    title <- "Main"
  } else if (stream[i] == "sgmale") {
    mcmc_results <- readRDS("calibration_sgmale/calibration/outputs/mcmc_results_tuned.rds")
    title <- "Upper"
  } else if (stream[i] == "ukmsm") {
    mcmc_results <- readRDS("calibration_ukmsm/calibration/outputs/mcmc_results_tuned.rds")
    title <- "Lower"
  }
  
  plots <- list()
  # Plotting for current stream
  plots <- plot_mcmc_traces(mcmc_results, title, plots)
  
  
  # Arrange legend
  plots <- combine_plots(plots)
  
  # Save the combined plots to PNG
  png(paste0("SI/mcmc_",title,".png"), width = 64, height = 10, units = "in", res = 320)
  grid.arrange(grobs = plots, ncol = 7)
  dev.off()
  
  # Save the combined plots to PDF
  pdf(paste0("SI/mcmc_", title, ".pdf"), width = 64, height = 10)
  grid.arrange(grobs = plots, ncol = 7)
  dev.off()
  
  

}

