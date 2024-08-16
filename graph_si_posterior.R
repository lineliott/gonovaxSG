# Define the order of metrics
metric_order <- c("beta_2004", "beta_2018", "prev_Asl", "prev_Ash", "epsilon", "eta_h", 
                  "omega", "nu", "sigma", "mu", "psi", "rho", "k_surveillance")

# Initialize an empty data frame to hold the combined results
combined_table <- data.frame(metric = metric_order,
                             stringsAsFactors = FALSE)

# Stream loop
stream <- c("sgmsm", "sgmale", "ukmsm")
for (i in seq_along(stream)) {
  
  if (stream[i] == "sgmsm") {
    mcmc_sample <- readRDS("calibration_sgmsm/calibration/vaccine/mcmc_sample.rds")
    title <- "Main"
  } else if (stream[i] == "sgmale") {
    mcmc_sample <- readRDS("calibration_sgmale/calibration/vaccine/mcmc_sample.rds")
    title <- "Upper"
  } else if (stream[i] == "ukmsm") {
    mcmc_sample <- readRDS("calibration_ukmsm/calibration/vaccine/mcmc_sample.rds")
    title <- "Lower"
  }
  
  par_names <- colnames(mcmc_sample$pars)
  ylab <- par_names
  ylab[1:4] <- paste(ylab[1:4], "(%)")
  fac <- ifelse(ylab == par_names, 1, 100)
  
  # Initialize table for current stream
  table <- data.frame(metric = par_names, value = rep("", length(par_names)), stringsAsFactors = FALSE)
  
  for (j in seq_along(par_names)) {
    qs <- quantile(mcmc_sample$pars[, j], c(0.5, 0.025, 0.975)) * fac[j]
    mean <- mean(mcmc_sample$pars[, j]) * fac[j]
    title_text <- sprintf("%.2f (%.2f, %.2f)", mean, qs[2], qs[3])
    table[table$metric == par_names[j], "value"] <- title_text
  }
  
  # Merge into combined_table based on metric_order
  combined_table <- merge(combined_table, table, by = "metric", all = TRUE)
  colnames(combined_table)[which(names(combined_table) %in% c("value.x", "value.y"))] <- c("Main", "Upper")
}

# Reorder combined_table according to metric_order
combined_table <- combined_table[match(metric_order, combined_table$metric), ]

# Save combined_table as CSV
write.csv(combined_table, file = "SI/posterior.csv", row.names = FALSE)

# Confirmation message
cat("Combined table saved successfully as SI/posterior.csv.\n")
