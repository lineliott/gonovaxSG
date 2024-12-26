# this script runs through the calibration procedure of incidence over years to gonorrhoea transmission model in singapore setting.
library('ggplot2')
library('dplyr')
source("script_calibration_plot_support.R")
source("script_calibration_support.R")
set.seed(111)

path <- "calibration_sgmsm_try/"
dir.create(path, FALSE, TRUE)
dir.create(paste0(path, "figs/"), FALSE, TRUE)
dir.create(paste0(path, "outputs/"), FALSE, TRUE)



version_check("gonovax", "0.4.16") 
# read data
stream <- "msm"
if (stream=="msm") {
  sg_data <- read_csv("MSM.csv")
} else if (stream == "male"){
  sg_data <- read_csv("male.csv")
} else if (stream == "msm_uk"){
  sg_data <- read_csv("male.csv")
  sg_data$diagnosed <- sg_data$diagnosed*0.7
}


sg_data$diagnosed <- round(as.integer(sg_data$diagnosed))
sg_data <- sg_data[1:15,] # till 2018


## set run_parameters
short_run <- TRUE

#kernel_scaling <- 0.5
#kernel_scaling <- 0.001
kernel_scaling <- 0.05


control <- mcmc_control(short_run)

parameters_info <- read_csv("parameters_info.csv")
parameters_prior <- read_csv("parameters_prior.csv")
parameters_proposal <- read_csv("parameters_proposal.csv", row.names = 1)

# delete k_grasp, not allow negative beta_2018,
{
  parameters_info <- subset(parameters_info, !(name == "k_grasp"))
  parameters_info <- subset(parameters_info, !(name == "phi_eta"))
  parameters_info[parameters_info["name"]=="beta2009",][["name"]] <- "beta_2004"
  parameters_prior[parameters_prior["name"]=="beta2009",][["name"]] <- "beta_2004"
  parameters_info[parameters_info["name"]=="phi_beta",][["name"]] <- "beta_2018"
  parameters_prior[parameters_prior["name"]=="beta2020",][["name"]] <- "beta_2018"
  parameters_info[parameters_info["name"]=="k_gumcad",][["name"]] <- "k_surveillance"
  parameters_prior[parameters_prior["name"]=="k_gumcad",][["name"]] <- "k_surveillance"
  
  parameters_prior <- subset(parameters_prior, !(name == "k_grasp"))
  parameters_proposal <- parameters_proposal[!(rownames(parameters_proposal) == "k_grasp"), ]
  parameters_proposal$k_grasp <- NULL
  # Update row and colnames
  rownames(parameters_proposal) <- gsub("beta2009", "beta_2004", rownames(parameters_proposal))
  colnames(parameters_proposal) <- gsub("beta2009", "beta_2004", colnames(parameters_proposal))
  rownames(parameters_proposal) <- gsub("k_gumcad", "k_surveillance", rownames(parameters_proposal))
  colnames(parameters_proposal) <- gsub("k_gumcad", "k_surveillance", colnames(parameters_proposal))
  rownames(parameters_proposal) <- gsub("phi_beta", "beta_2018", rownames(parameters_proposal))
  colnames(parameters_proposal) <- gsub("phi_beta", "beta_2018", colnames(parameters_proposal))
  
  parameters_info$min[parameters_info$name == "beta_2004"] <- 0
  parameters_info$max[parameters_info$name == "beta_2004"] <- 0.8
  parameters_info$initial[parameters_info$name == "beta_2004"] <- 0.8

  parameters_info$min[parameters_info$name == "beta_2018"] <- 0
  parameters_info$max[parameters_info$name == "beta_2018"] <- 0.8
  parameters_info$initial[parameters_info$name == "beta_2018"] <- 0.8

  
  # parameters_info$max[parameters_info$name == "prev_Ash"] <- 0.1#msm
  # parameters_info$initial[parameters_info$name == "prev_Ash"] <- 0.1#msm
  parameters_info$max[parameters_info$name == "prev_Ash"] <- 0.086#msm
  parameters_info$initial[parameters_info$name == "prev_Ash"] <- 0.086#msm
  
  # parameters_info$max[parameters_info$name == "prev_Ash"] <- 0.0095#male
  # parameters_info$initial[parameters_info$name == "prev_Ash"] <- 0.0095 #male
  # parameters_info$max[parameters_info$name == "prev_Ash"] <- 0.05 # msm_uk
  # parameters_info$initial[parameters_info$name == "prev_Ash"] <- 0.05 #msm_uk

  # parameters_info$max[parameters_info$name == "prev_Asl"] <- 0.01#msm
  # parameters_info$initial[parameters_info$name == "prev_Asl"] <- 0.01#msm
  parameters_info$max[parameters_info$name == "prev_Asl"] <- 0.015#msm
  parameters_info$initial[parameters_info$name == "prev_Asl"] <- 0.015#msm
  
  # parameters_info$max[parameters_info$name == "prev_Asl"] <- 0.016 # male
  # parameters_info$initial[parameters_info$name == "prev_Asl"] <- 0.016 #male
  # parameters_info$max[parameters_info$name == "prev_Asl"] <- 0.013 # msm_uk
  # parameters_info$initial[parameters_info$name == "prev_Asl"] <- 0.013 #msm_uk
  
  parameters_info$min[parameters_info$name == "eta_h"] <- 0
  parameters_info$max[parameters_info$name == "eta_h"] <- 4
  parameters_info$min[parameters_info$name == "omega"] <- 0
  parameters_info$max[parameters_info$name == "omega"] <- 1
  parameters_info$min[parameters_info$name == "nu"] <- 3
  parameters_info$initial[parameters_info$name == "nu"] <- 3
  # should be 3-8 but nosense
  # parameters_info$min[parameters_info$name == "nu"] <- 3
  # parameters_info$initial[parameters_info$name == "nu"] <- 3
  parameters_info$max[parameters_info$name == "nu"] <- 6
  parameters_info$min[parameters_info$name == "sigma"] <- 26
  # parameters_info$min[parameters_info$name == "sigma"] <- 60
  parameters_info$max[parameters_info$name == "sigma"] <- 180
  
  parameters_info$min[parameters_info$name == "mu"] <- 52
  parameters_info$max[parameters_info$name == "mu"] <- 365
  #parameters_info$min[parameters_info$name == "psi"] <- 0
  parameters_info$min[parameters_info$name == "psi"] <- 0.1
  parameters_info$max[parameters_info$name == "psi"] <- 0.2
  # parameters_info$max[parameters_info$name == "psi"] <- 0.4
  parameters_info$min[parameters_info$name == "rho"] <- 12
  parameters_info$max[parameters_info$name == "rho"] <- 90
  # parameters_info$min[parameters_info$name == "rho"] <- 26
  # parameters_info$max[parameters_info$name == "rho"] <- 52
  # parameters_info$initial[parameters_info$name == "rho"] <- 52
  
}

# Example usage with the provided posterior estimate
{
  uk_rho <- c(54.0, 43.6, 66.4)
  uk_mu <- c(218, 92.9, 521) 
  uk_nu <- c(3.08, 1.60, 6.03) 
  uk_psi <- c(0.150, 0.0737, 0.238)
  uk_sigma <- c(99.9, 56.1, 176)
  uk_omega <- c(0.475, 0.218, 0.876)
  uk_eta_h <- c(0.184, 0.107, 0.3)
  
  sg_rho <- transform_to_normal_distribution(uk_rho)
  sg_mu <- transform_to_normal_distribution(uk_mu)
  sg_nu <- transform_to_normal_distribution(uk_nu)
  sg_psi <- transform_to_normal_distribution(uk_psi)
  sg_sigma <- transform_to_normal_distribution(uk_sigma)
  sg_omega <- transform_to_normal_distribution(uk_omega)
  sg_eta_h <- transform_to_normal_distribution(uk_eta_h)
}

# priors adjustment 
{
  # Create a new row with "name" column value as "psi"
  new_row <- data.frame(name = "psi", type = "gaussain", beta_shape1 = NA, beta_shape2 = NA, gamma_shape = NA, gamma_scale = NA, lnorm_meanlog = NA, lnorm_sdlog = NA, as_duration = FALSE, lnorm_mean = NA,lnorm_q2.5 = NA, lnorm_q97.5= NA, gamma_mean=NA,gamma_q2.5=NA, gamma_q97.5=NA, beta_mean=NA, beta_q2.5=NA, beta_q97.5=NA)
  parameters_prior <- rbind(parameters_prior, new_row)
  new_row <- data.frame(name = "epsilon", type = "null", beta_shape1 = NA, beta_shape2 = NA, gamma_shape = NA, gamma_scale = NA, lnorm_meanlog = NA, lnorm_sdlog = NA, as_duration = FALSE, lnorm_mean = NA,lnorm_q2.5 = NA, lnorm_q97.5= NA, gamma_mean=NA,gamma_q2.5=NA, gamma_q97.5=NA, beta_mean=NA, beta_q2.5=NA, beta_q97.5=NA)
  parameters_prior <- rbind(parameters_prior, new_row)
  new_row <- data.frame(name = "prev_Ash", type = "null", beta_shape1 = NA, beta_shape2 = NA, gamma_shape = NA, gamma_scale = NA, lnorm_meanlog = NA, lnorm_sdlog = NA, as_duration = FALSE, lnorm_mean = NA,lnorm_q2.5 = NA, lnorm_q97.5= NA, gamma_mean=NA,gamma_q2.5=NA, gamma_q97.5=NA, beta_mean=NA, beta_q2.5=NA, beta_q97.5=NA)
  parameters_prior <- rbind(parameters_prior, new_row)
  new_row <- data.frame(name = "prev_Asl", type = "null", beta_shape1 = NA, beta_shape2 = NA, gamma_shape = NA, gamma_scale = NA, lnorm_meanlog = NA, lnorm_sdlog = NA, as_duration = FALSE, lnorm_mean = NA,lnorm_q2.5 = NA, lnorm_q97.5= NA, gamma_mean=NA,gamma_q2.5=NA, gamma_q97.5=NA, beta_mean=NA, beta_q2.5=NA, beta_q97.5=NA)
  parameters_prior <- rbind(parameters_prior, new_row)
  # Append empty columns "guassian_mean" and "guassian_sd" at index 3 and 4
  new_columns <- data.frame(gaussian_mean = NA, gaussian_sd = NA) 
  parameters_prior <- cbind(parameters_prior[, 1:2], new_columns, parameters_prior[, 3:ncol(parameters_prior)])
  
  # assgian UK posterior as prior
  parameters_prior[parameters_prior$name == "beta_2018",3:20] <- NA 
  parameters_prior[parameters_prior$name == "beta_2018", "type"] <- "null"
  parameters_prior[parameters_prior$name == "beta_2018", "as_duration"] <- FALSE
  
  
  
  parameters_prior[parameters_prior$name == "nu", 3:20] <- NA 
  parameters_prior[parameters_prior$name == "nu", "type"] <- "gaussain"
  parameters_prior[parameters_prior$name == "nu", "gaussian_mean"] <- sg_nu$gaussian_mean
  parameters_prior[parameters_prior$name == "nu", "gaussian_sd"] <- sg_nu$gaussian_sd
  parameters_prior[parameters_prior$name == "nu", "as_duration"] <- TRUE
  
  parameters_prior[parameters_prior$name == "psi",3:20] <- NA 
  parameters_prior[parameters_prior$name == "psi", "type"] <- 'gaussain'
  parameters_prior[parameters_prior$name == "psi", "gaussian_mean"] <- sg_psi$gaussian_mean
  parameters_prior[parameters_prior$name == "psi", "gaussian_sd"] <- sg_psi$gaussian_sd
  parameters_prior[parameters_prior$name == "psi", "as_duration"] <- FALSE
  
  parameters_prior[parameters_prior$name == "omega", 3:20] <- NA 
  parameters_prior[parameters_prior$name == "omega", "type"] <- 'gaussain'
  parameters_prior[parameters_prior$name == "omega", "gaussian_mean"] <- sg_omega$gaussian_mean
  parameters_prior[parameters_prior$name == "omega", "gaussian_sd"] <- sg_omega$gaussian_sd
  parameters_prior[parameters_prior$name == "omega", "as_duration"] <- FALSE
  
  parameters_prior[parameters_prior$name == "mu", 3:20] <- NA 
  parameters_prior[parameters_prior$name == "mu", "type"] <- 'gaussain'
  parameters_prior[parameters_prior$name == "mu", "gaussian_mean"] <- sg_mu$gaussian_mean
  parameters_prior[parameters_prior$name == "mu", "gaussian_sd"] <- sg_mu$gaussian_sd
  parameters_prior[parameters_prior$name == "mu", "as_duration"] <- TRUE
  
  parameters_prior[parameters_prior$name == "rho", 3:20] <- NA 
  parameters_prior[parameters_prior$name == "rho", "type"] <- 'gaussain'
  parameters_prior[parameters_prior$name == "rho", "gaussian_mean"] <- sg_rho$gaussian_mean
  parameters_prior[parameters_prior$name == "rho", "gaussian_sd"] <- sg_rho$gaussian_sd
  parameters_prior[parameters_prior$name == "rho", "as_duration"] <- TRUE
  
  parameters_prior[parameters_prior$name == "sigma", 3:20] <- NA 
  parameters_prior[parameters_prior$name == "sigma", "type"] <- 'gaussain'
  parameters_prior[parameters_prior$name == "sigma", "gaussian_mean"] <- sg_sigma$gaussian_mean
  parameters_prior[parameters_prior$name == "sigma", "gaussian_sd"] <- sg_sigma$gaussian_sd
  parameters_prior[parameters_prior$name == "sigma", "as_duration"] <- TRUE
  
  parameters_prior[parameters_prior$name == "eta_h", 3:20] <- NA 
  parameters_prior[parameters_prior$name == "eta_h", "type"] <- 'gaussain'
  parameters_prior[parameters_prior$name == "eta_h", "gaussian_mean"] <- sg_eta_h$gaussian_mean
  parameters_prior[parameters_prior$name == "eta_h", "gaussian_sd"] <- sg_eta_h$gaussian_sd
  parameters_prior[parameters_prior$name == "eta_h", "as_duration"] <- FALSE
  
  
  
  
  order <- c("beta_2004","beta_2018","prev_Asl","prev_Ash","epsilon",
             "eta_h","omega","nu","sigma",
             "mu","psi","rho","k_surveillance")
  parameters_info <- parameters_info %>%
    mutate(name = factor(name, levels = order)) %>%
    arrange(name)
  parameters_info$name <- as.character(parameters_info$name)
  parameters_prior <- parameters_prior %>%
    mutate(name = factor(name, levels = order)) %>%
    arrange(name) %>% 
    filter(!is.na(name)) 
  parameters_prior<-as.data.frame(parameters_prior)
}

if (kernel_scaling == 0) {
  ## Used when tuning kernels only
  kernel_scaling <- 2.38 ^ 2 / nrow(parameters_info)
}
message("  - scaling proposal kernel by ", signif(kernel_scaling, 3))
parameters_proposal <- parameters_proposal * kernel_scaling

message("Preparing")
priors <- prepare_priors(parameters_prior)
pars <- prepare_parameters(parameters_info, priors, parameters_proposal)

###################################################################################
message("Running chains - this will take a while! ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

future::plan(future::multicore, workers = control$n_chains)

n_steps = rep(control$n_mcmc, control$n_chains)
mcmc_raw <- 
  furrr::future_pmap(.l = list(n_steps = rep(control$n_mcmc, control$n_chains)),
                     .f = mcmc_eliott, 
                     pars = pars, 
                     compare = compare_eliott, 
                     #n_chains = control$n_chains, 
                     #progress = TRUE,
                     .options = furrr::furrr_options(seed = TRUE))

message("mcmc finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
mcmc_results <- mcmc_raw
mcmc_metrics <- calculate_mcmc_metrics(mcmc_results)
mcmc_results <- gonovax::mcmc_combine(samples = mcmc_results)
mcmc_sample <- gonovax::mcmc_sample(mcmc_results, control$n_sample, control$burnin)
mcmc_sample$transform <- transform_eliott
mcmc_results$metrics <- mcmc_metrics

message("Writing outputs")


saveRDS(mcmc_raw, paste0(path, "outputs/mcmc_raw_0.rds"))
saveRDS(mcmc_results, paste0(path, "outputs/mcmc_results_0.rds"))
saveRDS(mcmc_sample, paste0(path, "outputs/mcmc_sample_0.rds"))
write.csv(cov(mcmc_sample$pars), paste0(path, "outputs/parameters_proposal_0.csv"))
parameters_info <- calculate_parameters_info(mcmc_results, parameters_info)
write.csv(parameters_info, paste0(path, "outputs/parameters_info_0.csv", row.names = FALSE))


message("Creating plots")

write_png(paste0(path,"figs/mcmc_traces_0.png"), asp = 1, pointsize = 8, {
  par(mfcol = c(5, 4), mar = c(3, 3, 3, 0))
  plot_traces(mcmc_results, mcmc_metrics)
})

write_png(paste0(path, "figs/mcmc_posteriors_0.png"), asp = 3 / 5, {
  par(mfrow = c(3, 5), mar = c(3, 3, 2, 1))
  plot_posteriors(mcmc_sample, priors)
})

p <- boxplot_fits_eliott(data, mcmc_sample)
ggsave(paste0(path, "figs/mcmc_fit_0_box.png"), plot = p, width = 12, height = 9)

# calibration trace
plot_mcmc_traces_ggplot(mcmc_results, 
                        mcmc_metrics, 
                        paste0(path, "figs/mcmc_traces_0"))

########################################################################
# tuned run
message("final tuned run")
proposal_matrix <- cov(mcmc_sample$pars)
saveRDS(parameters_info, paste0(path, "outputs/parameters_info_tuned.rds"))
saveRDS(parameters_prior, paste0(path, "outputs/parameters_prior_tuned.rds"))
saveRDS(parameters_proposal, paste0(path, "outputs/parameters_proposal_tuned.rds"))
message("  - scaling proposal kernel by ", signif(kernel_scaling, 3))
parameters_proposal <- proposal_matrix * kernel_scaling

for (j in 1 : length(colnames(mcmc_sample[["pars"]]))){
  object <- colnames(mcmc_sample[['pars']])[j]
  initial <- median(mcmc_sample[['pars']][,j])
  parameters_info[parameters_info["name"] == object,][['initial']] <- initial
}

pars <- prepare_parameters(parameters_info, priors, parameters_proposal)

#compare_trans <- function(pars) compare_eliott(pars, transform_eliott)
short_run <- FALSE
control <- mcmc_control(short_run)
message("Running chains (final tuned) - this will take a while!", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

future::plan(future::multicore, workers = control$n_chains)

mcmc_raw <- 
  furrr::future_pmap(.l = list(n_steps = rep(control$n_mcmc, control$n_chains)),
                     
                     .f = mcmc_eliott, pars = pars, compare = compare_eliott, 
                     .progress = TRUE,
                     .options = furrr::furrr_options(seed = TRUE))
message("tuned mcmc finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

mcmc_results <- mcmc_raw
mcmc_metrics <- calculate_mcmc_metrics(mcmc_results)
mcmc_results <- gonovax::mcmc_combine(samples = mcmc_results)
mcmc_sample <- gonovax::mcmc_sample(mcmc_results, control$n_sample, control$burnin)
mcmc_sample$transform <- transform_eliott
mcmc_results$metrics <- mcmc_metrics

message("Writing tuned outputs")

saveRDS(mcmc_raw, paste0(path, "outputs/mcmc_raw_tuned.rds"))
saveRDS(mcmc_results, paste0(path, "outputs/mcmc_results_tuned.rds"))
saveRDS(mcmc_sample, paste0(path, "outputs/mcmc_sample_tuned.rds"))
write.csv(cov(mcmc_sample$pars), paste0(path, "outputs/parameters_proposal_tuned.csv"))
parameters_info <- calculate_parameters_info(mcmc_results, parameters_info)
write.csv(parameters_info, paste0(path, "outputs/parameters_info_tuned.csv", row.names = FALSE))
######################################################################################
message("Creating tuned plots")

mcmc_raw <- readRDS(paste0(path, "outputs/mcmc_raw_tuned.rds"))
mcmc_results <- readRDS(paste0(path, "outputs/mcmc_results_tuned.rds"))
mcmc_sample <- readRDS(paste0(path, "outputs/mcmc_sample_tuned.rds"))
mcmc_metrics <- mcmc_results$metrics

write_png(paste0(path, "figs/mcmc_posteriors_tuned.png"), asp = 3 / 5, {
  par(mfrow = c(3, 5), mar = c(3, 3, 2, 1))
  plot_posteriors(mcmc_sample, priors)
})

p <- boxplot_fits_eliott(data, mcmc_sample)
ggsave(paste0(path,"figs/mcmc_fit_tuned_box.pdf"), plot = p, width = 24, height = 16)

# calibration trace
plot_mcmc_traces_ggplot(mcmc_results, 
                        mcmc_metrics, 
                        paste0(path, "figs/mcmc_traces_tuned"))

