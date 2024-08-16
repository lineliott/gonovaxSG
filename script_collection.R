# source("global_util.R")
source("support.R")
#source("global_run_gonovax.R")
#source("global_mcmc_support.R")
source("script_calibration_support.R")
source("script_calibration_plot.R")
set.seed(111)

calculate_pathway_cost <- function(unit_cost, p_toc) {
  
  cost_toc <- p_toc * unit_cost$toc
  
  cost_a <- unit_cost$test + unit_cost$treat + cost_toc
  cost_s <- unit_cost$test + unit_cost$treat * (1 - unit_cost$discount_treat_sympt) +
    unit_cost$manage_symptoms + cost_toc
  
  data.frame(unit_cost_manage_symptomatic = cost_s,
             unit_cost_manage_asymptomatic = cost_a,
             unit_cost_screen_uninfected = unit_cost$test,
             qaly_loss_per_diag_s = unit_cost$qaly_loss_per_diag_s)
}


dir.create("calibration/vaccine", FALSE, TRUE)
# control parameters
control <- mcmc_control(short_run = FALSE)

#mcmc_results <- readRDS("MALE/meeting5/outputs/mcmc_results_1.rds")
#mcmc_results <- readRDS("calibration_01/outputs/mcmc_results_2.rds")
mcmc_results <- readRDS("calibration/outputs/mcmc_results_tuned.rds")

n_par <- 1000 #sampling quantity
sample <- gonovax::mcmc_sample(mcmc_results, n_par, control$burnin)
#sample$transform <- gonovax::transform
sample$transform <- transform_eliott

# uptake parameters
# HPV parameters reported by PHE, PW email 2021-01-05
# MeNZb params from Petousis-Harris
data <-
  list(uptake = list(offered  = list(first_dose = 76033, second_dose = 14612),
                     accepted = list(first_dose = 32562, second_dose = 11267)),
       menzb_efficacy = list(mean = 0.31, l = 0.21, u = 0.39))


# assume each dose uptake ~ Beta(accepted, offered - accepted)
vaccine_params <-
  with(data$uptake, Map(sample_uptake, offered, accepted, n_par))

# calculate overall uptake
vaccine_params$total <- vaccine_params$first_dose * vaccine_params$second_dose
names(vaccine_params) <- sprintf("uptake_%s", names(vaccine_params))
mean_ci(vaccine_params$uptake_total)

menzb_efficacy_pars <- with(data$menzb_efficacy, fit_beta(mean, l, u))
vaccine_params$menzb_efficacy <- rbeta(n_par, menzb_efficacy_pars$alpha,
                                       menzb_efficacy_pars$beta)
mean_ci(vaccine_params$menzb_efficacy)


## use 3.5% discount rate
data$discount_rate <- 0.035
data$bex_price <- 85

# test of cure parameters
# test of cure data from GRASP 2019 report (published 2020)
data$toc <- list(offered = 971, accepted = 552)
p_toc <- with(data$toc, sample_uptake(offered, accepted, n_par))
mean_ci(p_toc)

# cost-effectiveness
# allow for a standard error of mean ± 20% 
se <- 0.2 
data$unit_cost <- list(test = 88.35,
                       manage_symptoms = 21.72,
                       treat = 79.82,
                       toc = 44.28)
unit_cost <- Map(sample_gamma, data$unit_cost, n_par, se)

# check generated samples
testthat::expect_equivalent(sapply(unit_cost, mean) / unlist(data$unit_cost),
                            rep(1, length(unit_cost)), tolerance = 0.05)
testthat::expect_equivalent(sapply(unit_cost, sd) / unlist(data$unit_cost) / se,
                            rep(1, length(unit_cost)), tolerance = 0.05)

# 14% discount when treating symptomatic patients (i.e at first visit)
data$unit_cost$discount_treat_sympt <- 0.14
unit_cost$discount_treat_sympt <- data$unit_cost$discount_treat_sympt

#' Disutility of gonorrhoea symptoms was obtained from literature with
#' uncertainty represented using a Pert distribution with range ±20%;

data$qaly_loss_pa <- 0.16
qaly_loss_pa <- freedom::rpert(n = n_par,
                      x.mode = data$qaly_loss_p,
                      x.min = data$qaly_loss_pa * (1 - se),
                      x.max = data$qaly_loss_pa * (1 + se))

mean_ci(qaly_loss_pa)

#' QALY loss was determined by the assumed average duration of symptoms, the
#' time until obtaining care (1/μ) plus half the duration of treatment (1/2ρ)

duration <- 1 / sample$pars[ , "mu"] + 0.5 / sample$pars[, "rho"]
qaly_loss_per_diag_s <- qaly_loss_pa * duration


# combine parameters and save

cost_eff_params <- data.frame(p_toc, unit_cost, qaly_loss_pa, qaly_loss_per_diag_s)

incl_toc <- calculate_pathway_cost(cost_eff_params, p_toc = p_toc)
excl_toc <- calculate_pathway_cost(cost_eff_params, p_toc = 0)
pathway_cost <- list(incl_toc = incl_toc,
                     excl_toc = excl_toc)

names(incl_toc) <- paste0(names(incl_toc), "_incl_toc")
names(excl_toc) <- paste0(names(excl_toc), "_excl_toc")

# plot parameters
write_png("calibration/vaccine/cost_eff_params.png", {
 par(mfrow = c(5, 3))
 plot_param(cost_eff_params)
 plot.new()
 plot_param(excl_toc[-4], xlim = c(0, 300))

 plot_param(incl_toc[-4], xlim = c(0, 300))
})

write_png("calibration/vaccine/vaccine_params.png", {
  par(mfrow = c(2, 2))
  plot_param(vaccine_params)
})


## extract gono_params in both constant risk and increasing risk versions
gono_params <- list(constant_risk = apply(sample$pars, 1, sample$transform,
                                          fix_par_t = TRUE),
                    increasing_risk = apply(sample$pars, 1, sample$transform,
                                            fix_par_t = FALSE))

## extract pathway_cost in both incl_toc and excl_toc versions


parameters <- list(gono = gono_params,
                  health_economic = cost_eff_params,
                  pathway_cost = pathway_cost,
                  vaccine = as.data.frame(vaccine_params),
                  data = data)
{
  # Extract the "beta_2004" variable from the "constant_risk" column for all samples
  beta_2004_values <- unlist(lapply(gono_params[["constant_risk"]], function(x) x$beta_2004))
  
  # Plot the histogram
  hist(beta_2004_values, main = "Histogram of beta_2004", xlab = "Beta Values", col = "skyblue", border = "black")
  mean(beta_2004_values)
  
  # Extract the "beta_2018" variable from the "constant_risk" column for all samples
  beta_2018_values <- unlist(lapply(gono_params[["constant_risk"]], function(x) x$beta_2018))
  
  # Plot the histogram
  hist(beta_2018_values, main = "Histogram of beta_2018", xlab = "Phi Beta Values", col = "lightgreen", border = "black")
  mean(beta_2018_values)
}



saveRDS(parameters, "calibration/vaccine/parameters.rds")
saveRDS(sample, "calibration/vaccine/mcmc_sample.rds")
