# This script run gonorrhoea transmission projection given collected parameter samples

library(dplyr)
source("script_projection_support.R")
source("script_calibration_support.R")

version_check("gonovax", "0.4.16")
params <- readRDS("calibration/vaccine/parameters.rds")


# Title
if (is.null(strategy)) {
  title <- paste0("VbE")
} else if (strategy=="VoA"){
  title <- paste0("VoA")
} else if (strategy=="VoA(H)"){
  title <- paste0("VoA(H)")
} else if (strategy=="VoD+VoA(H)"){
  title <- paste0("VaR")
} else if (strategy=="VoS"){
  title <- paste0("VoS")
} else if (strategy=="VoD"){
  title <- paste0("VoD")
} else if (strategy=="VoD(H)"){
  title <- paste0("VoD(H)")
}
# short_run<-TRUE
n <- ifelse(short_run, 1000, nrow(params$health_economic))
idx <- seq_len(n)

if (fix_par_t) {
  params$gono <- params$gono$constant_risk[idx]
} else {
  params$gono <- params$gono$increasing_risk[idx]
}

params$health_economic <- params$health_economic[idx, ]
params$pathway_cost <- params$pathway_cost$incl_toc[idx, ]
params$vaccine <- params$vaccine[idx, ]


# run model from initialisation at end of 2004 to end of 2024, and extract
# starting conditions for use in the vaccination model
fit_years <- seq(2004, 2024)

novax_2024 <- lapply(params$gono, gonovax::run,
                     tt = gonovax_year_eliott(fit_years),
                     demographic_params = demographic_params_sg)
novax_2024_flows <- gonovax::extract_flows(novax_2024)


# run model with no vaccination from 2024 -> 2046 to extract a baseline
forecast_years <- seq(2024, 2046)
tt <- gonovax_year_eliott(forecast_years)

novax_2044 <- Map(gonovax::run,
                  gono_params = params$gono,

                  init_params = lapply(novax_2024, gonovax::restart_params),
                  MoreArgs = list(tt = tt,
                                  demographic_params = demographic_params_sg))

novax_2044_flows <- gonovax::extract_flows(novax_2044)



control <- run_grid_control(short_run, title, revax = TRUE) # define Vbe percentage
baseline <- rep(list(novax_2044_flows), length(control$eff) * length(control$dur))


# prepare for parallel
future::plan(future::multisession, workers = 8)

# run grid for all eff / dur combinations
if (vaccine_algorithm == "xvwv") {
  params$init <- lapply(novax_2024, gonovax::restart_params, n_vax = 3) # how many stratums
  # XVWV
  intervention <- run_grid_xvwv (gono_params = params$gono,
                                 init_params = params$init,
                                 cost_params = params$pathway_cost,
                                 uptake_second_dose = 1,
                                 model = control$model,
                                 eff = control$eff,
                                 dur = control$dur,
                                 vbe = control$vbe,                                    
                                 baseline = baseline,
                                 disc_rate = params$data$discount_rate,
                                 strategy = strategy,
                                 uptake_total = uptake_total,
                                 full_output = TRUE)
}


# Create the key using the formatted dur
# Check if control$dur is less than 10
 if (control$dur < 10) {
   # Ensure it is a two-digit number with leading zero
   dur_formatted <- sprintf("%05.2f", control$dur)
 } else {
   # If control$dur is 10 or greater, no need for formatting
   dur_formatted <- control$dur
 }
  
key <- paste0("eff", format(control$eff, nsmall = 2), "_dur", dur_formatted)
  


