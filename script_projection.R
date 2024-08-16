# orderly::orderly_develop_start(parameters = list(short_run = TRUE))
library(dplyr)
# source("global_util.R")
# source("global_plot_util.R")
source("script_projection_support.R")
# source("util.R")
# source("utils_assert.R")
source("script_calibration_support.R")

version_check("gonovax", "0.4.16")
params <- readRDS("calibration/vaccine/parameters.rds")
#params_previous <- readRDS("parameters.rds")


# Title
if (is.null(strategy)) {
  # if (control$vbe == 0) (
  #   title <- paste0("novax"))
  # else if (control$vbe != 0) (
  #   title <- paste0("VbE"))
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


#fix_par_t<-TRUE # parameters are kept constant at their inferred value at the final data point 
#fix_par_t<-FALSE # continue after last data point
if (fix_par_t) {
  params$gono <- params$gono$constant_risk[idx]
} else {
  params$gono <- params$gono$increasing_risk[idx]
}

params$health_economic <- params$health_economic[idx, ]
params$pathway_cost <- params$pathway_cost$incl_toc[idx, ]
params$vaccine <- params$vaccine[idx, ]

# eta_try
# for (i in 1:100){
#   params[["gono"]][[i]][["eta_h_t"]]<-params[["gono"]][[i]][["eta_h_t"]]+0.5
# }

## run model from initialisation at end of 2009 to end of 2021, and extract
## starting conditions for use in the vaccination model
#fit_years <- seq(2009, 2021)
fit_years <- seq(2004, 2024)

novax_2024 <- lapply(params$gono, gonovax::run,
                     tt = gonovax_year_eliott(fit_years),
                     demographic_params = demographic_params_sg)
novax_2024_flows <- gonovax::extract_flows(novax_2024)

# eliott
# novax_2024_diff <- lapply(novax_2024, function(sublist) {
#   sublist <- diff(sublist$tot_treated)
#   return(sublist)
# })

## run model with no vaccination from 2022 -> 2042 to extract a baseline
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

#future::plan(future::multisession, workers = cores)
future::plan(future::multisession, workers = 8)

## run grid for all eff / dur combinations
# if (vaccine_algorithm == "xvwr") {
#   params$init <- lapply(novax_2024, gonovax::restart_params, n_vax = 4) # how many stratums
#   # XVWR
#   intervention <- run_grid_xvwr (gono_params = params$gono,
#                                  init_params = params$init,
#                                  cost_params = params$pathway_cost,
#                                  primary_uptake = 0.43,
#                                  booster_uptake = 0.77,
#                                 # uptake_second_dose = 1,
#                                  model = control$model,
#                                  eff = control$eff,
#                                  dur = control$dur,
#                                  dur_revax = control$dur_revax,
#                                  vbe = control$vbe,
#                                  baseline = baseline,
#                                  disc_rate = params$data$discount_rate,
#                                  strategy = strategy,
#                                  # uptake_total = uptake_total,
#                                  full_output = TRUE)
# }
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




if (short_run){
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
  
} else { 
  
  if (graph == "heatmap") { 

  # # Create empty vectors to store results
  # eff <- numeric(length(intervention[["results"]]))
  # dur <- numeric(length(intervention[["results"]]))
  # averted_cases <- numeric()
  # averted_cases_per_dose <- numeric()
  # averted_cases_percentage <- numeric()
  # median_averted_case <- numeric(length(intervention[["results"]]))
  # 
  # median_perdose <- numeric(length(intervention[["results"]]))
  # median_percentage <- numeric(length(intervention[["results"]]))
  # front_perdose <- numeric(length(intervention[["results"]]))
  # front_percentage <- numeric(length(intervention[["results"]]))
  # tail_perdose <- numeric(length(intervention[["results"]]))
  # tail_percentage <- numeric(length(intervention[["results"]]))
  # 
  # 
  # # Iterate through each list element
  # for (i in seq_along(intervention[["results"]])) {
  #   result_key <- names(intervention[["results"]][i])
  #   # Split names into eff and dur
  #   eff_dur <- strsplit(result_key, "_dur")
  #   eff[i] <- sapply(eff_dur, function(x) as.numeric(sub("eff", "", x[1])))
  #   dur[i] <- sapply(eff_dur, function(x) as.numeric(sub("_dur", "", x[2])))
  #   
  #   # Retrieve wanted metrics by samples in 10 year expectation
  #   cases_averted_per_dose <- intervention[["results"]][[i]][["cases_averted_per_dose"]][10,]
  #   inc_cum_doses <- intervention[["results"]][[i]][["inc_cum_doses"]][10,]
  #   averted_case <- cases_averted_per_dose * inc_cum_doses
  #   averted_cases_percentage <- ((novax_2044_flows[["cum_treated"]][10,]-intervention[["results"]][[i]][["cum_treated"]][10,])/novax_2044_flows[["cum_treated"]][10,])
  #   if (any(is.na(averted_case))) {
  #     averted_case <- inc_cum_doses
  #   }
  #   
  #   if (title == "VbE") {
  #     singapore_male_entrance <- demographic_params_sg$male_enr
  #     vaccinated <- singapore_male_entrance*10*control$vbe
  #     incidence <- -intervention_2044[["results"]][[i]]$inc_cum_treated[elapsed,]
  #     cases_averted_per_dose <- incidence/(vaccinated*2) 
  #   }
  #   # averted_cases <- c(averted_cases, averted_case)
  #   # averted_cases_per_dose <- c(averted_cases_per_dose, cases_averted_per_dose)
  #   # averted_cases_percentage <- c(averted_cases_percentage, ((novax_2044_flows[["cum_treated"]][10,]-intervention[["results"]][[i]][["cum_treated"]][10,])/novax_2044_flows[["cum_treated"]][10,])*100)
  #   
  #   # Calculate stats
  #   perdose <- stats_summary(cases_averted_per_dose)
  #   percentage <- stats_summary(averted_cases_percentage)
  #   
  #   # median_averted_case[i] <- median(averted_case)
  #   #front_perdose[i] <- perdose [[1]]
  #   median_perdose[i] <- perdose [[2]]
  #   #tail_perdose[i] <- perdose [[3]]
  #   
  #   #front_percentage[i] <- percentage[[1]]
  #   median_percentage[i] <- percentage [[2]]
  #   #tail_percentage[i] <- percentage [[3]]
  # }
  # 
  # 
  # # Create the data frame with the modified object naming convention
  # results_summary <- data.frame(eff = eff, 
  #                               dur = dur, 
  #                               uptake = uptake_total,
  #                               median_perdose = median_perdose,
  #                               median_percentage = median_percentage,
  #                               strategy = title)
  # 
  } 
  
  if (graph == "trivariates") {
  
    # # Initialize lists to store data
    # eff <- numeric(length(intervention[["results"]]))
    # dur <- numeric(length(intervention[["results"]]))
    # uptake <- numeric(length(intervention[["results"]]))
    # result_details_list <- list()
    # 
    # # Iterate through each list element
    # for (i in seq_along(intervention[["results"]])) {
    #   result_key <- names(intervention[["results"]][i])
    #   # Split names into eff and dur
    #   eff_dur <- strsplit(result_key, "_dur")
    #   eff[i] <- sapply(eff_dur, function(x) as.numeric(sub("eff", "", x[1])))
    #   dur[i] <- sapply(eff_dur, function(x) as.numeric(sub("_dur", "", x[2])))
    #   
    #   # Retrieve wanted metrics
    #   cases_averted_per_dose <- intervention[["results"]][[i]][["cases_averted_per_dose"]][elapsed,]
    #   inc_cum_doses <- intervention[["results"]][[i]][["inc_cum_doses"]][elapsed,]
    #   inc_cum_treated <- intervention [["results"]][[i]][["inc_cum_treated"]][elapsed,]
    #   averted_case <- cases_averted_per_dose * inc_cum_doses
    #   averted_cases_percentage <- ((novax_2044_flows[["cum_treated"]][elapsed,]-intervention[["results"]][[i]][["cum_treated"]][elapsed,])/novax_2044_flows[["cum_treated"]][elapsed,])
    #   if (any(is.na(averted_case))) {
    #     averted_case <- inc_cum_doses
    #   }
    #   
    #   # Create a data frame for each result
    #   result_df <- data.frame(eff = eff[i], 
    #                           dur = dur[i], 
    #                           uptake = uptake_total,
    #                           Averted_cases = averted_case,
    #                           verify = inc_cum_treated,
    #                           Averted_Case_Per_Dose = cases_averted_per_dose,
    #                           Averted_Case_Percentage = averted_cases_percentage,
    #                           strategy = title)
    #   result_details_list[[i]] <- result_df
    # }
    # 
    # # Combine all data frames into one
    # results_details <- do.call(rbind, result_details_list)
    # 
    # # Rename the data frame according to the naming convention
    # object_name <- paste(title, "_", uptake_title, "_uptake_all_vax_details", sep = "")
    # assign(object_name, results_details)
    # 
    # print(paste0(object_name, " analysis finished"))
  }
}

