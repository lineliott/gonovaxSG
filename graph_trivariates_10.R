# This script produces figure 4 (boxplots: senstivity analysis)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lemon)
library(patchwork)
source("script_calibration_plot.R")

# decimal transformation function
scaleFUN <- function(x) sprintf("%.1f", x)

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
    eff <- c(0.22, 0.31, 0.47)
    dur <- c(1.5, 4, 7.5)
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
       model = model, 
       vbe = vbe_uptake)
}

# create file dir
path <- "calibration/graph/multivax/trivariates_10/"
dir.create(path, FALSE, TRUE)
strategy_list <- c("VbE", "VoD(H)", "VoD", "VoA", "VoA(H)", "VoD+VoA(H)")
uptake_rate_list <- list("low"= 0.1, "mid"= 0.33, "high"= 1.0)

for (i in strategy_list) {
  short_run <- FALSE
  graph <- "trivariates"
  fix_par_t <- TRUE
  vaccine_algorithm <- "xvwv"
  elapsed <- 10
  stream <- "msm"
  strategy <- if (i == "VbE") NULL else i
  for (j in 1:length(uptake_rate_list)){
    uptake_total <- uptake_rate_list[[j]]
    uptake_title <- names(uptake_rate_list[j])
    if (i == "VbE") {
      run_grid_control <- function(short_run, title, revax) {
        vbe_uptake <- 0 # no vbe background
        
        if (short_run) {
          eff <- c(0.31)
          dur <- c(4)
        } else {
          eff <- c(0.22, 0.31, 0.47)
          dur <- c(1.5, 4, 7.5)
        }
        if (revax) {
          model <- run_onevax_xvwv
        } else {
          stop("no-revaccination not supported")
        }
        
        if (title == "VbE") {
          vbe_uptake <- uptake_total # vbe strategy
        }
        
        list(eff = eff, 
             dur = dur,
             model = model, 
             vbe = vbe_uptake)
      }
    }
    # run result
    source("script_projection.R")
    
    # process
    
    # Initialize lists to store data
    eff <- numeric(length(intervention[["results"]]))
    dur <- numeric(length(intervention[["results"]]))
    uptake <- numeric(length(intervention[["results"]]))
    result_details_list <- list()
    
    # Iterate through each list element
    for (i in seq_along(intervention[["results"]])) {
      result_key <- names(intervention[["results"]][i])
      # Split names into eff and dur
      eff_dur <- strsplit(result_key, "_dur")
      eff[i] <- sapply(eff_dur, function(x) as.numeric(sub("eff", "", x[1])))
      dur[i] <- sapply(eff_dur, function(x) as.numeric(sub("_dur", "", x[2])))
      
      # Retrieve wanted metrics
      cases_averted_per_dose <- intervention[["results"]][[i]][["cases_averted_per_dose"]][elapsed,]
      inc_cum_doses <- intervention[["results"]][[i]][["inc_cum_doses"]][elapsed,]
      inc_cum_treated <- intervention [["results"]][[i]][["inc_cum_treated"]][elapsed,]
      averted_case <- cases_averted_per_dose * inc_cum_doses
      averted_cases_percentage <- ((novax_2044_flows[["cum_treated"]][elapsed,]-intervention[["results"]][[i]][["cum_treated"]][elapsed,])/novax_2044_flows[["cum_treated"]][elapsed,])
      if (any(is.na(averted_case))) {
        averted_case <- inc_cum_doses
      }
      if (title == "VbE") {
        singapore_male_entrance <- demographic_params_sg$male_enr
        vaccinated <- singapore_male_entrance*elapsed*control$vbe
        incidence <- -intervention[["results"]][[i]][["inc_cum_treated"]][elapsed,]
        cases_averted_per_dose <- incidence/(vaccinated*2)
      }
      
      
      # Create a data frame for each result
      result_df <- data.frame(eff = eff[i], 
                              dur = dur[i], 
                              uptake = uptake_total,
                              Averted_cases = averted_case,
                              verify = inc_cum_treated,
                              Averted_Case_Per_Dose = cases_averted_per_dose,
                              Averted_Case_Percentage = averted_cases_percentage,
                              strategy = title)
      result_details_list[[i]] <- result_df
    }
    
    # Combine all data frames into one
    results_details <- do.call(rbind, result_details_list)
    
    # Rename the data frame according to the naming convention
    object_name <- paste(title, "_", uptake_title, "_uptake_all_vax_details", sep = "")
    assign(object_name, results_details)
    
    print(paste0(object_name, " analysis finished"))
  }
}



saveRDS(VoA_low_uptake_all_vax_details, paste0(path,"VoA_low_uptake_all_vax_details.rds"))
saveRDS(`VaR_low_uptake_all_vax_details`, paste0(path,"VaR_low_uptake_all_vax_details.rds"))
saveRDS(VoD_low_uptake_all_vax_details, paste0(path,"VoD_low_uptake_all_vax_details.rds"))
saveRDS(`VoD(H)_low_uptake_all_vax_details`, paste0(path,"VoD(H)_low_uptake_all_vax_details.rds"))
#saveRDS(novax_low_uptake_all_vax_details, paste0(path,"novax_low_uptake_all_vax_details.rds"))
saveRDS(VbE_low_uptake_all_vax_details, paste0(path,"VbE_low_uptake_all_vax_details.rds"))
saveRDS(`VoA(H)_low_uptake_all_vax_details`, paste0(path,"VoA(H)_low_uptake_all_vax_details.rds"))

saveRDS(VoA_mid_uptake_all_vax_details, paste0(path,"VoA_mid_uptake_all_vax_details.rds"))
saveRDS(`VaR_mid_uptake_all_vax_details`, paste0(path,"VaR_mid_uptake_all_vax_details.rds"))
saveRDS(VoD_mid_uptake_all_vax_details, paste0(path,"VoD_mid_uptake_all_vax_details.rds"))
saveRDS(`VoD(H)_mid_uptake_all_vax_details`, paste0(path,"VoD(H)_mid_uptake_all_vax_details.rds"))
#saveRDS(novax_mid_uptake_all_vax_details, paste0(path,"novax_mid_uptake_all_vax_details.rds"))
saveRDS(VbE_mid_uptake_all_vax_details, paste0(path,"VbE_mid_uptake_all_vax_details.rds"))
saveRDS(`VoA(H)_mid_uptake_all_vax_details`, paste0(path,"VoA(H)_mid_uptake_all_vax_details.rds"))

saveRDS(VoA_high_uptake_all_vax_details, paste0(path,"VoA_high_uptake_all_vax_details.rds"))
saveRDS(`VaR_high_uptake_all_vax_details`, paste0(path,"VaR_high_uptake_all_vax_details.rds"))
saveRDS(VoD_high_uptake_all_vax_details, paste0(path,"VoD_high_uptake_all_vax_details.rds"))
saveRDS(`VoD(H)_high_uptake_all_vax_details`, paste0(path,"VoD(H)_high_uptake_all_vax_details.rds"))
#saveRDS(novax_high_uptake_all_vax_details, paste0(path,"novax_high_uptake_all_vax_details.rds"))
saveRDS(VbE_high_uptake_all_vax_details, paste0(path,"VbE_high_uptake_all_vax_details.rds"))
saveRDS(`VoA(H)_high_uptake_all_vax_details`, paste0(path,"VoA(H)_high_uptake_all_vax_details.rds"))
#####################################################################################################
# read everything in
VoA_low_uptake_all_vax_details <- readRDS(paste0(path,"VoA_low_uptake_all_vax_details.rds"))
`VaR_low_uptake_all_vax_details` <- readRDS(paste0(path,"VaR_low_uptake_all_vax_details.rds"))
VoD_low_uptake_all_vax_details <- readRDS(paste0(path,"VoD_low_uptake_all_vax_details.rds"))
`VoD(H)_low_uptake_all_vax_details` <- readRDS(paste0(path,"VoD(H)_low_uptake_all_vax_details.rds"))
#novax_low_uptake_all_vax_details <- readRDS(paste0(path,"novax_low_uptake_all_vax_details.rds"))
VbE_low_uptake_all_vax_details <- readRDS(paste0(path,"VbE_low_uptake_all_vax_details.rds"))
`VoA(H)_low_uptake_all_vax_details` <- readRDS(paste0(path,"VoA(H)_low_uptake_all_vax_details.rds"))

VoA_mid_uptake_all_vax_details <- readRDS(paste0(path,"VoA_mid_uptake_all_vax_details.rds"))
`VaR_mid_uptake_all_vax_details` <- readRDS(paste0(path,"VaR_mid_uptake_all_vax_details.rds"))
VoD_mid_uptake_all_vax_details <- readRDS(paste0(path,"VoD_mid_uptake_all_vax_details.rds"))
`VoD(H)_mid_uptake_all_vax_details` <- readRDS(paste0(path,"VoD(H)_mid_uptake_all_vax_details.rds"))
#novax_mid_uptake_all_vax_details <- readRDS(paste0(path,"novax_mid_uptake_all_vax_details.rds"))
VbE_mid_uptake_all_vax_details <- readRDS(paste0(path,"VbE_mid_uptake_all_vax_details.rds"))
`VoA(H)_mid_uptake_all_vax_details` <- readRDS(paste0(path,"VoA(H)_mid_uptake_all_vax_details.rds"))

VoA_high_uptake_all_vax_details <- readRDS(paste0(path,"VoA_high_uptake_all_vax_details.rds"))
`VaR_high_uptake_all_vax_details` <- readRDS(paste0(path,"VaR_high_uptake_all_vax_details.rds"))
VoD_high_uptake_all_vax_details <- readRDS(paste0(path,"VoD_high_uptake_all_vax_details.rds"))
`VoD(H)_high_uptake_all_vax_details` <- readRDS(paste0(path,"VoD(H)_high_uptake_all_vax_details.rds"))
#novax_high_uptake_all_vax_details <- readRDS(paste0(path,"novax_high_uptake_all_vax_details.rds"))
VbE_high_uptake_all_vax_details <- readRDS(paste0(path,"VbE_high_uptake_all_vax_details.rds"))
`VoA(H)_high_uptake_all_vax_details` <- readRDS(paste0(path,"VoA(H)_high_uptake_all_vax_details.rds"))



# Combine dataframes
combined_df <- bind_rows(VoA_low_uptake_all_vax_details,
                         `VaR_low_uptake_all_vax_details`, 
                         VoD_low_uptake_all_vax_details, 
                         `VoD(H)_low_uptake_all_vax_details`,
                         VbE_low_uptake_all_vax_details, 
                         `VoA(H)_low_uptake_all_vax_details`,
                         
                         VoA_mid_uptake_all_vax_details,
                         `VaR_mid_uptake_all_vax_details`, 
                         VoD_mid_uptake_all_vax_details, 
                         `VoD(H)_mid_uptake_all_vax_details`,
                         VbE_mid_uptake_all_vax_details, 
                         `VoA(H)_mid_uptake_all_vax_details`,
                         
                         VoA_high_uptake_all_vax_details,
                         `VaR_high_uptake_all_vax_details`, 
                         VoD_high_uptake_all_vax_details, 
                         `VoD(H)_high_uptake_all_vax_details`,
                         VbE_high_uptake_all_vax_details,
                         `VoA(H)_high_uptake_all_vax_details`,
                         
                         .id = "id")


combined_df <- combined_df %>%
  mutate(eff = case_when(
    eff == min(eff) ~ "M",   # mild
    eff == median(eff) ~ "N",# normal
    eff == max(eff) ~ "S"   # strong
  ))

combined_df <- combined_df %>%
  mutate(dur = case_when(
    dur == min(dur) ~ "1.5 years",
    dur == median(dur) ~ "4 years",
    dur == max(dur) ~ "7.5 years"
  ))

combined_df <- combined_df %>%
  mutate(uptake = case_when(
    uptake == min(uptake) ~ "Low = 10%",
    uptake == median(uptake) ~ "Central = 33%",
    uptake == max(uptake) ~ "High = 100%"
  ))





strategies <- c(#"novax",
  "VbE", "VoD(H)", "VoD", "VoA", "VoA(H)", "VaR" )


fill_order<- c("Low = 10%", "Central = 33%", "High = 100%")
text_font_size <- 30


for (i in strategies) {
  # Filter data for the current strategy
  df <- combined_df %>%
    filter(strategy == i)
  
  # Plot Averted Cases Percentage
  plot1 <- ggplot(df, aes(x = factor(eff), y = Averted_Case_Percentage, 
                          fill = factor(uptake, level = fill_order),
                          color = factor(uptake, level = fill_order))) +
    stat_boxplot_custom(qs = c(0.025, 0.025, 0.5, 0.975, 0.975),
                        outlier.shape = NA,
                        alpha = 0.5)+
    # geom_boxplot(outlier.shape = NA,
    #              alpha = 1.0,
    #              color = "black") +
    scale_fill_manual(values = c("Low = 10%" = rgb(243,162,097, maxColorValue = 255), 
                                 "Central = 33%" = rgb(042,157,140, maxColorValue = 255), 
                                 "High = 100%" = rgb(075,116,178, maxColorValue = 255)),
                      name = "Uptake") +
    scale_color_manual(values = c("Low = 10%" = "darkred",
                                  "Central = 33%" = "darkgreen",
                                  "High = 100%" = "deepskyblue4"),
                       name = "Uptake") +
    ggtitle(i)+
    facet_wrap(~ dur, scales = "free_x", nrow = 1, strip.position = "top") +
    # mean annotation
    stat_summary(fun.y = mean, 
                 position = position_dodge2(width = 0.75),
                 geom ="point", 
                 shape = 18, 
                 size = 5,
                 aes(color = factor(uptake, level = fill_order))) +
    labs(x = NULL, y = ifelse(i == strategies[1], "Averted cases by percentage", ""), fill = "Uptake") +
    scale_y_continuous(limits = c(0, 1), n.breaks = 4, label = scaleFUN) +
    theme(          
      # Title adjustment
      plot.title = element_text(size = text_font_size, hjust = 0.5),
      # Hide panel borders and remove grid lines
      panel.background = element_blank(),
      # Add axis line
      axis.line = element_line(colour = "grey"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = text_font_size),
      axis.text.y = element_text(size = text_font_size),
      
      # Facet wrap adjustment
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      
      # Set the legend position based on the branch
      legend.position = if (i == strategies[1]) c(0.3, 0.8) else "none",
      legend.text = element_text(size = text_font_size),  # Increase legend text size
      legend.title = element_text(size = text_font_size),  # Increase legend title size
      plot.margin = margin(t = 10, r = 10, b = 20, l = 10))
  
  # Plot Averted Cases Per Dose
  plot2 <- ggplot(df, aes(x = factor(eff), 
                          y = Averted_Case_Per_Dose, 
                          fill = factor(uptake, level = fill_order),
                          color = factor(uptake, level = fill_order))) +
    stat_boxplot_custom(qs = c(0.025, 0.025, 0.5, 0.975, 0.975),
                        outlier.shape = NA,
                        alpha = 0.5)+
    scale_fill_manual(values = c("Low = 10%" = rgb(243,162,097, maxColorValue = 255), 
                                 "Central = 33%" = rgb(042,157,140, maxColorValue = 255), 
                                 "High = 100%" = rgb(075,116,178, maxColorValue = 255)),
                      name = "Uptake") +
    scale_color_manual(values = c("Low = 10%" = "darkred",
                                  "Central = 33%" = "darkgreen",
                                  "High = 100%" = "deepskyblue4"),
                       name = "Uptake") +
    facet_wrap(~ dur, scales = "fixed", nrow = 1, strip.position = "bottom") +
    # mean annotation
    stat_summary(fun.y = mean, 
                 position = position_dodge2(width = 0.75),
                 geom ="point", 
                 shape = 18, 
                 size = 5,
                 aes(color = factor(uptake, level = fill_order))) +
    labs(x = "Vaccine efficacy-protection duration", y = ifelse(i == strategies[1], "Averted cases per dose", ""), fill = "Uptake") +
    scale_y_continuous(limits = c(0, 1.0), n.breaks = 3, label=scaleFUN) +
    theme(      
      # Hide panel borders and remove grid lines
      panel.background = element_blank(),
      # Add axis line
      axis.line = element_line(colour = "grey"),
      # X axis title adjustment
      axis.title.x.bottom = element_text(size = text_font_size, hjust = 0.5),
      axis.text.x = element_text(size = text_font_size),
      axis.title.y = element_text(size = text_font_size),
      axis.text.y = element_text(size = text_font_size),
      # Facet wrap adjustment
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = text_font_size),  

      # Set the legend position based on the branch
      legend.position = if (i == strategies[1]) c(0.3, 0.8) else "none",
      legend.text = element_text(size = text_font_size),  # Increase legend text size
      legend.title = element_text(size = text_font_size),  # Increase legend title size
      plot.margin = margin(t = 10, r = 10, b = 20, l = 10))
  
  
  # Combine three metrices in a row
  gA <- ggplotGrob(plot1)
  gB <- ggplotGrob(plot2)
  grid::grid.newpage()
  big_plot <- rbind(gA, gB)
  
  object_name <- paste(i,"boxplot", sep = "")
  assign(object_name, big_plot)
}


# Combine strategies into a big plot
big_plot <- cbind(VbEboxplot, VoDboxplot, `VoD(H)boxplot`, VoAboxplot,`VoA(H)boxplot`,`VaRboxplot`)

# Save the plot
# ggsave(file = paste0(path,"combined_plot.png"), big_plot, width = 48, height = 16, dpi = 320)

ggsave(file = paste0(path,"trivariates_plot.pdf"), big_plot, width = 48, height = 16)


