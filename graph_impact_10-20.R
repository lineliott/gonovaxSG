# This script produces figure 2 (incidence projection & impact radar graph)
library(ggplot2)
library(tidyr)
library(reshape2)
library(lemon)
library(fmsb)
library(patchwork)


run_grid_control <- function(short_run, title, revax) {
  vbe_uptake <- 0 # no vbe background
  
  if (short_run) {
    eff <- c(0.31)
    dur <- c(4)
  } else {
    eff <- seq(0.05, 1, 0.05)
    dur <- seq_len(20)
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

# Reshape data to long format
convert_to_long <- function(df, title) {
  df$year <- as.numeric(rownames(df))
  df_long <- tidyr::gather(data = df, key = "variable", value = "value", -year)
  return(df_long)
}

# decimal transformation function
scaleFUN <- function(x) sprintf("%.2f", x)


# create file dir
path <- "calibration/graph/singlevax/impact_10-20/"
dir.create(path, FALSE, TRUE)
strategy_list <- c("VbE", "VoD(H)", "VoD", "VoA", "VoA(H)", "VoD+VoA(H)")
for (i in strategy_list) {
  short_run <- TRUE
  uptake_total<- 0.33
  graph <- "timeseries"
  vaccine_algorithm <- "xvwv"
  fix_par_t <- TRUE
  stream <- "msm"
  strategy <- if (i == "VbE") NULL else i
  # run result
  source("script_projection.R")
  ##############################################################################
  # boxplot
  # intervention all time projection
  intervention_2044 <- intervention$results[[key]]
  intervention_single_flow <- intervention_2044
  # glue together
  intervention_single_flow <- lapply(names(novax_2024_flows), function(name) {
    rbind(novax_2024_flows[[name]], intervention_single_flow[[name]])
  })
  names(intervention_single_flow) <- names(novax_2024_flows)
  # as data frame 
  flow_intervention_treated <- as.data.frame(intervention_single_flow$treated)
  rownames(flow_intervention_treated) <- c(2004:2045)
  
  # novax all time projection
  novax_single_flow <- lapply(names(novax_2024_flows), function(name) {
    rbind(novax_2024_flows[[name]], novax_2044_flows[[name]])
  })
  names(novax_single_flow) <- names(novax_2024_flows)
  
  flow_novax_treated <- as.data.frame(novax_single_flow$treated)
  rownames(flow_novax_treated) <- c(2004:2045)
  
  
  flow_novax_treated <- convert_to_long(flow_novax_treated, title)
  flow_novax_treated$type <- "baseline"
  flow_intervention_treated <- convert_to_long(flow_intervention_treated, title)
  flow_intervention_treated$type <- "intervention"
  flow_combined <- rbind(flow_novax_treated, flow_intervention_treated)
  flow_combined$strategy <- title
  
  object_name <- paste(path,title,"_flow_combined.rds", sep = "")
  saveRDS(flow_combined, object_name)
  
  ####################################################################################
  # radar_10
  elapsed <- 10
  novax_treated <- novax_2044_flows$cum_treated[elapsed+1,]
  intervention_treated <- intervention_2044$cum_treated[elapsed,]
  averted_cases_percentage <- (novax_treated - intervention_treated)/novax_treated
  
  averted_cases_per_dose <- intervention_2044$cases_averted_per_dose[elapsed,]
  vaccinated <- intervention_2044$cum_vaccinated[elapsed,]
  primaryratio <- intervention_2044$inc_cum_primary[elapsed,]/intervention_2044$cum_vaccinated[elapsed,]
  if (i == "VbE"){
    singapore_male_entrance <- demographic_params_sg$male_enr
    vaccinated <- singapore_male_entrance*elapsed*control$vbe
    averted_cases <- -intervention_2044$inc_cum_treated[elapsed,]
    cases_averted_per_dose <- averted_cases/(vaccinated*2)
    primaryratio <- intervention_2044$cum_vaccinated[elapsed,]/intervention_2044$cum_vaccinated[elapsed,]
  }
  radar_10 <- data.frame(
    #novax_treated = novax_treated,
    #intervention_treated = intervention_treated,
    percentage = averted_cases_percentage,
    perdose = averted_cases_per_dose,
    vaccinated = vaccinated,
    primaryratio = primaryratio,
    strategy = title
  )
  object_name <- paste(path, title,"_radar_10.rds", sep = "")
  saveRDS(radar_10, object_name)
  
##################################################################################
# radar_20
  elapsed <- 20
  novax_treated <- novax_2044_flows$cum_treated[elapsed+1,]
  intervention_treated <- intervention_2044$cum_treated[elapsed,]
  averted_cases_percentage <- (novax_treated - intervention_treated)/novax_treated
  
  averted_cases_per_dose <- intervention_2044$cases_averted_per_dose[elapsed,]
  vaccinated <- intervention_2044$cum_vaccinated[elapsed,]
  primaryratio <- intervention_2044$inc_cum_primary[elapsed,]/intervention_2044$cum_vaccinated[elapsed,]
if (i == "VbE"){
  singapore_male_entrance <- demographic_params_sg$male_enr
  vaccinated <- singapore_male_entrance*elapsed*control$vbe
  averted_cases <- -intervention_2044$inc_cum_treated[elapsed,]
  cases_averted_per_dose <- averted_cases/(vaccinated*2)
  primaryratio <- intervention_2044$cum_vaccinated[elapsed,]/intervention_2044$cum_vaccinated[elapsed,]
}
radar_20 <- data.frame(
  #novax_treated = novax_treated,
  #intervention_treated = intervention_treated,
  percentage = averted_cases_percentage,
  perdose = averted_cases_per_dose,
  vaccinated = vaccinated,
  primaryratio = primaryratio,
  strategy = title
)
object_name <- paste(path, title,"_radar_20.rds", sep = "")
saveRDS(radar_20, object_name)

print(paste0(i, " analysis finished"))
}
###################################################################################
VbE_flow_combined <- readRDS(paste0(path,"VbE_flow_combined.rds"))
VoH_flow_combined <- readRDS(paste0(path,"VoD_flow_combined.rds"))
VoA_flow_combined <- readRDS(paste0(path,"VoA_flow_combined.rds"))
`VoA(H)_flow_combined` <- readRDS(paste0(path,"VoA(H)_flow_combined.rds"))
`VoD(H)_flow_combined`<- readRDS(paste0(path,"VoD(H)_flow_combined.rds"))
VaR_flow_combined <- readRDS(paste0(path,"VaR_flow_combined.rds"))


flow_combined <- rbind(VbE_flow_combined,
                       VoH_flow_combined,
                       VoA_flow_combined,
                       `VoA(H)_flow_combined`,
                       `VoD(H)_flow_combined`,
                       VaR_flow_combined)
# subset interested time period
flow_combined <- subset(flow_combined, year >= 2024 & year <= 2044)

# Calculate mean of value by strategy and year
mean_combined <- subset(flow_combined, type == "intervention") %>%
  group_by(strategy, year) %>%
  summarise(mean_value = ifelse(mean(value) < 0.001, mean(value), NaN))

mean_combined$type <- "intervention"


VbE_radar_10<- readRDS(paste0(path,"VbE_radar_10.rds"))
VoH_radar_10 <- readRDS(paste0(path,"VoD_radar_10.rds"))
VoA_radar_10 <- readRDS(paste0(path,"VoA_radar_10.rds"))
`VoA(H)_radar_10` <- readRDS(paste0(path,"VoA(H)_radar_10.rds"))
`VoD(H)_radar_10`<- readRDS(paste0(path,"VoD(H)_radar_10.rds"))
VaR_radar_10 <- readRDS(paste0(path,"VaR_radar_10.rds"))


radar_combined_10 <- rbind(VbE_radar_10,
                        VoH_radar_10,
                        VoA_radar_10,
                        `VoA(H)_radar_10`,
                        `VoD(H)_radar_10`,
                        VaR_radar_10)


VbE_radar_20<- readRDS(paste0(path,"VbE_radar_20.rds"))
VoH_radar_20 <- readRDS(paste0(path,"VoD_radar_20.rds"))
VoA_radar_20 <- readRDS(paste0(path,"VoA_radar_20.rds"))
`VoA(H)_radar_20` <- readRDS(paste0(path,"VoA(H)_radar_20.rds"))
`VoD(H)_radar_20`<- readRDS(paste0(path,"VoD(H)_radar_20.rds"))
VaR_radar_20 <- readRDS(paste0(path,"VaR_radar_20.rds"))


radar_combined_20 <- rbind(VbE_radar_20,
                           VoH_radar_20,
                           VoA_radar_20,
                           `VoA(H)_radar_20`,
                           `VoD(H)_radar_20`,
                           VaR_radar_20)

# Specify the strategy order
strategy_order <- c("VbE", "VoD(H)","VoD", "VoA", "VoA(H)", "VaR")
flow_combined$strategy <- factor(flow_combined$strategy, levels = strategy_order)
# mean_combined$strategy <- factor(mean_combined$strategy, levels = strategy_order)
# scatter_combined$strategy <- factor(scatter_combined$strategy, levels = strategy_order)

plot1 <- ggplot(flow_combined, 
                aes(x = factor(year), 
                    y = value/1000, 
                    fill = factor(type),
                    col= factor(type))) + 
  stat_boxplot_custom(qs = c(0.025, 0.025, 0.5, 0.975, 0.975),
                      outlier.shape = NA,
                      alpha = 0.7,
                      #color = "#00000000",
                      position = "identity")+
  ylab("Annual number of diagnosed (thousand)")+
  stat_summary(fun = mean, 
               geom ="point", 
               shape = 18, 
               size = 6) +
  # geom_point(data = subset(mean_combined, !is.na(mean_combined$mean_value)), 
  #            aes(x = factor(year), y = mean_value, fill = "Eliminated"),
  #            shape = 23,
  #            size = 4)+
  scale_color_manual(name = "Projection",
                     labels = c(#"Eliminated"
                       "Baseline", "Intervention"),
                     values= c(baseline = "brown",
                               intervention = "deepskyblue4")) +
  scale_fill_manual(name = "Projection",
                    labels = c(#"Eliminated"
                      "Baseline", "Intervention"),
                    values = c(baseline = "darksalmon",
                               #Eliminated = "darkgreen",
                               intervention = "cadetblue1")) +
  scale_x_discrete(breaks = c(2024, 2034, 2044)) + 
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 1), labels=scaleFUN) +
  coord_capped_cart(bottom='both', left='both') +
  facet_rep_wrap(~strategy, ncol = 6, scales = "free") + 
  theme(
    axis.title.y = element_text(size = 30),
    axis.title.x = element_blank(),
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
    panel.border=element_blank(), axis.line=element_line()
  )


# Specify the strategy order
strategy_order <- c("VbE", "VoD(H)","VoD", "VoA", "VoA(H)", "VaR")



normalize_df <- function(df) {
  # Extract the strategy column
  strategy_col <- df$strategy
  
  # Remove the strategy column for normalization
  df <- df[, !(names(df) %in% c("strategy"))]
  
  # Normalize each column by its maximum value
  normalized_df <- as.data.frame(lapply(df, function(x)((x-min(x)) / (max(x)-min(x)))))
  
  # Add the strategy column back
  normalized_df <- cbind(strategy = strategy_col, normalized_df)
  
  return(normalized_df)
}

# Normalize radar
radar_normalized_10 <- normalize_df(radar_combined_10)
radar_normalized_20 <- normalize_df(radar_combined_20)
# Create a sample ID column
radar_normalized_10$sample_id <- seq_len(nrow(radar_normalized_10))
radar_normalized_20$sample_id <- seq_len(nrow(radar_normalized_20))
# Rename colnames
radar_normalized_10 <- radar_normalized_10 %>%
  rename(A.P. = percentage) %>%
  rename(A.D. = perdose) %>%
  rename(N.V. = vaccinated) %>%
  rename(P.U. = primaryratio)
radar_normalized_20 <- radar_normalized_20 %>%
  rename(A.P. = percentage) %>%
  rename(A.D. = perdose) %>%
  rename(N.V. = vaccinated) %>%
  rename(P.U. = primaryratio)


# Specify the strategy order
strategy_order <- c("VbE", "VoD(H)","VoD", "VoA", "VoA(H)", "VaR")
radar_normalized_10$strategy <- factor(radar_normalized_10$strategy, levels = strategy_order)
radar_normalized_20$strategy <- factor(radar_normalized_20$strategy, levels = strategy_order)

# Convert to long data frame
radar_long_10 <-
  radar_normalized_10 %>%
  melt(id.vars=c("sample_id", "strategy")) %>%
  rbind(subset(., variable == names(radar_normalized_10)[2]))

radar_long_20 <-
  radar_normalized_20 %>%
  melt(id.vars=c("sample_id", "strategy")) %>%
  rbind(subset(., variable == names(radar_normalized_20)[2]))

# Calculate mean values for each variable
means_10 <- radar_long_10 %>%
  group_by(variable,strategy) %>%
  summarize(mean_value = mean(value, na.rm = TRUE))
means_20 <- radar_long_20 %>%
  group_by(variable,strategy) %>%
  summarize(mean_value = mean(value, na.rm = TRUE))


mean_line_10 <- radar_normalized_10 %>%
  group_by(strategy) %>%
  summarize(across(everything(), mean, na.rm = TRUE)) %>%
  melt(id.vars=c("sample_id", "strategy")) %>%
  rbind(subset(., variable == names(radar_normalized_10)[2]))
mean_line_20 <- radar_normalized_20 %>%
  group_by(strategy) %>%
  summarize(across(everything(), mean, na.rm = TRUE)) %>%
  melt(id.vars=c("sample_id", "strategy")) %>%
  rbind(subset(., variable == names(radar_normalized_20)[2]))



# create new coord : inherit coord_polar
coord_radar <- 
  function(theta='x', start=0, direction=1){
    # input parameter sanity check
    match.arg(theta, c('x','y'))
    
    ggproto(
      NULL, CoordPolar, 
      theta=theta, r=ifelse(theta=='x','y','x'),
      start=start, direction=sign(direction),
      is_linear=function() TRUE)
  }


# Plot with ggplot2
plot2 <- ggplot(radar_long_10, aes(x = factor(variable), y = value, group = sample_id)) +
  geom_path(aes(color = "samples", alpha = "samples"), linewidth = 2) +
  geom_point(data = means_10, 
             aes(x = factor(variable), y = mean_value, color = "mean"), 
             size = 6, shape = 18, inherit.aes = FALSE) +
  geom_path(data = mean_line_10,
            aes(x = factor(variable), y = value, group = sample_id, color = "mean"),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Score(10 y')",
                     values = c("samples" = "aquamarine4", "mean" = "darkgreen"),
                     labels = c("samples" = "Sample", "mean" = "Mean")) +
  scale_alpha_manual(values = c("samples"= 0.01), guide = "none") +
  #labs(x = "Sample ID", y = "Values", title = "Line Graph by Sample ID") +  
  coord_radar() +
  facet_rep_wrap(~strategy, ncol = 6, scales = "fixed") +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 30, hjust = 0.5),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_line(color = "grey"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(50, "pt"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    panel.border = element_blank()
  )

plot3 <- ggplot(radar_long_20, aes(x = factor(variable), y = value, group = sample_id)) +
  geom_path(aes(color = "samples", alpha = "samples"), linewidth = 2) +
  geom_point(data = means_20, 
             aes(x = factor(variable), y = mean_value, color = "mean"), 
             size = 6, shape = 18, inherit.aes = FALSE) +
  geom_path(data = mean_line_20,
            aes(x = factor(variable), y = value, group = sample_id, color = "mean"),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Score (20 y')",
                     values = c("samples" = "aquamarine4", "mean" = "darkgreen"),
                     labels = c("samples" = "Sample", "Mean" = "mean")) +
  scale_alpha_manual(values = c("samples"= 0.01), guide = "none") +
  #labs(x = "Sample ID", y = "Values", title = "Line Graph by Sample ID") +  
  coord_radar() +
  facet_rep_wrap(~strategy, ncol = 6, scales = "fixed") +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 30, hjust = 0.5),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_line(color = "grey"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(50, "pt"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    panel.border = element_blank()
  )
# Save the plot
#ggsave(file = paste0(path,"radar.png"), plot = plot2, width = 48, height = 8, dpi = 320)

# Combine strategies into a big plot
big_plot <- plot1 + plot2 + plot3 + plot_layout(ncol=1)

# Save the plot
#ggsave(file = paste0(path, "impact_combined.png"), plot = big_plot, width = 48, height = 16, dpi = 320)

ggsave(file = paste0(path,"impact_plot.pdf"), big_plot, width = 48, height = 24)

