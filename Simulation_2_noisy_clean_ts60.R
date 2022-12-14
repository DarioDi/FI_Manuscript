#Simulation 2 - a_PF is changing parameter and time step is 600 (model_parameters_2) ####

library(deSolve)
library(ggplot2)
library(rootSolve)
library(mgcv)
library(dplyr)
library(tibbletime)
library(beepr)

source("helper_functions.R")
source("model_parameters_2_ts60.R")

#set seed to ensure reproducibility

set.seed(10)

#set parameters for changing observation error and rate of change

sim2_parameter_ts60 = expand.grid(obs_error = seq(0.00, 0.3, length.out = 7),
                                  rate_of_change = seq(-0.0002, -0.0011, length.out = 10),
                                  replicate = 1:100)

sim2_parameter_ts60$replicate <- 1:nrow(sim2_parameter_ts60)

n_sim2_ts60 = nrow(sim2_parameter_ts60)

for(i in seq_len(length.out = n_sim2_ts60)){
  
  #run the simulation
  
  sim_current_ts60 <- run_simulation(parameters = sim2_parameter_ts60[i,],
                                     time_seq = time_series)
  
  #fit population density of the simulation using GAM
  
  sim_predictions_ts60 <- extract_gam_predictions(parameters = sim2_parameter_ts60[i,], 
                                                  sim = sim_current_ts60,
                                                  time_seq = time_series)
  
  #calculate fisher information of the system
  
  sim_fisher_current_ts60 <- calc_fisher_current(parameters = sim2_parameter_ts60[i,],
                                                 predictions = sim_predictions_ts60,
                                                 time_seq = time_series)
  
  #calculate fisher information of the system using log transform of the s
  #ystem population densities
  
  sim_log_fisher_current_ts60 <- calc_fisher_current(parameters = sim2_parameter_ts60[i,],
                                                     predictions = log(sim_predictions_ts60),
                                                     time_seq = time_series)
  
  #find the min fisher information timepoint for when fisher information d
  #etects the regime shift
  
  sim_fisher_min_ts60 <- (which.min(sim_fisher_current_ts60)*10)
  sim2_parameter_ts60[i, "sim_fisher_min_ts60"] = sim_fisher_min_ts60
  
  #find the min fisher information timepoint for when fisher information 
  #of the log transform data detects the regime shift
  
  sim_log_fisher_min_ts60 <- (which.min(sim_log_fisher_current_ts60)*10)
  sim2_parameter_ts60[i, "sim_log_fisher_min_ts60"] = sim_log_fisher_min_ts60
  
  #find the regime shift of the simulation
  
  regime_shift_time_ts60 <- calc_regime_shift(parameters = sim2_parameter_ts60[i,],
                                              times = time_series)
  sim2_parameter_ts60[i, "regime_shift_time_ts60"] = regime_shift_time_ts60
  
}

beep(sound = 8)

#table for simulation data collected

View(sim2_parameter_ts60)

#set up the summary of the data so that the mean for fisher information,
#log fisher information and regime shift are calculated from the replicates
#of the simulations

sim2_parameter_summary_ts60 = sim2_parameter_ts60 %>%
  group_by(obs_error, rate_of_change)%>%
  summarize(fisher_min_mean_ts60 = mean(sim_fisher_min_ts60),
            log_fisher_min_mean_ts60 = mean(sim_log_fisher_min_ts60),
            regime_shift_time_mean_ts60 = mean(regime_shift_time_ts60),
            fisher_diff_mean_ts60 = mean(regime_shift_time_ts60-sim_fisher_min_ts60),
            fisher_log_diff_mean_ts60 = mean(regime_shift_time_ts60-sim_log_fisher_min_ts60))

View(sim2_parameter_summary_ts60)

#plot the difference between fisher information and the regime shift

ggplot(sim2_parameter_summary_ts60, aes(y = obs_error, x = rate_of_change, 
                                        fill = fisher_diff_mean_ts60)) +
  ggtitle("Fisher Information Difference, Time step of 60")+
  xlab("rate of change of parameter a_PF")+
  ylab("observation error")+
  geom_tile()+
  scale_fill_gradient2()

#plot the difference between the log fisher information and the regime shift

ggplot(sim2_parameter_summary_ts60, aes(y = obs_error, x = rate_of_change,
                                        fill = fisher_log_diff_mean_ts60)) +
  ggtitle("Log Fisher Information Difference, Time step of 60")+
  xlab("rate of change of parameter a_PF")+
  ylab("observation error")+
  geom_tile()+
  scale_fill_gradient2()

beep(sound = 8)