This repo contains essential documents used for gonorrhoea vaccination programme modelling project in Singapore setting. Please see files 
  (A) with "script" initials for analysis scripts 
  (B) with "graph/table" initials for visualization/summary generation.

The pipline is simply built as in: 
  (1)data - (2)calibration - (3)parameter collection - (4)projection - (5) visualization
  
(1) data: "MSM.csv" contains the Singapore gonorrhoea incidence over years for main scenario (MSM).
(2) calibration: "script_calibration_sgmsm.R" is the main script for running through a Bayesian inference clibration procedure; "script_calibration_support.R" and "script_calibration_plot_support.R" are the supporting functions along with the calibration procedure.
(3) parameter collection: "script_collection.R" is the main script for collecting parameter posteriors in the form that is easy for later projection.
(4) projection: "script_projection.R" is the main script for projecting transmission trajectory given collected parameters and specified vaccine profile, and it is called by result generating scripts; "script_projection_support.R" contains helper functions.
(5) visualization: "graph_impact_10-20.R", "graph_bivariate_10.R", "graph_trivariates_10.R" are the vaccination programme imact generating scripts; "graph_si_fitting.R", "graph_si_mcmc.R", "support.R" are calibration visualization generating scripts.
