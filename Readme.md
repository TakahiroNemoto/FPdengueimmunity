# Generating Fig.1 to Fig.6 in the article "Reconstructing long-term dengue virus immunity in French Polynesia" published from PLOS Neglected Tropical Diseases

# Data

In Data folder, there are three types of data: case data, population data, and seroprevalence survey data. 
- case data : cases_dist_IDV (a matrix storing the number of reported cases for each time period and for each age in IDV), cases_dist_x (a vector for defining each age group), cases_dist_y (a vector for defining each time period), Serotype_number_IDV (a vector indicating circulating serotypes for each year in IDV) and Age_report_rate_IDV (a vector storing the rate at which an age information is provided for a given patient)
- demographic data : pop_dist_IDV (a matrix storing the population for each time period and for each age), pop_dist_x (age groups) and pop_dist_y (time periods)
- seroprevalence survey data : Seroposi_T (a vector storing the number of all the seropositives for each age in the survey conducted in the year T), TrialNum_T (a vector storing the number of participants of the survey conducted in the year T for each age), Seroposi_i_T (a vector storing the number of seropositives for each age and for each serotype i in the survey conducted in the year T) and Seroposi_twice_T (a vector storing the number of seropositives against multiple serotypes for each age in the year T), where i = 1,2,3,4 and T=2014,2015 


# Code

In Code folder, 
- Simulation.R : it runs a MCMC inference. To do so, first fill in setwd() the path to the FPdengueimmunity folder. If Test_consistency == 0, it performs the inference using the FP data in Data folder (for Fig.1-Fig.5). The result appears in Result/With_data/ as a RData file. If Test_consistency == 1, it first generates synthetic data using the parameters of the inference result obtained with Test_consistency == 0 and then performs an inference with this synthetic data (for Fig.6). The result appears in Result/With_synthetic_data/ as a RData file. Before doing this second simulation, do Figure.R with which_data == 'FP' to generate "Parameters_y_50_.txt" in Figure/With_data folder. 
- Figure.R : after doing Simulation.R, executing this script will make figures used in the article. To do so, first fill in setwd() the path to the FPdengueimmunity folder. If which_data == 'FP', this makes figures based on the simulation result in Result/With_data/ (Fig.1-Fig.5). If which_data == 'Synthetic', this makes figures based on the simulation result in Result/With_synthetic_data/ (Fig.6).

# Figure
The figures appear in Figure/With_data/ when which_data == 'FP', while they appear in Figure/With_synthetic_data/ when which_data == 'Synthetic'.

# Result
The results of simulations appear in Result/With_data/ when Test_consistency == 0, while they appear in Result/With_synthetic_data/ when Test_consistency == 1.


# Version of R
R 3.6.0 : https://cran.r-project.org/src/base/R-3/R-3.6.0.tar.gz





