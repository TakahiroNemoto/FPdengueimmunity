library(boot)
library(stats)
library(Rcpp)
#library(devtools)

rm(list = ls());

## move to the main folder
setwd('.../FPdengueimmunity')
getwd(); 

sourceCpp("./Code/sero_mcmc.cpp")
source("./Code/Functions.R")

#### test MCMC using synthetic data 
Test_consistency <- 0 # 1 -> yes, 0 -> no ... if yes, true parameters are chosen according to a inference result (for Fig.6). 
                      # In this case, it's required to do FigureProduction.R first to create a file of inferred parameters (Parameters_y_50_.txt) in ./Figure/With_data/


Islands_name = "_IDV" #
cases_dist <- as.matrix(t(read.table(paste('./Data/cases_dist',Islands_name,sep = "")))) # case distribution (matrix) for each time period and for each age
cases_dist_age <- as.vector(read.table('./Data/cases_dist_x'))
cases_dist_year <- as.vector(read.table('./Data/cases_dist_y')) # it should have at longest one year interval: less than one year is fine.
pop_dist <- as.matrix(t(read.table(paste('./Data/pop_dist',Islands_name,sep = "")))) # population distribution (matrix) for each time period and for each age
pop_dist_age <- as.vector(read.table('./Data/pop_dist_x'))
pop_dist_year <- as.vector(read.table('./Data/pop_dist_y')) # it should have at longest one year interval: less than one year is fine.

######### In the case data, age information is sometimes not provided. Age_reporting is the proportion of the data where the age information is provided.
######### (This will be set to 1 when NAN and 0.01 when 0 for a technical reason.)
Age_reporting <- t(read.table(paste('./Data/Age_report_rate',Islands_name,sep = "")))[,]
for (i in c(1:length(Age_reporting))){
  if (is.nan(Age_reporting[i])==TRUE){
    Age_reporting[i] <- 1
  } else if (Age_reporting[i]==0){
    Age_reporting[i] <- 0.01
  }
}
Age_reporting <- c(Age_reporting,1,1,1,1) # Just for 2015 - 2018 (the year between 2015 and 2018 won't be used  in the final version.)

#### Circulating serotype for each time interval
serotype_number <- t(read.table(paste('./Data/Serotype_number',Islands_name,sep = "")))[,]
serotype_number <- c(serotype_number,0,0,0,0) # for 2015 - 2018 (it won't be used in the final version.)

roof_age <- cases_dist_age[,1]
Surveilanceyear <- c(cases_dist_year[,][1:(length(cases_dist_year[,])-1)],2015,2016,2017,2018)
Surveilanceyear_roof <- c(cases_dist_year[,][1:(length(cases_dist_year[,])-1)],2015,2016,2017,2018,2019)

AA <- (length(roof_age) - 1)
TT <- length(Surveilanceyear) 
L <- max(roof_age)-1
Year_ <- c((min(Surveilanceyear)-max(roof_age)+1):((min(Surveilanceyear)-1)),Surveilanceyear)


################ This is a preparation for cross protections (see the main article for the cross protection) ############################ 
dummy_ <- Preparation_crossprotection(Surveilanceyear,serotype_number)
Surveilanceyear_extended <- dummy_[[1]] # to calculate the cross protections, the surveilance year is defined from 1973
serotype_number_extended <- dummy_[[2]] # to calculate the cross protections, serotype is defined from 1973
AddBegining <- dummy_[[3]] # 6 = 1979 - 1973
Change_serotype <- dummy_[[4]] # if there is a different serotype within two years, the corresponding element takes 1 otherwise 0.
How_many_within_kkyears <- dummy_[[5]]  # how many time intervals within two years. 
Time_interval_set <- dummy_[[6]]  # the lengths of each time interval within two years.


############# prepare the seroprevalence data 
Seroposi_2014 = read.table(file='./Data/Seroposi_2014')[,]
TrialNum_2014 = read.table(file='./Data/TrialNum_2014')[,]
Seroposi_1_2014 = read.table(file='./Data/Seroposi_1_2014')[,]
Seroposi_2_2014 = read.table(file='./Data/Seroposi_2_2014')[,]
Seroposi_3_2014 = read.table(file='./Data/Seroposi_3_2014')[,]
Seroposi_4_2014 = read.table(file='./Data/Seroposi_4_2014')[,]
Seroposi_twice_2014 = read.table(file='./Data/Seroposi_twice_2014')[,]
#
Seroposi_2015 = read.table(file='./Data/Seroposi_2015')[,]
Seroposi_1_2015 = read.table(file='./Data/Seroposi_1_2015')[,]
Seroposi_2_2015 = read.table(file='./Data/Seroposi_2_2015')[,]
Seroposi_3_2015 = read.table(file='./Data/Seroposi_3_2015')[,]
Seroposi_4_2015 = read.table(file='./Data/Seroposi_4_2015')[,]
Seroposi_twice_2015 = read.table(file='./Data/Seroposi_twice_2015')[,]
TrialNum_2015 = read.table(file='./Data/TrialNum_2015')[,]


############# initial condition for the parameters

# FOI
# initial condition for FOI from 1979 to 2018. The last four components represent FOI in the four epidemics before 1979 (1944, 1964-1969, 1971-1975, 1976-1978). (We set FOI is constant over the same epidemic before 1978.)
dummy_lambda <-  c(rep(0.03,length(Year_[which(Year_==1979):which(Year_==2018)])),0.03,0.03,0.03,0.03)

# Reporting rate
Primary_PHI <- rep(1.,4) ## Relative reporting rate for primary infections. Using the notation in the article, Primary_PHI[2]: phi(2,1)/phi(1,1), Primary_PHI[3]: phi(3,1)/phi(1,1), Primary_PHI[4]: phi(4,1)/phi(1,1)
Primary_PHI[1] <- 0.21   ## Primary_PHI[1] is the largest value of Background_PHI (see below)
Secondary_PHI <- rep(1.,12) ## Relative reporting rate for secondary infections. Using the notation in the article, Secondary_PHI[1]: phi(1,2)/phi(1,1), Secondary_PHI[4]: phi(2,2)/phi(1,1), Secondary_PHI[7]: phi(3,2)/phi(1,1), Secondary_PHI[10]: phi(4,2)/phi(1,1). The other components of Secondary_PHI are not used in the final version.
Background_PHI <- c(5.,0.02,0.03,0.05,0.03,0.02,0.03, # parameters that defines the background phi T(t) in the article. They are "steep, ca, cb, cc, ya, yb, yc".
                    # T(t) = ca*(tanh(steep*(Surveilanceyear - 1985)) + 1.)/2. + cb*(tanh(steep*(Surveilanceyear - 1990)) + 1.)/2. + cc*(tanh(steep*(Surveilanceyear - 1995)) + 1.)/2. - ca - cb - cc
                    #       + ya*(tanh(steep*(Surveilanceyear - 2000)) + 1.)/2. + yb*(tanh(steep*(Surveilanceyear - 2004)) + 1.)/2. + yc*(tanh(steep*(Surveilanceyear - 2008)) + 1.)/2. + (Primary_PHI[1] - ya - yb - yc)
                    1., # parameter for the negative binomial distribution (k in the main article.)
                    1., # parameter for age factor. A(1-4) in the notation of the article
                    1., # parameter for age factor. A(10-14) in the notation of the article
                    1.  # parameter for age factor. A(15+) in the notation of the article
) 
dummy_PHI <- c(Background_PHI,Primary_PHI,Secondary_PHI)

# New parameter 
b_miss <- 0.1 # This parameter is not used in the final version. (Please ignore it.)
dummy_new <- c(b_miss)

# total initial parameters
xtest <- c(dummy_lambda,dummy_PHI,dummy_new)

lambda_length <- length(dummy_lambda)


## This is when the seroprevalence study was conducted (measured from the begining of the year.). 
Seroposi_timeshiftpercent_2014 <- 0.415 ## meaning that the survey is conducted at 2014.415
Seroposi_timeshiftpercent_2015 <- 0.749 ## meaning that the survey is conducted at 2015.749
How_many_points_for_integral_approximation <- 4 # The step number of the Euler integral to evaluate the integrals.


## This is (for Fig.6) to perform the inference based on synthetic data.
if (Test_consistency==1) {
  xtest_test <- read.table("./Figure/With_data/Parameters_y_50_.txt",sep="")[,]
  xtest_test[!is.finite(xtest_test)] <- 1.
  
  ## perform a simulation based on the parameter (xtest_test) to create a synthetic data
  cases_dist <- matrix(Simulation_casedata(xtest_test), nrow = AA, ncol = length(Simulation_casedata(xtest_test))/AA)
  
  ## perform a simulation based on the parameter (xtest_test) to generate seroprevalence data  
  dummy_ <- Simulation_serodata(xtest_test)
  Seroposi_2014[1:AA] <- dummy_[[1]]
  Seroposi_2015[1:AA] <- dummy_[[2]]
}



############## MCMC inference below
nbIteration<-450000 ## Number of iteration
burnIn<-30000 ## Initial relaxation period
#nbIteration<-40000
#burnIn<-20000
nbParam<-length(xtest)
param<-xtest
sdProposal<-rep(0.01,nbParam) #  MC step size (this will be adjusted during the inference)


calculateLL <- function(x){
  x <- abs(x)
  return(Log_Likelihood_casedata(x))
} 

logLik<-rep(0,nbIteration) # in which log likelihood is stored.
accept<-matrix(0,ncol=nbParam,nrow=nbIteration) # in which how many accepted is stored. Based on this, MC step size is modified.
parameters<-matrix(0,ncol=nbParam,nrow=nbIteration) # in which parameters are stored.

iteration<-1
parameters[iteration,]<-param

Changefreq <- 500 ## each Changefreq MC steps, sdProposal is modified during initial relaxation stage.
logLik[iteration]<-calculateLL(parameters[iteration,])
for(iteration in 2:nbIteration)
{ 
  if (iteration %% 50==0){print(parameters[iteration-1,(lambda_length+1):(lambda_length+length(dummy_PHI))])} # each 50 step, show some parameters.
  if (iteration %% 50==0) {print(c(iteration,logLik[iteration-1],calculateLL(parameters[iteration-1,])))}  # each 50 step, show loglikelihood
  if ((iteration %% Changefreq == 0)){ 
    if ((iteration < burnIn)){sdProposal<-sdProposal*exp((apply(accept[(iteration+1-Changefreq):(iteration-1),],2,mean)-0.24))} # change sdProposal so that the acceptance rate will become 0.24
    print(apply(accept[(iteration+1-Changefreq):(iteration-1),],2,mean))
  }
  parameters[iteration,]<-parameters[iteration-1,]
  logLik[iteration]<-logLik[iteration-1]
  for(iParam in 1:nbParam)
  {oldParam<-parameters[iteration,iParam]
   dummytoremovedivergence <- sdProposal[iParam]*rnorm(1)
   newParam<-oldParam*exp(dummytoremovedivergence)
  
   parameters[iteration,iParam]<-newParam
   newLogLik<-calculateLL(parameters[iteration,])

   if(log(runif(1))<newLogLik-logLik[iteration]+dummytoremovedivergence){ 
     logLik[iteration]<-newLogLik # accept the change
     accept[iteration,iParam]<-1
   }
   else{
     parameters[iteration,iParam]<-oldParam # reject the change
   }
 }
}

# data save
if (Test_consistency==0){
  save.image(file = "./Result/With_data/my_work_space.RData")
} else if (Test_consistency==1){
  save.image(file = "./Result/With_synthetic_data/my_work_space.RData")
}

