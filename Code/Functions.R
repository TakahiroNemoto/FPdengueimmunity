
Preparation_crossprotection <- function(Surveilanceyear,
                                        serotype_number){
  ################ Below is a preparation for cross protections (see the main article for the cross protection) ############################ 
  kk <- 2 ## this is the duration of the year the cross protection is taken into account.
  Surveilanceyear_extended <- c(1973,1974,1975,1976,1977,1978,Surveilanceyear)
  serotype_number_extended <- c(2,2,2,1,1,1,serotype_number)
  AddBegining <- length(Surveilanceyear_extended)-length(Surveilanceyear)
  Change_serotype <- Surveilanceyear
  How_many_within_kkyears <- Surveilanceyear
  # Surveilanceyear_extended : to calculate the cross protections, the surveillance year is defined from 1973
  # serotype_number_extended : to calculate the cross protections, serotype is defined from 1973
  # AddBegining  : 6 = 1979 - 1973
  # Change_serotype  : if there is a different serotype within two years, the corresponding element takes 1 otherwise 0.
  # How_many_within_kkyears : how many time intervals within two years. 
  
  
  for (i in c(1:length(Surveilanceyear_extended))){
    if (Surveilanceyear_extended[i] >= min(Surveilanceyear)){
      The_different_sero <- 0
      out_ <- 0
      j <- 0
      while (out_ == 0){
        j <- j+1
        if (serotype_number_extended[i] != serotype_number_extended[i-j]){
          The_different_sero <- 1
        }
        if (Surveilanceyear_extended[i-j] <= Surveilanceyear_extended[i] - kk){
          out_ <- 1
        }
      }
      Change_serotype[i-AddBegining] <- The_different_sero
      How_many_within_kkyears[i-AddBegining] <- j
    }
#    print(c(Surveilanceyear_extended[i],serotype_number_extended[i]))
  }
  
  Time_interval_set <- matrix(0,nrow=length(Surveilanceyear),ncol=max(How_many_within_kkyears))
  # Time_interval_set : the lengths of each time interval within two years.
  
  for (i in c(1:length(Surveilanceyear_extended))){
    if (Surveilanceyear_extended[i] >= min(Surveilanceyear)){
      The_different_sero <- 0
      out_ <- 0
      j <- 0
      while (out_ == 0){
        j <- j+1
        Time_interval_set[i-AddBegining,j] <- Surveilanceyear_extended[i-j+1] - Surveilanceyear_extended[i-j]
        if (serotype_number_extended[i] != serotype_number_extended[i-j]){
          The_different_sero <- 1
        }
        if (Surveilanceyear_extended[i-j] <= Surveilanceyear_extended[i] - kk){
          out_ <- 1
          Time_interval_set[i-AddBegining,j] <- Surveilanceyear_extended[i-j+1] - (Surveilanceyear_extended[i] - kk)
        }
      }
    }
#    print(c(Surveilanceyear_extended[i],serotype_number_extended[i]))
  }
  
  # # to check
  # for (i in c(1:length(Surveilanceyear))){
  #   print(c(Surveilanceyear[i],Surveilanceyear[i]-2,serotype_number[i],Change_serotype[i],How_many_within_kkyears[i]))
  #   print(Time_interval_set[i,])
  # }
  
  Change_serotype[(length(Surveilanceyear)-3):length(Surveilanceyear)] <- 0
  
  return(list(Surveilanceyear_extended, serotype_number_extended, AddBegining, Change_serotype, How_many_within_kkyears,Time_interval_set))
  
}


#################################################################
#################################################################
#################################################################
#### Calculate loglikelihood function for a given parameter x.
Log_Likelihood_casedata <- function(x){
  
  ########  get FOI for each serotype from the parameter x 
  dummy_sub_lambda <- Substitution_lambda(x,serotype_number,Year_,lambda_length,
                                          TT, length(Year_), which(Year_==1944),
                                          which(Year_==1978))
  lambda1 <- dummy_sub_lambda[1,]
  lambda2 <- dummy_sub_lambda[2,]
  lambda3 <- dummy_sub_lambda[3,]
  lambda4 <- dummy_sub_lambda[4,]
  
  ########  get the integral of lambda from each lambda
  integrate_lambda1 <- integrate_lambda_func_Integral(lambda1,Year_,length(Surveilanceyear), L) 
  integrate_lambda2 <- integrate_lambda_func_Integral(lambda2,Year_,length(Surveilanceyear), L) 
  integrate_lambda3 <- integrate_lambda_func_Integral(lambda3,Year_,length(Surveilanceyear), L) 
  integrate_lambda4 <- integrate_lambda_func_Integral(lambda4,Year_,length(Surveilanceyear), L) 
  
  ########  compute x, y_i (in the notation of the article) 
  x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4))  # x
  mono1 <- x0 * (exp(integrate_lambda1) - 1.) # y1
  mono2 <- x0 * (exp(integrate_lambda2) - 1.) # y2
  mono3 <- x0 * (exp(integrate_lambda3) - 1.) # y3
  mono4 <- x0 * (exp(integrate_lambda4) - 1.) # y4
  
  ##### these parameters are used to define the time factor T(t) for reporting rates
  steep <- 5.
  ca <- x[(lambda_length+2)]
  cb <- x[(lambda_length+3)]
  cc <- x[(lambda_length+4)]
  ya <- x[(lambda_length+5)]
  yb <- x[(lambda_length+6)]
  yc <- x[(lambda_length+7)]
  negbin_disp <- x[(lambda_length+8)] # parameter k in the negative binomial distribution 
  age_factor1 <- x[(lambda_length+9)] # age factor A(1-4)
  age_factor2 <- x[(lambda_length+10)] # age factor A(5-9)
  age_factor3 <- x[(lambda_length+11)] # age factor A(10-14)
  
  # primary infection reporting rates: phi(i,1)/phi(1,1) using the notation of the article
  PHI_Primary <- x[(lambda_length+length(Background_PHI)+1):(lambda_length+length(Background_PHI)+length(Primary_PHI))]
  PHI_Primary[1] <- 1. # this is phi(1,1)/phi(1,1)

  # secondary infection reporting rates:  phi(i,2)/phi(1,1) using the notation of the article
  PHI_Secondary <- x[(lambda_length+length(Background_PHI)+length(Primary_PHI)+1):(lambda_length+length(Background_PHI)+length(Primary_PHI)+length(Secondary_PHI))]
  PHI_Secondary[1:3] = PHI_Secondary[1] # phi(1,2)/phi(1,1)
  PHI_Secondary[4:6] = PHI_Secondary[4] # phi(2,2)/phi(1,1)
  PHI_Secondary[7:9] = PHI_Secondary[7] # phi(3,2)/phi(1,1)
  PHI_Secondary[10:12] = PHI_Secondary[10] # phi(4,2)/phi(1,1)
  
  largest_PHI_background <- x[lambda_length+length(Background_PHI)+1]
  Notaccepted <- 0
  # PHI_background is the time factor T(t) defined using tanh functions
  PHI_background <- ca*(tanh(steep*(Surveilanceyear - 1985)) + 1.)/2. + cb*(tanh(steep*(Surveilanceyear - 1990)) + 1.)/2. + cc*(tanh(steep*(Surveilanceyear - 1995)) + 1.)/2. - ca - cb - cc
  PHI_background <- PHI_background + ya*(tanh(steep*(Surveilanceyear - 2000)) + 1.)/2. + yb*(tanh(steep*(Surveilanceyear - 2004)) + 1.)/2. + yc*(tanh(steep*(Surveilanceyear - 2008)) + 1.)/2. + (largest_PHI_background - ya - yb - yc)
  if (ca + cb + cc + ya + yb + yc > largest_PHI_background){
    Notaccepted <- 1
  } else if (steep>10){
    Notaccepted <- 1
  } else if (negbin_disp > 3){
    Notaccepted <- 1
  }
  
  ## calculate the susceptible to secondary infections by taking into account the cross immunity (that last two years). 
  ## mono_dummy[[i]]: y_i - del x_i (i=1,2,3,4) See the eq (5) in the article.
  mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                     which(Change_serotype==1),
                                     How_many_within_kkyears,
                                     dummy_sub_lambda,
                                     Time_interval_set,
                                     serotype_number_extended,
                                     AddBegining)  
  
  #  Perform the Euler method to compute the expected value in the negative binomial distribution eq.(7) (without the population term P(t,a)) in the article
  #  with a step number "How_many_points_for_integral_approximation" In the article this parameter is 4.
  if (How_many_points_for_integral_approximation == 1){
    poi_rate <- poi_rate_func(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                              PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                              TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3)
  } else if (How_many_points_for_integral_approximation == 2){
    poi_rate_1 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.5)
    
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.5) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_2 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.5)
    poi_rate <- poi_rate_1 + poi_rate_2
  } else if (How_many_points_for_integral_approximation == 4){
    poi_rate_1 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    #############    
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.25) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.25) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.25) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.25) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_2 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    #############
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.5) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_3 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    ############    
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.75) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.75) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.75) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.75) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_4 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    ############    
    
    poi_rate <- poi_rate_1 + poi_rate_2 + poi_rate_3 + poi_rate_4
  }     

  #  Multiplying the population term this becomes eq.(7)in the article
  poi_rate_cohort <- poi_rate_cohort_func(poi_rate[,1:which(Surveilanceyear==2014)],roof_age,pop_dist,which(Surveilanceyear == 2014),AA)
  PoissonRateVec <- as.vector(poi_rate_cohort)
  
  
  # Computing the log likelihood for the case data corresponding to the first term in eq.(10)
  Resul_temp <- sum(dnbinom(as.vector(cases_dist), size=(PoissonRateVec)**negbin_disp, mu=PoissonRateVec,log=TRUE))
  
  
  # below is to compute the log likelihood for seroprevalence survey
  # This is to shift x0 to the year 2014.42 and 2015.75 (which are the time at which the survey was conducted). 
  integrate_lambda_shifted_July1_2014 <- integrate_lambda_func_Integral_SHIFTED(lambda1+lambda2+lambda3+lambda4,Year_,length(Surveilanceyear), L,Seroposi_timeshiftpercent_2014)[,which(Surveilanceyear == 2014)] 
  integrate_lambda_shifted_Oct1_2015 <- integrate_lambda_func_Integral_SHIFTED(lambda1+lambda2+lambda3+lambda4,Year_,length(Surveilanceyear), L,Seroposi_timeshiftpercent_2015)[,which(Surveilanceyear == 2015)]
  x0_shifted_2014 <- exp(-integrate_lambda_shifted_July1_2014) 
  x0_shifted_2015 <- exp(-integrate_lambda_shifted_Oct1_2015) 
  
  
  # This parameter is not used finally.
  b_miss <- x[length(dummy_lambda) + length(dummy_PHI) + 1]
  
  # Computing the log likelihood for the seroprevalence data corresponding to the second term in eq.(10)
  Sero_2014_fit = Seroposi_2014[6:16]
  DataSero_2014_fit <- 1-x0_shifted_2014[6:16] 
  Resul_temp <- Resul_temp + sum(dbinom(Sero_2014_fit,TrialNum_2014[6:16],DataSero_2014_fit,log=TRUE))
  Sero_2015_fit = Seroposi_2015[12:70]
  DataSero_2015_fit <- 1-x0_shifted_2015[12:70]
  Resul_temp <- Resul_temp + sum(dbinom(Sero_2015_fit,TrialNum_2015[12:70],DataSero_2015_fit,log=TRUE))
  
  if (is.na(Resul_temp)){Resul_temp <- -Inf}
  
  
  if (Notaccepted == 1){
    Resul_temp <- -Inf
  } else if (max(lambda1) > 1){
    Resul_temp <- -Inf
  } else if (max(lambda2) > 1){
    Resul_temp <- -Inf
  } else if (max(lambda3) > 1){
    Resul_temp <- -Inf
  } else if (max(lambda4) > 1){
    Resul_temp <- -Inf
  } else if (b_miss < 0){
    Resul_temp <- -Inf
  } else if (b_miss > 1){
    Resul_temp <- -Inf
  } else if (max(max(PHI_Secondary)*PHI_background) > 1.){
    Resul_temp <- -Inf
  } else if (max(max(PHI_Primary)*PHI_background) > 1.){
    Resul_temp <- -Inf
  }
  
  return(Resul_temp)
}



#################################################################
#################################################################
#################################################################
### simulate the case number from the parameter x
Simulation_casedata <- function(x){
  
  ########  get FOI for each serotype from the parameter x 
  dummy_sub_lambda <- Substitution_lambda(x,serotype_number,Year_,lambda_length,
                                          TT, length(Year_), which(Year_==1944),
                                          which(Year_==1978))
  lambda1 <- dummy_sub_lambda[1,]
  lambda2 <- dummy_sub_lambda[2,]
  lambda3 <- dummy_sub_lambda[3,]
  lambda4 <- dummy_sub_lambda[4,]
  
  ########  get the integral of lambda from each lambda
  integrate_lambda1 <- integrate_lambda_func_Integral(lambda1,Year_,length(Surveilanceyear), L) 
  integrate_lambda2 <- integrate_lambda_func_Integral(lambda2,Year_,length(Surveilanceyear), L) 
  integrate_lambda3 <- integrate_lambda_func_Integral(lambda3,Year_,length(Surveilanceyear), L) 
  integrate_lambda4 <- integrate_lambda_func_Integral(lambda4,Year_,length(Surveilanceyear), L) 
  
  ########  compute x, y_i (in the notation of the article) 
  x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
  mono1 <- x0 * (exp(integrate_lambda1) - 1.)
  mono2 <- x0 * (exp(integrate_lambda2) - 1.)
  mono3 <- x0 * (exp(integrate_lambda3) - 1.)
  mono4 <- x0 * (exp(integrate_lambda4) - 1.)
  
  ##### these parameters are used to define the time factor T(t) for reporting rates
  steep <- 5.
  ca <- x[(lambda_length+2)]
  cb <- x[(lambda_length+3)]
  cc <- x[(lambda_length+4)]
  ya <- x[(lambda_length+5)]
  yb <- x[(lambda_length+6)]
  yc <- x[(lambda_length+7)]
  negbin_disp <- x[(lambda_length+8)]
  age_factor1 <- x[(lambda_length+9)]
  age_factor2 <- x[(lambda_length+10)]
  age_factor3 <- x[(lambda_length+11)]
  largest_PHI_background <- x[lambda_length+length(Background_PHI)+1]
  
  # PHI_background is the time factor T(t) defined using tanh functions
  PHI_background <- ca*(tanh(steep*(Surveilanceyear - 1985)) + 1.)/2. + cb*(tanh(steep*(Surveilanceyear - 1990)) + 1.)/2. + cc*(tanh(steep*(Surveilanceyear - 1995)) + 1.)/2. - ca - cb - cc
  PHI_background <- PHI_background + ya*(tanh(steep*(Surveilanceyear - 2000)) + 1.)/2. + yb*(tanh(steep*(Surveilanceyear - 2004)) + 1.)/2. + yc*(tanh(steep*(Surveilanceyear - 2008)) + 1.)/2. + (largest_PHI_background - ya - yb - yc)
  # primary infection reporting rates: phi(i,1)/phi(1,1) using the notation of the article
  PHI_Primary <- x[(lambda_length+length(Background_PHI)+1):(lambda_length+length(Background_PHI)+length(Primary_PHI))]
  PHI_Primary[1] <- 1. # this is phi(1,1)/phi(1,1)
  # secondary infection reporting rates:  phi(i,2)/phi(1,1) using the notation of the article
  PHI_Secondary <- x[(lambda_length+length(Background_PHI)+length(Primary_PHI)+1):(lambda_length+length(Background_PHI)+length(Primary_PHI)+length(Secondary_PHI))]
  PHI_Secondary[1:3] = PHI_Secondary[1]
  PHI_Secondary[4:6] = PHI_Secondary[4]
  PHI_Secondary[7:9] = PHI_Secondary[7]
  PHI_Secondary[10:12] = PHI_Secondary[10]
  
  ## calculate the susceptible to secondary infections by taking into account the cross immunity (that last two years). 
  ## mono_dummy[[i]]: y_i - del x_i (i=1,2,3,4) See the eq (5) in the article.
  mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                     which(Change_serotype==1),
                                     How_many_within_kkyears,
                                     dummy_sub_lambda,
                                     Time_interval_set,
                                     serotype_number_extended,
                                     AddBegining)  
  
  
  #  Perform the Euler method to compute the expected value in the negative binomial distribution eq.(7) (without the population term P(t,a)) in the article
  #  with a step number "How_many_points_for_integral_approximation" In the article this parameter is 4.
  if (How_many_points_for_integral_approximation == 1){
    poi_rate <- poi_rate_func(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                              PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                              TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3)
  } else if (How_many_points_for_integral_approximation == 2){
    poi_rate_1 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.5)
    
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.5) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_2 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.5)
    poi_rate <- poi_rate_1 + poi_rate_2
  } else if (How_many_points_for_integral_approximation == 4){
    poi_rate_1 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    #############    
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.25) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.25) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.25) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.25) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_2 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    #############
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.5) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.5) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_3 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    ############    
    integrate_lambda1 <- integrate_lambda_func_Integral_SHIFTED(lambda1,Year_,length(Surveilanceyear), L, 0.75) 
    integrate_lambda2 <- integrate_lambda_func_Integral_SHIFTED(lambda2,Year_,length(Surveilanceyear), L, 0.75) 
    integrate_lambda3 <- integrate_lambda_func_Integral_SHIFTED(lambda3,Year_,length(Surveilanceyear), L, 0.75) 
    integrate_lambda4 <- integrate_lambda_func_Integral_SHIFTED(lambda4,Year_,length(Surveilanceyear), L, 0.75) 
    
    x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
    
    mono1 <- x0 * (exp(integrate_lambda1) - 1.)
    mono2 <- x0 * (exp(integrate_lambda2) - 1.)
    mono3 <- x0 * (exp(integrate_lambda3) - 1.)
    mono4 <- x0 * (exp(integrate_lambda4) - 1.)  
    mono_dummy <- mono_tilde_calculate(L,TT,x0,mono1,mono2,mono3,mono4,
                                       which(Change_serotype==1),
                                       How_many_within_kkyears,
                                       dummy_sub_lambda,
                                       Time_interval_set,
                                       serotype_number_extended,
                                       AddBegining)  
    
    poi_rate_4 <- poi_rate_func_with_rate(lambda1,lambda2,lambda3,lambda4,x0,mono_dummy[[1]],mono_dummy[[2]],mono_dummy[[3]],mono_dummy[[4]],
                                          PHI_background*Age_reporting,PHI_Primary,PHI_Secondary,serotype_number,
                                          TT, L,Surveilanceyear_roof,age_factor1,age_factor2,age_factor3,0.25)
    ############    
    
    poi_rate <- poi_rate_1 + poi_rate_2 + poi_rate_3 + poi_rate_4
  }     
  
  #  Multiplying the population term this becomes eq.(7)in the article
  poi_rate_cohort <- poi_rate_cohort_func(poi_rate[,1:which(Surveilanceyear==2014)],roof_age,pop_dist,which(Surveilanceyear == 2014),AA)
  PoissonRateVec <- as.vector(poi_rate_cohort)

  # simulation  
  Resul_temp <- rnbinom(as.vector(cases_dist), size=(PoissonRateVec)**negbin_disp, mu=PoissonRateVec)
  return(Resul_temp)
}




#################################################################
#################################################################
#################################################################
### compute x0, mono1, ..., mono4
x0_mono1_mono2_mono3_mono4 <- function(x){
  
  ########  get FOI for each serotype from the parameter x 
  dummy_sub_lambda <- Substitution_lambda(x,serotype_number,Year_,lambda_length,
                                          TT, length(Year_), which(Year_==1944),
                                          which(Year_==1978))
  lambda1 <- dummy_sub_lambda[1,]
  lambda2 <- dummy_sub_lambda[2,]
  lambda3 <- dummy_sub_lambda[3,]
  lambda4 <- dummy_sub_lambda[4,]
  
  ########  get the integral of lambda from each lambda
  integrate_lambda1 <- integrate_lambda_func_Integral(lambda1,Year_,length(Surveilanceyear), L) 
  integrate_lambda2 <- integrate_lambda_func_Integral(lambda2,Year_,length(Surveilanceyear), L) 
  integrate_lambda3 <- integrate_lambda_func_Integral(lambda3,Year_,length(Surveilanceyear), L) 
  integrate_lambda4 <- integrate_lambda_func_Integral(lambda4,Year_,length(Surveilanceyear), L) 
  
  ########  compute x, y_i (in the notation of the article) 
  x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
  mono1 <- x0 * (exp(integrate_lambda1) - 1.)
  mono2 <- x0 * (exp(integrate_lambda2) - 1.)
  mono3 <- x0 * (exp(integrate_lambda3) - 1.)
  mono4 <- x0 * (exp(integrate_lambda4) - 1.)
  
  return(c(as.vector(x0),as.vector(mono1),as.vector(mono2),as.vector(mono3),as.vector(mono4)))
}


#################################################################
#################################################################
#################################################################
### compute seropositive
Serovositivity_casedata <- function(x){
  ########  get FOI for each serotype from the parameter x 
  dummy_sub_lambda <- Substitution_lambda(x,serotype_number,Year_,lambda_length,
                                          TT, length(Year_), which(Year_==1944),
                                          which(Year_==1978))
  lambda1 <- dummy_sub_lambda[1,]
  lambda2 <- dummy_sub_lambda[2,]
  lambda3 <- dummy_sub_lambda[3,]
  lambda4 <- dummy_sub_lambda[4,]
  
  ########  get the integral of lambda from each lambda
  integrate_lambda1 <- integrate_lambda_func_Integral(lambda1,Year_,length(Surveilanceyear), L) 
  integrate_lambda2 <- integrate_lambda_func_Integral(lambda2,Year_,length(Surveilanceyear), L) 
  integrate_lambda3 <- integrate_lambda_func_Integral(lambda3,Year_,length(Surveilanceyear), L) 
  integrate_lambda4 <- integrate_lambda_func_Integral(lambda4,Year_,length(Surveilanceyear), L) 
  
  ########  compute x, y_i (in the notation of the article) 
  x0 <- exp(-(integrate_lambda1 + integrate_lambda2 + integrate_lambda3 + integrate_lambda4)) 
  
  
  # This is to shift x0 to the year 2014.42 and 2015.75 (which are the time at which the survey was conducted). 
  integrate_lambda_shifted_July1_2014 <- integrate_lambda_func_Integral_SHIFTED(lambda1+lambda2+lambda3+lambda4,Year_,length(Surveilanceyear), L,Seroposi_timeshiftpercent_2014)[,which(Surveilanceyear == 2014)] 
  integrate_lambda_shifted_Oct1_2015 <- integrate_lambda_func_Integral_SHIFTED(lambda1+lambda2+lambda3+lambda4,Year_,length(Surveilanceyear), L,Seroposi_timeshiftpercent_2015)[,which(Surveilanceyear == 2015)]
  x0_shifted_2014 <- exp(-integrate_lambda_shifted_July1_2014) 
  x0_shifted_2015 <- exp(-integrate_lambda_shifted_Oct1_2015) 
  
  x0[,which(Surveilanceyear == 2014)] <- x0_shifted_2014
  x0[,which(Surveilanceyear == 2015)] <- x0_shifted_2015
  
  return(as.vector(x0))
}




#################################################################
#################################################################
#################################################################
### Simulate seroprevalence data
Simulation_serodata <- function(x){
  xtest_test <- x
  dummy_x0_monos_ <- x0_mono1_mono2_mono3_mono4(xtest_test)
  x0_dum_test <- matrix(dummy_x0_monos_[1:(L*TT)], nrow = L, ncol = TT)
  
  ########  get FOI for each serotype from the parameter x 
  dummy_sub_lambda <- Substitution_lambda(xtest_test,serotype_number,Year_,lambda_length,
                                          TT, length(Year_), which(Year_==1944),which(Year_==1978))
  lambda1 <- dummy_sub_lambda[1,]
  lambda2 <- dummy_sub_lambda[2,]
  lambda3 <- dummy_sub_lambda[3,]
  lambda4 <- dummy_sub_lambda[4,]
  
  # This is to shift x0 to the year 2014.42 and 2015.75 (which are the time at which the survey was conducted). 
  integrate_lambda_shifted_July1_2014 <- integrate_lambda_func_Integral_SHIFTED(lambda1+lambda2+lambda3+lambda4,Year_,length(Surveilanceyear), L,Seroposi_timeshiftpercent_2014)[,which(Surveilanceyear == 2014)] 
  integrate_lambda_shifted_Oct1_2015 <- integrate_lambda_func_Integral_SHIFTED(lambda1+lambda2+lambda3+lambda4,Year_,length(Surveilanceyear), L,Seroposi_timeshiftpercent_2015)[,which(Surveilanceyear == 2015)]
  x0_shifted_2014 <- exp(-integrate_lambda_shifted_July1_2014) 
  x0_shifted_2015 <- exp(-integrate_lambda_shifted_Oct1_2015) 
  
  x0_dum_test[,which(Surveilanceyear == 2014)] <- x0_shifted_2014
  x0_dum_test[,which(Surveilanceyear == 2015)] <- x0_shifted_2015
  
  sero_2014_len <- length(x0_dum_test[,which(Surveilanceyear == 2014)]) 
  sero_2015_len <- length(x0_dum_test[,which(Surveilanceyear == 2015)]) 
  return(list(rbinom(sero_2014_len,TrialNum_2014[1:sero_2014_len],1-x0_dum_test[,which(Surveilanceyear == 2014)]),
              rbinom(sero_2015_len,TrialNum_2015[1:sero_2015_len],1-x0_dum_test[,which(Surveilanceyear == 2015)])))
}
