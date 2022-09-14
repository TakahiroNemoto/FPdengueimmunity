library(boot)
library(stats)
library(Rcpp)
library(devtools)
library(binom)
library("readxl")


rm(list = ls());

## move to the main folder
setwd('.../FPdengueimmunity')
getwd(); 

#### Choose the data used for simulations: 
which_data = 'FP' ## FP data
#which_data = 'Synthetic' ## synthetic data

if (which_data == 'FP') {
  data_get_place = "./Result/With_data"
  fig_save_place = "./Figure/With_data"
} else if (which_data == 'Synthetic') {
  data_get_place = "./Result/With_synthetic_data"
  fig_save_place = "./Figure/With_synthetic_data"
}


#### Upload the simulation result
load(file = paste(data_get_place,"/my_work_space.RData",sep=""))
sourceCpp("./Code/sero_mcmc.cpp")
source("./Code/Functions.R")


#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
## Fig.1 
## (1) Case number vs year, (2) Case number vs Age, (3) Case number vs Age, (4) Seroprevalence
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################

#################################
##  panel "b"
##  Case number vs year
#################################
## data download 
x_Data_for_Figure1b = read.table('./Data/x_Data_for_Figure1b.txt', header = FALSE, sep = "", dec = ".")
y_total_Data_for_Figure1b = read.table('./Data/y_total_Data_for_Figure1b.txt', header = FALSE, sep = "", dec = ".")
y_DE1_Data_for_Figure1b = read.table('./Data/y_DE1_Data_for_Figure1b.txt', header = FALSE, sep = "", dec = ".")
y_DE2_Data_for_Figure1b = read.table('./Data/y_DE2_Data_for_Figure1b.txt', header = FALSE, sep = "", dec = ".")
y_DE3_Data_for_Figure1b = read.table('./Data/y_DE3_Data_for_Figure1b.txt', header = FALSE, sep = "", dec = ".")
y_DE4_Data_for_Figure1b = read.table('./Data/y_DE4_Data_for_Figure1b.txt', header = FALSE, sep = "", dec = ".")

pdf(paste(fig_save_place,"/casenumber_vs_year_equalwidth.pdf",sep=''),width = 8., height = 4)
par(mfrow=c(1,1),
    oma = c(1,2,0,0) + 0.1,
    mar = c(1,2,0,0) + 0.1)
dt_ = x_Data_for_Figure1b[1:nrow(x_Data_for_Figure1b)-1,1]
interval_ = 12*5
lab <- rep(NA, length(dt_))
lab[seq(1, length(dt_), by=interval_)] <-round(dt_[seq(1, length(dt_), by=interval_)])
barplot(height=rbind(y_DE1_Data_for_Figure1b[,1],y_DE2_Data_for_Figure1b[,1],y_DE3_Data_for_Figure1b[,1],y_DE4_Data_for_Figure1b[,1]),
        width=x_Data_for_Figure1b[2:nrow(x_Data_for_Figure1b),1] - x_Data_for_Figure1b[1:nrow(x_Data_for_Figure1b)-1,1],
        space=0,beside=FALSE,col=c("red","blue","green","gold"),
        border = NA,names.arg=lab,
        axisnames = TRUE,las=1,cex.axis =1.6,cex.names = 1.6)
legend( x="topleft", 
        legend=c("DENV-1","DENV-2","DENV-3","DENV-4"), fill=c("red","blue","green","gold"),
        inset=.04,
        col=c("red","blue","green","gold"),cex=1.6)#, lwd=1, lty=c(NA,NA,NA,NA),
dev.off()


#################################
##  panel "c"
##  Case number vs Age
#################################
### Calculate aggregated case numbers for each age for each serotype
case_sum_set_total <-  c(1:nrow(cases_dist))*0
case_sum_set_s1 <- c(1:nrow(cases_dist))*0
case_sum_set_s2 <- c(1:nrow(cases_dist))*0
case_sum_set_s3 <- c(1:nrow(cases_dist))*0
case_sum_set_s4 <- c(1:nrow(cases_dist))*0
for (i in 1:nrow(cases_dist)){
  case_sum_set_total[i] = sum(cases_dist[i,])
  for (j in 1:ncol(cases_dist)){
    if (serotype_number[j]==1){
      case_sum_set_s1[i] = case_sum_set_s1[i] + cases_dist[i,j]
    } else if (serotype_number[j]==2){
      case_sum_set_s2[i] = case_sum_set_s2[i] + cases_dist[i,j]
    } else if (serotype_number[j]==3){
      case_sum_set_s3[i] = case_sum_set_s3[i] + cases_dist[i,j]
    } else if (serotype_number[j]==4){
      case_sum_set_s4[i] = case_sum_set_s4[i] + cases_dist[i,j]
    }
  }
}

type__ = 2 ### 1 : bar graph, 2 : line graph
pdf(paste(fig_save_place,"/casenumber_vs_age.pdf",sep=''),width = 6., height = 5)
par(mfrow=c(1,1),
    oma = c(1,2,0,0) + 0.1,
    mar = c(1,2,0,0) + 0.1)
if (type__ == 1){
  barplot(height=rbind(case_sum_set_s1,case_sum_set_s2,case_sum_set_s3,case_sum_set_s4),
          width=c(1:nrow(cases_dist))*0 + 1,
          space=0,beside=FALSE,col=c("red","blue","green","gold"),
          border = TRUE,names.arg=c(1:nrow(cases_dist)),
          axisnames = TRUE)
  
  legend( x="topright", 
          legend=c("DENV-1","DENV-2","DENV-3","DENV-4"), fill=c("red","blue","green","gold"),
          inset=.04,
          col=c("red","blue","green","gold"),cex=1.7)
} else if (type__ == 2){
  plot(c(1:nrow(cases_dist)),case_sum_set_s1,type='l',col='red',lwd=2,cex.axis =1.7,las=1)
  lines(c(1:nrow(cases_dist)),case_sum_set_s2,type='l',col='blue',lwd=2)
  lines(c(1:nrow(cases_dist)),case_sum_set_s3,type='l',col='green',lwd=2)
  lines(c(1:nrow(cases_dist)),case_sum_set_s4,type='l',col='gold',lwd=2)
  
  legend( x="topright", 
          legend=c("DENV-1","DENV-2","DENV-3","DENV-4"), 
          inset=.04,
          col=c("red","blue","green","gold"),cex=1.7, lwd=2, lty=1) 
}

dev.off()



#################################
##  panel "d"
##  Seroprevalence
#################################
pdf(paste(fig_save_place,"/Seroposi_2014_2015.pdf",sep=''),width = 6., height = 6)
par(mfrow=c(2,1),
    oma = c(1.5,1.5,0,0) + 0.5,
    mar = c(1.5,1.5,0,0) + 0.5)

## seroprevalence study in 2014
n_2014 <- which(Surveilanceyear == 2014)
age_x <- c(1:length(Seroposi_2014))[6:16]
Sero_2014_fit = Seroposi_1_2014 + Seroposi_2_2014 + Seroposi_3_2014 + Seroposi_4_2014 + Seroposi_twice_2014
Binom_infer <- binom.bayes(Sero_2014_fit[6:16],TrialNum_2014[6:16],conf.level=0.95)
plot(age_x,Binom_infer$mean,type="l",xlim=c(6,16),ylim=c(0,1),lty=2,
     lwd=2,xlab="age",ylab="proportion",col="red",
     las=1,cex.axis =1.7)
arrows(age_x, Binom_infer$lower, age_x, Binom_infer$upper, 
       length=0.05, angle=90, code=3,col="red")
legend(x="bottomright",
       c("2014"), 
       lty=c(2,1), 
       inset=.04,
       pch=c(NA,NA),
       lwd=c(2,2),
       col=c("red"),
       cex=1.7,bg='white')

## seroprevalence study in 2015
n_2015 <- which(Surveilanceyear == 2015)
age_x <- c(1:length(Seroposi_2015))
Sero_2015_fit = Seroposi_1_2015 + Seroposi_2_2015 + Seroposi_3_2015 + Seroposi_4_2015 + Seroposi_twice_2015
coarse_num <- 5
Seroposi_1_2015_coarse <- seq(coarse_num,80,coarse_num)
TrialNum_2015_coarse <- seq(coarse_num,80,coarse_num)
age_x_coarse <- seq(coarse_num,80,coarse_num) - coarse_num/2
count_ <- 0
for (i in seq(coarse_num,80,coarse_num)){
  count_ <- count_ + 1
  Seroposi_1_2015_coarse[count_] <- sum(Sero_2015_fit[(i-coarse_num+1):i])
  TrialNum_2015_coarse[count_] <- sum(TrialNum_2015[(i-coarse_num+1):i])
}
Binom_infer <- binom.bayes(Seroposi_1_2015_coarse,TrialNum_2015_coarse,conf.level=0.95)
plot(age_x_coarse,Binom_infer$mean,type="l",xlim=c(0,80),ylim=c(0,1),lty=2,
     lwd=2,xlab="age",ylab="proportion",col="blue",
     las=1,cex.axis =1.7)
arrows(age_x_coarse, Binom_infer$lower, age_x_coarse, Binom_infer$upper, 
       length=0.05, angle=90, code=3,col="blue")
legend(x="bottomright",
       c("2015"), 
       lty=c(2,1), 
       inset=.04,
       pch=c(NA,NA),
       lwd=c(2,2),
       col=c("blue"),
       cex=1.7,bg='white')

dev.off()




#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
## (Fig.2) Reporting rates
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################

############## all parameters (estimate percentiles) #######################################
parameter50 <- numeric(length(xtest)) ## parameter (median)
parameter2p5 <- numeric(length(xtest)) ## parameter (2.5th-percentile)
parameter97p5 <- numeric(length(xtest)) ## parameter (97.5th-percentile)
if (which_data == 'Synthetic'){
  parameter_correct <- xtest_test ## the correct parameters in the case of synthetic data
}
for (i in c(1:length(xtest))){
  parameter50[i] <- quantile(parameters[burnIn:nbIteration,i],probs=0.5) ## estimate the percentiles
  parameter2p5[i] <- quantile(parameters[burnIn:nbIteration,i],probs=0.025)
  parameter97p5[i] <- quantile(parameters[burnIn:nbIteration,i],probs=0.975)
}
pdf(paste(fig_save_place,"/Parameters.pdf",sep=""))
par(mfrow=c(1,1))
plot(c(1:length(xtest)),parameter50,type="points",ylim = c(0,1),xlim=c(0,length(xtest)),xlab = "Parameter index",ylab="Value")
arrows(c(1:length(xtest)), parameter2p5, c(1:length(xtest)), parameter97p5, length=0.05, angle=90, code=3)
if (which_data == 'Synthetic'){
  points(c(1:length(xtest)),parameter_correct, col="red" , pch=4)
}
dev.off()
write.table(parameter50, paste(fig_save_place,"/Parameters_y_50_.txt",sep=""), row.names = FALSE, col.names = FALSE)
write.table(parameter2p5, paste(fig_save_place,"/Parameters_y_2p5_.txt",sep=""), row.names = FALSE, col.names = FALSE)
write.table(parameter97p5, paste(fig_save_place,"/Parameters_y_97p5_.txt",sep=""), row.names = FALSE, col.names = FALSE)
#################################################################



############## PHI (estimate percentiles of reporting probabilities) ##############################################
PHI_background_50 = c(1:length(Surveilanceyear))  # median of T(t)  (Look at the notation in SI.)
PHI_background_2p5 = c(1:length(Surveilanceyear)) 
PHI_background_97p5 = c(1:length(Surveilanceyear))
PHI_Primary_50 = c(1:length(Primary_PHI))  # median of phi(i,1) / phi(1,1)  (Look at the notation in SI.)
PHI_Primary_2p5 = c(1:length(Primary_PHI))
PHI_Primary_97p5 = c(1:length(Primary_PHI))
PHI_Secondary_50 = c(1:length(Secondary_PHI)) # median of phi(i,2) / phi(1,1)
PHI_Secondary_2p5 = c(1:length(Secondary_PHI))
PHI_Secondary_97p5 = c(1:length(Secondary_PHI))
PHI_Secondary_v2_50 = c(1:length(Secondary_PHI)) # median of phi(i,2) / phi(1,2) 
PHI_Secondary_v2_2p5 = c(1:length(Secondary_PHI))
PHI_Secondary_v2_97p5 = c(1:length(Secondary_PHI))
Age_k_50 = c(1:4) # median of ... Age_k_50[1]: k in the negative binomial parameter, Age_k_50[2]: Age factor 1-4, Age_k_50[3]: Age factor 10-14, Age_k_50[4]: Age factor 15+
Age_k_2p5 = c(1:4)
Age_k_97p5 = c(1:4)
if (which_data == 'Synthetic'){ ## the parameters for the true values
  PHI_background_correct = c(1:length(Surveilanceyear))
  PHI_Primary_correct = c(1:length(Primary_PHI))
  PHI_Secondary_correct = c(1:length(Secondary_PHI))
  PHI_crossp_correct = c(1:4)
}


PHI_background_data <- matrix(0,nrow=nbIteration-burnIn,ncol=length(Surveilanceyear))
for (i in c((burnIn+1):nbIteration)){
  steep <- 5.
  ca <- parameters[i,(lambda_length+2)]
  cb <- parameters[i,(lambda_length+3)]
  cc <- parameters[i,(lambda_length+4)]
  ya <- parameters[i,(lambda_length+5)]
  yb <- parameters[i,(lambda_length+6)]
  yc <- parameters[i,(lambda_length+7)]
  ### Be careful ##### I used the value from PHI_Primary to fix the largest value of PHI_background ########
  largest_PHI_background <- parameters[i,lambda_length+length(Background_PHI)+1]
  ### Be careful ##### I used the value from PHI_Primary to fix the largest value of PHI_background ########
  PHI_background_data[i-burnIn,] <- ca*(tanh(steep*(Surveilanceyear - 1985)) + 1.)/2. + cb*(tanh(steep*(Surveilanceyear - 1990)) + 1.)/2. + cc*(tanh(steep*(Surveilanceyear - 1995)) + 1.)/2. - ca - cb - cc
  PHI_background_data[i-burnIn,] <- PHI_background_data[i-burnIn,] + ya*(tanh(steep*(Surveilanceyear - 2000)) + 1.)/2. + yb*(tanh(steep*(Surveilanceyear - 2004)) + 1.)/2. + yc*(tanh(steep*(Surveilanceyear - 2008)) + 1.)/2. + (largest_PHI_background - ya - yb - yc)
}
for (i in c(1:length(Surveilanceyear))){
  PHI_background_50[i] <- quantile(PHI_background_data[,i],probs=0.5)
  PHI_background_2p5[i] <- quantile(PHI_background_data[,i],probs=0.025)
  PHI_background_97p5[i] <- quantile(PHI_background_data[,i],probs=0.975)
}
if (which_data == 'Synthetic'){
  steep <- 5.
  ca <- parameter_correct[(lambda_length+2)]
  cb <- parameter_correct[(lambda_length+3)]
  cc <- parameter_correct[(lambda_length+4)]
  ya <- parameter_correct[(lambda_length+5)]
  yb <- parameter_correct[(lambda_length+6)]
  yc <- parameter_correct[(lambda_length+7)]
  ### Be careful ##### I used the value from PHI_Primary to fix the largest value of PHI_background ########
  largest_PHI_background <- parameter_correct[lambda_length+length(Background_PHI)+1]
  ### Be careful ##### I used the value from PHI_Primary to fix the largest value of PHI_background ########
  PHI_background_correct <- ca*(tanh(steep*(Surveilanceyear - 1985)) + 1.)/2. + cb*(tanh(steep*(Surveilanceyear - 1990)) + 1.)/2. + cc*(tanh(steep*(Surveilanceyear - 1995)) + 1.)/2. - ca - cb - cc
  PHI_background_correct <- PHI_background_correct + ya*(tanh(steep*(Surveilanceyear - 2000)) + 1.)/2. + yb*(tanh(steep*(Surveilanceyear - 2004)) + 1.)/2. + yc*(tanh(steep*(Surveilanceyear - 2008)) + 1.)/2. + (largest_PHI_background - ya - yb - yc)
}

for (i in c(1:length(Primary_PHI))){
  PHI_Primary_50[i] = quantile(parameters[burnIn:nbIteration,i+lambda_length+length(Background_PHI)],probs=0.5)
  PHI_Primary_2p5[i] = quantile(parameters[burnIn:nbIteration,i+lambda_length+length(Background_PHI)],probs=0.025)
  PHI_Primary_97p5[i] = quantile(parameters[burnIn:nbIteration,i+lambda_length+length(Background_PHI)],probs=0.975)
  PHI_Primary_50[1] = 1.
  PHI_Primary_2p5[1] = 1.
  PHI_Primary_97p5[1] = 1.
  if (which_data == 'Synthetic'){
    PHI_Primary_correct[i]<- parameter_correct[i+lambda_length+length(Background_PHI)]
  }
}
for (i in c(1:length(Secondary_PHI))){
  PHI_Secondary_50[i] = quantile(parameters[burnIn:nbIteration,i+lambda_length+length(Background_PHI)+length(Primary_PHI)],probs=0.5)
  PHI_Secondary_2p5[i] = quantile(parameters[burnIn:nbIteration,i+lambda_length+length(Background_PHI)+length(Primary_PHI)],probs=0.025)
  PHI_Secondary_97p5[i] = quantile(parameters[burnIn:nbIteration,i+lambda_length+length(Background_PHI)+length(Primary_PHI)],probs=0.975)
  dummy_ = parameters[burnIn:nbIteration,i+lambda_length+length(Background_PHI)+length(Primary_PHI)]/parameters[burnIn:nbIteration,1+lambda_length+length(Background_PHI)+length(Primary_PHI)]
  PHI_Secondary_v2_50[i] = quantile(dummy_,probs=0.5)
  PHI_Secondary_v2_2p5[i] = quantile(dummy_,probs=0.025)
  PHI_Secondary_v2_97p5[i] = quantile(dummy_,probs=0.975)
  if (which_data == 'Synthetic'){
    PHI_Secondary_correct[i]<- parameter_correct[i+lambda_length+length(Background_PHI)+length(Primary_PHI)]
  }
}

crossp_data <- matrix(0,nrow=nbIteration-burnIn,ncol=4)
for (i in c((burnIn+1):nbIteration)){
  crossp1 <- parameters[i,(lambda_length+8)]
  crossp2 <- parameters[i,(lambda_length+9)]
  crossp3 <- parameters[i,(lambda_length+10)]
  crossp4 <- parameters[i,(lambda_length+11)]
  crossp_data[i-burnIn,] <- c(crossp1,crossp2,crossp3,crossp4)
}

for (i in c(1:4)){
  Age_k_50[i] <- quantile(crossp_data[,i],probs=0.5)
  Age_k_2p5[i] <- quantile(crossp_data[,i],probs=0.025)
  Age_k_97p5[i] <- quantile(crossp_data[,i],probs=0.975)
}
if (which_data == 'Synthetic'){
  crossp1 <- parameter_correct[(lambda_length+8)]
  crossp2 <- parameter_correct[(lambda_length+9)]
  crossp3 <- parameter_correct[(lambda_length+10)]
  crossp4 <- parameter_correct[(lambda_length+11)]
  PHI_crossp_correct <- c(crossp1,crossp2,crossp3,crossp4)
}  





### for panel a
pdf(paste(fig_save_place,"/PHI_background.pdf",sep=""),width = 5.5, height = 4)
par(mfrow=c(1,1))
dum_ <- which(Surveilanceyear<2015)
plot(Surveilanceyear[dum_],PHI_background_50[dum_]*PHI_Primary_50[1],type="l",
     ylim = c(0.005,0.072),xlim=c(1979,2014),xlab = "",ylab="",las=1,cex.axis=1.2)#,ylab="Primary reporting rate of DENV1",main=Model_description[3])
polygon(c(Surveilanceyear[dum_],rev(Surveilanceyear[dum_])),c(PHI_background_2p5[dum_]*PHI_Primary_50[1],rev(PHI_background_97p5[dum_]*PHI_Primary_50[1])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),border=NA)
dev.off()

primary_x <- c(1:length(PHI_Primary_50))
### for panel c
pdf(paste(fig_save_place,"/PHI_Primary_Secondary_par_eachDE1.pdf",sep=""),width = 6.*0.85, height = 5*0.85)
plot(primary_x-0.1,PHI_Primary_50/PHI_Primary_50[1],
     xlab = "",ylab="",type="p",xlim = c(0.9,4.1),ylim = c(0.01,6.),
     xaxt = "n",col="black",las=1,
     cex.axis =1.25)
points(c(1,2,3,4)+0.1,PHI_Secondary_v2_50[c(1,4,7,10)]/PHI_Secondary_v2_50[1],
     xlab = "",ylab="",type="p",ylim = c(0.01,6.),xaxt = "n",col="red")
arrows(primary_x-0.1, PHI_Primary_2p5, primary_x-0.1, PHI_Primary_97p5, length=0.05, angle=90, code=3,col="black")
arrows(c(1,2,3,4)+0.1, PHI_Secondary_v2_2p5[c(1,4,7,10)], c(1,2,3,4)+0.1, PHI_Secondary_v2_97p5[c(1,4,7,10)], length=0.05, angle=90, code=3,col="red")
axis(1, at=1:4, labels=c("DENV-1","DENV-2","DENV-3","DENV-4"),cex.axis = 1.25)  
abline(h = 1, lty = 3)
legend(x="topleft",
       c("Primary",
         "Secondary"), 
       pch=c(1,1),
       col=c("black","red"),
       cex=1.25)
dev.off()


## for panel b
pdf(paste(fig_save_place,"/PHI_Secondary_par_priDE1.pdf",sep=""),width = 3.2, height = 3.7)
plot(c(1,2),c(1,PHI_Secondary_50[c(1)]),type="p",
     xlab = "",ylab="",xaxt = "n",xlim=c(0.5,2.5),ylim = c(0.01,5.3),las=1)
arrows(c(1,2), c(1,PHI_Secondary_2p5[c(1)]), c(1,2), c(1,PHI_Secondary_97p5[c(1)]), length=0.05, angle=90, code=3)
abline(h = 1, lty = 2)
axis(1, at=c(1,2), labels=c("Primary","Secondary"))  
dev.off()


## for panel d
pdf(paste(fig_save_place,"/Age_dependent_par_5_10.pdf",sep=""),width = 5., height = 5)
age_fac_50 <- c(Age_k_50[2],1,Age_k_50[3],Age_k_50[4])
age_fac_2p5 <- c(Age_k_2p5[2],1,Age_k_2p5[3],Age_k_2p5[4])
age_fac_97p5 <- c(Age_k_97p5[2],1,Age_k_97p5[3],Age_k_97p5[4])
plot(c(1:4),age_fac_50,xlab = "", #"Age",
     type="l",ylab="",ylim = c(0.3,1.2),xaxt = "n",main="",las=1,cex.axis=1.5)#,log="y")
polygon(c(c(1:4),rev(c(1:4))),c(age_fac_2p5,rev(age_fac_97p5)),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),border=NA)
axis(1, at=1:4, labels=c("1-4","5-9","10-14","15-80"),cex.axis=1.5)
dev.off()


## k for negative binomial parameter
####################
pdf(paste(fig_save_place,"/k_value_negbin.pdf",sep=""),width = 2., height = 4)
plot(c(1),Age_k_50[1],main="",ylab="",xlab="",ylim=c(0.55,0.8),xaxt='n',las=1)
arrows(c(1), Age_k_2p5[1], c(1), Age_k_97p5[1], 
       length=0.05, angle=90, code=3,col="black")
dev.off()


############# Calculate the mean value of the reporting rates
Pri_rep = PHI_Primary_50
Sec_rep = PHI_Secondary_v2_50[c(1,4,7,10)] * PHI_Secondary_50[1]
print(mean(Sec_rep) / mean(Pri_rep))





#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
## Fig.3 
## (i) Case number, (ii) FOI and (iii) immunity. 
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################


############### Simulated case data: using the parameters obtained during MCMC, simulate the case data ###################### 
Resul_data <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = length(as.vector(cases_dist)))
count_ <- 0
for (i in seq(burnIn,nbIteration,50)) {
  count_ <- count_ + 1
  print(i)
  x_infered <- parameters[i,]
  Resul_data[count_,] <- Simulation_casedata(x_infered) ## simulate the case counts from estimated parameters
}
Resul_data_50 <- c(1:length(Resul_data[1,]))
Resul_data_2p5 <- c(1:length(Resul_data[1,]))
Resul_data_97p5 <- c(1:length(Resul_data[1,]))
for (i in c(1:length(Resul_data[1,]))){
  Resul_data_50[i] <- quantile(Resul_data[,i],probs = 0.5)
  Resul_data_2p5[i] <- quantile(Resul_data[,i],probs = 0.025)
  Resul_data_97p5[i] <- quantile(Resul_data[,i],probs = 0.975)
}

Full_roof_age_middle <- 1:(length(roof_age)-1)
for (i in  1:(length(roof_age)-1)){
  Full_roof_age_middle[i] <- c(floor((roof_age[i]+1+roof_age[i+1])/2))
}
Dummmy_case_dist <- cases_dist
Simulated_count_50 = matrix(Resul_data_50,ncol=ncol(cases_dist),nrow=nrow(cases_dist))
Simulated_count_2p5 = matrix(Resul_data_2p5,ncol=ncol(cases_dist),nrow=nrow(cases_dist))
Simulated_count_97p5 = matrix(Resul_data_97p5,ncol=ncol(cases_dist),nrow=nrow(cases_dist))

Simulated_count_sum_50 = Surveilanceyear[1:which(Surveilanceyear==2014)]
Simulated_count_sum_2p5 = Surveilanceyear[1:which(Surveilanceyear==2014)]
Simulated_count_sum_97p5 = Surveilanceyear[1:which(Surveilanceyear==2014)]
for (i in 1:which(Surveilanceyear==2014)){
  Simulated_count_sum_50[i] <-  sum(Simulated_count_50[,i])
  Simulated_count_sum_2p5[i] <-  sum(Simulated_count_2p5[,i])
  Simulated_count_sum_97p5[i] <-  sum(Simulated_count_97p5[,i])
}




######## Extract the information of x0, mono1, mono2, mono3, mono4 (where x0 is x and monoi is y_i (i=1,2,3,4) in the notation of the article) for all the obtained parameters during MCMC and then estimate the median of these values.
Resul_data_x0 <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = L*TT)
Resul_data_m1 <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = L*TT)
Resul_data_m2 <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = L*TT)
Resul_data_m3 <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = L*TT)
Resul_data_m4 <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = L*TT)
Resul_data_twice <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = L*TT)
count_ <- 0
for (i in seq(burnIn,nbIteration,50)) {
  count_ <- count_ + 1
  print(i)
  x_infered <- parameters[i,]
  dummy_x0_mono1_mono2_mono3_mono4 <- x0_mono1_mono2_mono3_mono4(x_infered) ## from the parameters, obtain x0, mono1, mono2, ..., mono4.
  Resul_data_x0[count_,] <- dummy_x0_mono1_mono2_mono3_mono4[1:(1*L*TT)]
  Resul_data_m1[count_,] <- dummy_x0_mono1_mono2_mono3_mono4[(1+1*L*TT):(2*L*TT)]
  Resul_data_m2[count_,] <- dummy_x0_mono1_mono2_mono3_mono4[(1+2*L*TT):(3*L*TT)]
  Resul_data_m3[count_,] <- dummy_x0_mono1_mono2_mono3_mono4[(1+3*L*TT):(4*L*TT)]
  Resul_data_m4[count_,] <- dummy_x0_mono1_mono2_mono3_mono4[(1+4*L*TT):(5*L*TT)]
  morethantwice__ <- 1-dummy_x0_mono1_mono2_mono3_mono4[1:(1*L*TT)] - dummy_x0_mono1_mono2_mono3_mono4[(1+1*L*TT):(2*L*TT)] - dummy_x0_mono1_mono2_mono3_mono4[(1+2*L*TT):(3*L*TT)] - dummy_x0_mono1_mono2_mono3_mono4[(1+3*L*TT):(4*L*TT)] - dummy_x0_mono1_mono2_mono3_mono4[(1+4*L*TT):(5*L*TT)]
  Resul_data_twice[count_,] <- morethantwice__
}
x0_data_50 <- c(1:(L*TT))
mono1_data_50 <- c(1:(L*TT))
mono2_data_50 <- c(1:(L*TT))
mono3_data_50 <- c(1:(L*TT))
mono4_data_50 <- c(1:(L*TT))
twice_data_50 <- c(1:(L*TT))
for (i in c(1:(L*TT))){
  x0_data_50[i] <- quantile(Resul_data_x0[,i],probs = 0.5)
  mono1_data_50[i] <- quantile(Resul_data_m1[,i],probs = 0.5)
  mono2_data_50[i] <- quantile(Resul_data_m2[,i],probs = 0.5)
  mono3_data_50[i] <- quantile(Resul_data_m3[,i],probs = 0.5)
  mono4_data_50[i] <- quantile(Resul_data_m4[,i],probs = 0.5)
  twice_data_50[i] <- quantile(Resul_data_twice[,i],probs = 0.5)
}

x0_50 = matrix(x0_data_50,nrow = L, ncol = TT)
mono1_50 = matrix(mono1_data_50,nrow = L, ncol = TT)
mono2_50 = matrix(mono2_data_50,nrow = L, ncol = TT)
mono3_50 = matrix(mono3_data_50,nrow = L, ncol = TT)
mono4_50 = matrix(mono4_data_50,nrow = L, ncol = TT)
twice_50 = matrix(twice_data_50,nrow = L, ncol = TT)

### Extend the population data to 2015- to calculate the weighted average. 
pop_dist_detail <- pop_dist
for (i in 1:length(Surveilanceyear)){
  if (Surveilanceyear[i] > 2014){
    iidum <- which(Surveilanceyear == 2014)
    pop_dist_detail <- cbind(pop_dist_detail,pop_dist[,iidum])
  }
}



############ FOI ################################################ 
## Using FOI from the obtained inference results (parameters), calculate percentiles.
Year_plot <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_1_50 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_1_2p5 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_1_97p5 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_2_50 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_2_2p5 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_2_97p5 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_3_50 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_3_2p5 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_3_97p5 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_4_50 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_4_2p5 <- Year_[which(Year_==1944):which(Year_==2018)]
FOI_4_97p5 <- Year_[which(Year_==1944):which(Year_==2018)]
if (which_data == 'Synthetic'){
  FOI_1_correct <- Year_[which(Year_==1944):which(Year_==2018)]
  FOI_2_correct <- Year_[which(Year_==1944):which(Year_==2018)]
  FOI_3_correct <- Year_[which(Year_==1944):which(Year_==2018)]
  FOI_4_correct <- Year_[which(Year_==1944):which(Year_==2018)]
}

counter_ <- 0
for (i in c(which(Year_plot == 1979):which(Year_plot == 2018))){
  counter_ <- counter_ + 1
  if (serotype_number[counter_] == 1){
    FOI_1_50[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.5)
    FOI_1_2p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.025)
    FOI_1_97p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.975)
    FOI_2_50[i] <- 0.
    FOI_2_2p5[i] <- 0.
    FOI_2_97p5[i] <- 0.
    FOI_3_50[i] <- 0.
    FOI_3_2p5[i] <- 0.
    FOI_3_97p5[i] <- 0.
    FOI_4_50[i] <- 0.
    FOI_4_2p5[i] <- 0.
    FOI_4_97p5[i] <- 0.
    if (which_data == 'Synthetic'){
      FOI_1_correct[i] <- parameter_correct[counter_] 
      FOI_2_correct[i] <- 0.
      FOI_3_correct[i] <- 0.
      FOI_4_correct[i] <- 0.
    }
  } else  if (serotype_number[counter_] == 2){
    FOI_1_50[i] <- 0.
    FOI_1_2p5[i] <- 0.
    FOI_1_97p5[i] <- 0.
    FOI_2_50[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.5)
    FOI_2_2p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.025)
    FOI_2_97p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.975)
    FOI_3_50[i] <- 0.
    FOI_3_2p5[i] <- 0.
    FOI_3_97p5[i] <- 0.
    FOI_4_50[i] <- 0.
    FOI_4_2p5[i] <- 0.
    FOI_4_97p5[i] <- 0.
    if (which_data == 'Synthetic'){
      FOI_1_correct[i] <- 0.
      FOI_2_correct[i] <- parameter_correct[counter_] 
      FOI_3_correct[i] <- 0.
      FOI_4_correct[i] <- 0.
    }
  } else  if (serotype_number[counter_] == 3){
    FOI_1_50[i] <- 0.
    FOI_1_2p5[i] <- 0.
    FOI_1_97p5[i] <- 0.
    FOI_2_50[i] <- 0.
    FOI_2_2p5[i] <- 0.
    FOI_2_97p5[i] <- 0.
    FOI_3_50[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.5)
    FOI_3_2p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.025)
    FOI_3_97p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.975)
    FOI_4_50[i] <- 0.
    FOI_4_2p5[i] <- 0.
    FOI_4_97p5[i] <- 0.
    if (which_data == 'Synthetic'){
      FOI_1_correct[i] <- 0.
      FOI_2_correct[i] <- 0.
      FOI_3_correct[i] <- parameter_correct[counter_] 
      FOI_4_correct[i] <- 0.
    }
  } else if (serotype_number[counter_] == 4){
    FOI_1_50[i] <- 0.
    FOI_1_2p5[i] <- 0.
    FOI_1_97p5[i] <- 0.
    FOI_2_50[i] <- 0.
    FOI_2_2p5[i] <- 0.
    FOI_2_97p5[i] <- 0.
    FOI_3_50[i] <- 0.
    FOI_3_2p5[i] <- 0.
    FOI_3_97p5[i] <- 0.
    FOI_4_50[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.5)
    FOI_4_2p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.025)
    FOI_4_97p5[i] <- quantile(parameters[burnIn:nbIteration,counter_],probs=0.975)
    if (which_data == 'Synthetic'){
      FOI_1_correct[i] <- 0.
      FOI_2_correct[i] <- 0.
      FOI_3_correct[i] <- 0.
      FOI_4_correct[i] <- parameter_correct[counter_] 
    }
  } else if (serotype_number[counter_] == 0){ ## this condition is for 2015-2018. (Only 2015 is relevant.) 
    FOI_1_50[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-3],probs=0.5)
    FOI_1_2p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-3],probs=0.025)
    FOI_1_97p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-3],probs=0.975)
    FOI_2_50[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-2],probs=0.5)
    FOI_2_2p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-2],probs=0.025)
    FOI_2_97p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-2],probs=0.975)
    FOI_3_50[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-1],probs=0.5)
    FOI_3_2p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-1],probs=0.025)
    FOI_3_97p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)-1],probs=0.975)
    FOI_4_50[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)],probs=0.5)
    FOI_4_2p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)],probs=0.025)
    FOI_4_97p5[i] <- quantile(parameters[burnIn:nbIteration,length(Surveilanceyear)],probs=0.975)
    if (which_data == 'Synthetic'){
      FOI_1_correct[i] <- parameter_correct[length(Surveilanceyear)-3]
      FOI_2_correct[i] <- parameter_correct[length(Surveilanceyear)-2]
      FOI_3_correct[i] <- parameter_correct[length(Surveilanceyear)-1]
      FOI_4_correct[i] <- parameter_correct[length(Surveilanceyear)] 
    }
  }
}
## Before 1978, FOI = 0 
FOI_1_50[1:which(Year_plot == 1978)] <- 0.
FOI_1_2p5[1:which(Year_plot == 1978)] <- 0.
FOI_1_97p5[1:which(Year_plot == 1978)] <- 0.
FOI_2_50[1:which(Year_plot == 1978)] <- 0.
FOI_2_2p5[1:which(Year_plot == 1978)] <- 0.
FOI_2_97p5[1:which(Year_plot == 1978)] <- 0.
FOI_3_50[1:which(Year_plot == 1978)] <- 0.
FOI_3_2p5[1:which(Year_plot == 1978)] <- 0.
FOI_3_97p5[1:which(Year_plot == 1978)] <- 0.
FOI_4_50[1:which(Year_plot == 1978)] <- 0.
FOI_4_2p5[1:which(Year_plot == 1978)] <- 0.
FOI_4_97p5[1:which(Year_plot == 1978)] <- 0.

## except when there is an epidemic.
FOI_1_50[which(Year_plot == 1944)] <- quantile(parameters[burnIn:nbIteration,lambda_length],probs=0.5)
FOI_1_2p5[which(Year_plot == 1944)] <- quantile(parameters[burnIn:nbIteration,lambda_length],probs=0.025)
FOI_1_97p5[which(Year_plot == 1944)] <- quantile(parameters[burnIn:nbIteration,lambda_length],probs=0.975)
FOI_3_50[which(Year_plot == 1964):which(Year_plot == 1969)] <- quantile(parameters[burnIn:nbIteration,lambda_length-1],probs=0.5)
FOI_3_2p5[which(Year_plot == 1964):which(Year_plot == 1969)] <- quantile(parameters[burnIn:nbIteration,lambda_length-1],probs=0.025)
FOI_3_97p5[which(Year_plot == 1964):which(Year_plot == 1969)] <- quantile(parameters[burnIn:nbIteration,lambda_length-1],probs=0.975)
FOI_2_50[which(Year_plot == 1971):which(Year_plot == 1975)] <- quantile(parameters[burnIn:nbIteration,lambda_length-2],probs=0.5)
FOI_2_2p5[which(Year_plot == 1971):which(Year_plot == 1975)] <- quantile(parameters[burnIn:nbIteration,lambda_length-2],probs=0.025)
FOI_2_97p5[which(Year_plot == 1971):which(Year_plot == 1975)] <- quantile(parameters[burnIn:nbIteration,lambda_length-2],probs=0.975)
FOI_1_50[which(Year_plot == 1976):which(Year_plot == 1978)] <- quantile(parameters[burnIn:nbIteration,lambda_length-3],probs=0.5)
FOI_1_2p5[which(Year_plot == 1976):which(Year_plot == 1978)] <- quantile(parameters[burnIn:nbIteration,lambda_length-3],probs=0.025)
FOI_1_97p5[which(Year_plot == 1976):which(Year_plot == 1978)] <- quantile(parameters[burnIn:nbIteration,lambda_length-3],probs=0.975)

if (which_data == 'Synthetic'){
  FOI_1_correct[1:which(Year_plot == 1978)] <- 0.
  FOI_2_correct[1:which(Year_plot == 1978)] <- 0.
  FOI_3_correct[1:which(Year_plot == 1978)] <- 0.
  FOI_4_correct[1:which(Year_plot == 1978)] <- 0.
  
  FOI_1_correct[which(Year_plot == 1944)] <- parameter_correct[lambda_length]
  FOI_3_correct[which(Year_plot == 1964):which(Year_plot == 1969)] <- parameter_correct[lambda_length-1]
  FOI_2_correct[which(Year_plot == 1971):which(Year_plot == 1975)] <- parameter_correct[lambda_length-2]
  FOI_1_correct[which(Year_plot == 1976):which(Year_plot == 1978)] <- parameter_correct[lambda_length-3] 
}
Year_plot_1 <- Year_plot
Year_plot <- Year_plot_1[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_1_50 <- FOI_1_50[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_1_2p5 <- FOI_1_2p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_1_97p5 <- FOI_1_97p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_2_50 <- FOI_2_50[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_2_2p5 <- FOI_2_2p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_2_97p5 <- FOI_2_97p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_3_50 <- FOI_3_50[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_3_2p5 <- FOI_3_2p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_3_97p5 <- FOI_3_97p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_4_50 <- FOI_4_50[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_4_2p5 <- FOI_4_2p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
FOI_4_97p5 <- FOI_4_97p5[which(Year_plot_1==1979):which(Year_plot_1==2014)]
if (which_data == 'Synthetic'){
  FOI_1_correct <- FOI_1_correct[which(Year_plot_1==1979):which(Year_plot_1==2014)]
  FOI_2_correct <- FOI_2_correct[which(Year_plot_1==1979):which(Year_plot_1==2014)]
  FOI_3_correct <- FOI_3_correct[which(Year_plot_1==1979):which(Year_plot_1==2014)]
  FOI_4_correct <- FOI_4_correct[which(Year_plot_1==1979):which(Year_plot_1==2014)]
}

### (From the data) calculate the number of total reported cases for each serotype between 1979 and 2014.
case_sum_set_s1 <- numeric(length(1:which(Surveilanceyear==2014)))
case_sum_set_s2 <- numeric(length(1:which(Surveilanceyear==2014)))
case_sum_set_s3 <- numeric(length(1:which(Surveilanceyear==2014)))
case_sum_set_s4 <- numeric(length(1:which(Surveilanceyear==2014)))
for (i in 1:which(Surveilanceyear==2014)){
  if (serotype_number[i]==1){
    case_sum_set_s1[i] <- sum(cases_dist[,i])
  } else if (serotype_number[i]==2){
    case_sum_set_s2[i] <- sum(cases_dist[,i])
  } else if (serotype_number[i]==3){
    case_sum_set_s3[i] <- sum(cases_dist[,i])
  } else if (serotype_number[i]==4){
    case_sum_set_s4[i] <- sum(cases_dist[,i])
  }
}
###




## make figure 
pdf(paste(fig_save_place,"/FitFig1_FOI_Immunity.pdf",sep=""),width = 8.*0.8, height = 6.4*0.8)
par(mfrow=c(3,1),
    oma = c(4,3,0,0) + 0.1,
    mar = c(2,2,2,2) + 0.1)
#################
### Plot Case number
#################
colset <- c("red","blue","green","gold")
diff_surveilance = Surveilanceyear[2:which(Surveilanceyear==2015)]-Surveilanceyear[1:which(Surveilanceyear==2014)]
plot(Surveilanceyear[1:which(Surveilanceyear==2014)],Simulated_count_sum_50/diff_surveilance,type="l",xlab="year",
     xlim=c(1979,2014),ylim=c(0,5000),cex.axis=1.2,las=1)#,main="Case number")
polygon(c(Surveilanceyear[1:which(Surveilanceyear==2014)],rev(Surveilanceyear[1:which(Surveilanceyear==2014)])),c(Simulated_count_sum_2p5/diff_surveilance,rev(Simulated_count_sum_97p5/diff_surveilance)),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),border=NA)
case_sum_set <- Surveilanceyear[which(Surveilanceyear<=2014)]
for (i in 1:which(Surveilanceyear==2014)){
  case_sum_set[i] <- sum(cases_dist[,i])
  points(Surveilanceyear[i],case_sum_set[i]/(Surveilanceyear[i+1]-Surveilanceyear[i]),xlab="year",pch=19,cex=1.,
         ylab="Reported Number / year",col=colset[serotype_number[i]])#,lwd=2,main=paste(round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""))#type="lines",
}

legend( x="topleft", 
        legend=c("DENV-1 (Data)","DENV-2 (Data)",
                 "DENV-3 (Data)","DENV-4 (Data)",
                 "Simulations"), 
        col=c("red","blue","green","gold","black"), lwd=1, lty=c(NA,NA,NA,NA,1), 
        pch=c(19,19,19,19,NA), merge=FALSE )


##############################
######## Plot Force of infection
##############################
Ymaxvalue <- 1
plot(Year_plot,FOI_1_50,type="l",
     xlim=c(1979,2014),ylim=c(0,Ymaxvalue),col="red",
     xlab = "",ylab="",cex.axis=1.2,las=1)
polygon(c(Year_plot,rev(Year_plot)),c(FOI_1_2p5,rev(FOI_1_97p5)),col=rgb(red=1.,green=0.,blue=0.,alpha=0.3),border=NA)
#
lines(Year_plot,FOI_2_50,type="l",xlim=c(1945,2018),col="blue",
      ylim=c(0,Ymaxvalue),xlab = "",
      ylab="")
polygon(c(Year_plot,rev(Year_plot)),c(FOI_2_2p5,rev(FOI_2_97p5)),col=rgb(red=0.,green=0.,blue=1.,alpha=0.3),border=NA)
#####
lines(Year_plot,FOI_3_50,type="l",col="green",
      xlim=c(1945,2018),ylim=c(0,Ymaxvalue),xlab = "",
      ylab="")
polygon(c(Year_plot,rev(Year_plot)),c(FOI_3_2p5,rev(FOI_3_97p5)),col=rgb(red=0.,green=1.,blue=0.,alpha=0.3),border=NA)
######
lines(Year_plot,FOI_4_50,type="l",xlim=c(1945,2018),ylim=c(0,Ymaxvalue),col="gold",
      xlab = "",ylab="")
polygon(c(Year_plot,rev(Year_plot)),c(FOI_4_2p5,rev(FOI_4_97p5)),col=rgb(red=1.,green=1.,blue=0.,alpha=0.5),border=NA)
legend( x="topleft", 
        legend=c("DENV-1","DENV-2",
                 "DENV-3","DENV-4"), 
        col=c("red","blue","green","gold"), lwd=1, lty=c(1,1,1,1,1), 
        pch=c(NA,NA,NA,NA,NA), merge=FALSE )


################
##### Plot Immunity
################
count_ <- 0
for (ave_range in list(c(1:80))){
  count_ <- count_ + 1
  ### All average
  x0_50_mean <- c()
  mono1_50_mean <- c()
  mono2_50_mean <- c()
  mono3_50_mean <- c()
  mono4_50_mean <- c()
  for (i in 1:ncol(x0_50)){
    x0_50_mean[i] <- weighted.mean(x0_50[ave_range,i],pop_dist_detail[ave_range,i])
    mono1_50_mean[i] <- weighted.mean(mono1_50[ave_range,i],pop_dist_detail[ave_range,i])
    mono2_50_mean[i] <- weighted.mean(mono2_50[ave_range,i],pop_dist_detail[ave_range,i])
    mono3_50_mean[i] <- weighted.mean(mono3_50[ave_range,i],pop_dist_detail[ave_range,i])
    mono4_50_mean[i] <- weighted.mean(mono4_50[ave_range,i],pop_dist_detail[ave_range,i])
  }
  rl = which(Surveilanceyear==2014)
  plot(-10, -10,type="p",col="black",xlim=c(1979,2014),ylim=c(0,1),lwd=3,xlab="Year",ylab="Proportion",cex.axis=1.2,las=1)
  l_a <- x0_50_mean[1:rl]
  print(l_a)
  polygon(c(Surveilanceyear[1:rl],rev(Surveilanceyear[1:rl])),c(rep(0,length(l_a)),rev(l_a)),col=rgb(red=0.8,green=0.8,blue=0.8,alpha=1.),border=NA)
  l_b <- x0_50_mean[1:rl] + mono1_50_mean[1:rl]
  polygon(c(Surveilanceyear[1:rl],rev(Surveilanceyear[1:rl])),c(l_a,rev(l_b)),col=rgb(red=1.,green=0.,blue=0.,alpha=1.),border=NA)
  l_c <- l_b + mono2_50_mean[1:rl]
  polygon(c(Surveilanceyear[1:rl],rev(Surveilanceyear[1:rl])),c(l_b,rev(l_c)),col=rgb(red=0.,green=0.,blue=1.,alpha=1.),border=NA)
  l_d <- l_c + mono3_50_mean[1:rl]
  polygon(c(Surveilanceyear[1:rl],rev(Surveilanceyear[1:rl])),c(l_c,rev(l_d)),col=rgb(red=0.,green=1.,blue=0.,alpha=1.),border=NA)
  l_e <- l_d + mono4_50_mean[1:rl]
  polygon(c(Surveilanceyear[1:rl],rev(Surveilanceyear[1:rl])),c(l_d,rev(l_e)),col=rgb(red=1.,green=1.,blue=0.,alpha=1.),border=NA)
  print(1-l_e)
  polygon(c(Surveilanceyear[1:rl],rev(Surveilanceyear[1:rl])),c(l_e,rep(1,length(l_a))),col=rgb(red=0.,green=0.,blue=0.,alpha=1.),border=NA)
  if (count_ == 1){
    legend(x="bottomleft",
           c("never","once (DENV-1)",
             "once (DENV-2)",
             "once (DENV-3)",
             "once (DENV-4)",
             "twice"), 
           fill=c("gray", "red", "blue","green","gold","black"),  cex=1, bty = "n",
           bg="transparent",ncol=3)
  }
}
dev.off()






pdf(paste(fig_save_place,"/Histogram_reportedNumbers.pdf",sep=""),width = 8.*0.8, height = 6.4*0.5)
hist(case_sum_set,breaks=100,xlab="Number of reported cases for each period",ylab="Frequency",cex.axis=1.2,las=1)
dev.off()

#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
## Fig.4
## (i) Age distribution, (ii) Seroprevalence
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################


########################
### Plot age distribution 
########################
for (i in 1:which(Surveilanceyear==2014)){
  pdf(paste(fig_save_place,"/Age_dist_",Surveilanceyear[i],"_.pdf",sep=""),width = 4.5*0.85, height = 3.*0.85)
  par(mfrow=c(1,1))
  plot(Full_roof_age_middle,cases_dist[,i],xlab="",
       pch=19,cex=1.,ylab="",col=colset[serotype_number[i]],lwd=2,main="",
       ylim = c(0,max(max(cases_dist[,i]),max(Simulated_count_97p5[,i]))),
       xlim = c(2,60),las=1)
  polygon(c(Full_roof_age_middle,rev(Full_roof_age_middle)),c(Simulated_count_2p5[,i],rev(Simulated_count_97p5[,i])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),border=NA)
  points(Full_roof_age_middle,cases_dist[,i],xlab="age",
         pch=19,cex=1.,ylab="Reported Number",col=colset[serotype_number[i]],lwd=2,main=paste(round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),#type="lines",
         ylim = c(0,max(max(cases_dist[,i]),max(Simulated_count_50[,i]))))
  lines(Full_roof_age_middle,Simulated_count_50[,i],lwd=2,col="black")
  if (Surveilanceyear[i] == 1979){
    legend("topright",paste("", round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  } else if (Surveilanceyear[i] == 1988.9){
    legend("topright",paste("",round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  } else if (Surveilanceyear[i] == 1989.5){
    legend("topright",paste("",round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  } else if (Surveilanceyear[i] == 1996.67){
    legend("topright",paste("",round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  } else if (Surveilanceyear[i] == 2001.08){
    legend("topright",paste("",round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  } else if (Surveilanceyear[i] == 2006.75){
    legend("topright",paste("",round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  } else if (Surveilanceyear[i] == 2009.16){
    legend("topright",paste("",round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  } else if (Surveilanceyear[i] == 2014){
    legend("topright",paste("",round(Surveilanceyear[i], 1),"-",round(Surveilanceyear[i+1], 1)," DENV",serotype_number[i],sep=""),bty = "n")
  }
  dev.off()
}





# To estimate x0_shifted (Since the seroprevalence survey is not conducted on 1st of January 2014 nor 1st of January 2015, we calculate x0(2014.42) and x0(2015.75)) to compare the inferred result with the data.
Resul_data_x0_sero <- matrix(0, nrow = length(parameters[seq(burnIn,nbIteration,50),1]), ncol = L*TT)
count_ <- 0
for (i in seq(burnIn,nbIteration,50)) {
  count_ <- count_ + 1
  print(i)
  x_infered <- parameters[i,]
  dummy_x0_mono1_mono2_mono3_mono4 <- Serovositivity_casedata(x_infered)
  Resul_data_x0_sero[count_,] <- dummy_x0_mono1_mono2_mono3_mono4[1:(1*L*TT)]
}
x0_sero_data_50 <- c(1:(L*TT))
x0_sero_data_2p5 <- c(1:(L*TT))
x0_sero_data_97p5 <- c(1:(L*TT))
for (i in c(1:(L*TT))){
  x0_sero_data_50[i] <- quantile(Resul_data_x0_sero[,i],probs = 0.5)
  x0_sero_data_2p5[i] <- quantile(Resul_data_x0_sero[,i],probs = 0.025)
  x0_sero_data_97p5[i] <- quantile(Resul_data_x0_sero[,i],probs = 0.975)
}

x0_sero_2p5 = matrix(x0_sero_data_2p5,nrow = L, ncol = TT)
x0_sero_50 = matrix(x0_sero_data_50,nrow = L, ncol = TT)
x0_sero_97p5 = matrix(x0_sero_data_97p5,nrow = L, ncol = TT)


########################
### Plot seroprevalence 
########################
# 2014 seroprevalence
age_max <- length(mono1_50[,1])
age_set <- c(1:age_max)
n_2014 <- which(Surveilanceyear == 2014)
age_x <- c(1:length(Seroposi_2014))
pdf(paste(fig_save_place,"/Seroposi_compare_2014_All1or0.pdf",sep=""),width = 7.*0.7, height = 5*0.7)
Sero_2014_fit = Seroposi_1_2014 + Seroposi_2_2014 + Seroposi_3_2014 + Seroposi_4_2014 + Seroposi_twice_2014
plot(age_set,1. - x0_sero_50[,n_2014],type="l",xlim=c(6,16),ylim=c(0,1),lwd=2,col="black",las=1,
     ylab="",xlab="",main="")
polygon(c(age_set,rev(age_set)),c(1. - x0_sero_2p5[,n_2014],rev(1. - x0_sero_97p5[,n_2014])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.3),border=NA)
Binom_infer <- binom.bayes(Sero_2014_fit,TrialNum_2014,conf.level=0.95)
points(age_x,Binom_infer$mean,type="l",xlim=c(6,14),ylim=c(0,1),lty=2,lwd=2,xlab="age",ylab="proportion",col="red")
arrows(age_x, Binom_infer$lower, age_x, Binom_infer$upper, 
       length=0.05, angle=90, code=3,col="red")
legend(x="topleft",
       c("Measurements",
         "Model estimation"), 
       lty=c(2,1), 
       pch=c(NA,NA),
       lwd=c(2,2),
       col=c("red","black"),
       cex=1.,bg='white')
dev.off()


# 2015 seroprevalence
n_2015 <- which(Surveilanceyear == 2015)
age_x <- c(1:length(Seroposi_2015))
pdf(paste(fig_save_place,"/Seroposi_compare_2015_All1or0.pdf",sep=""),width = 7.*0.7, height = 5*0.7)
Sero_2015_fit = Seroposi_1_2015 + Seroposi_2_2015 + Seroposi_3_2015 + Seroposi_4_2015 + Seroposi_twice_2015
plot(age_set,1. - x0_sero_50[,n_2015],type="l",xlim=c(10,60),
     ylim=c(0,1),lwd=2,col="black",las=1,main="",ylab="",xlab="")
polygon(c(age_set,rev(age_set)),c(1.-x0_sero_2p5[,n_2015],rev(1.- x0_sero_97p5[,n_2015])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.3),border=NA)
coarse_num <- 5
Seroposi_1_2015_coarse <- seq(coarse_num,80,coarse_num)
TrialNum_2015_coarse <- seq(coarse_num,80,coarse_num)
age_x_coarse <- seq(coarse_num,80,coarse_num) - coarse_num/2
count_ <- 0
for (i in seq(coarse_num,80,coarse_num)){
  count_ <- count_ + 1
  Seroposi_1_2015_coarse[count_] <- sum(Sero_2015_fit[(i-coarse_num+1):i])
  TrialNum_2015_coarse[count_] <- sum(TrialNum_2015[(i-coarse_num+1):i])
}
Binom_infer <- binom.bayes(Seroposi_1_2015_coarse,TrialNum_2015_coarse,conf.level=0.95)
points(age_x_coarse,Binom_infer$mean,type="l",xlim=c(6,14),ylim=c(0,1),lty=2,lwd=2,xlab="age",ylab="proportion",col="red")
arrows(age_x_coarse, Binom_infer$lower, age_x_coarse, Binom_infer$upper, 
       length=0.05, angle=90, code=3,col="red")
legend(x="bottomright",
       c("Measurements",
         "Model estimation"), 
       lty=c(2,1), 
       pch=c(NA,NA),
       lwd=c(2,2),
       col=c("red","black"),
       cex=1.,bg='white')
dev.off()





#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
##### Figure 5           
##### (1) FOI vs proportion of suceptible
##### (2) logistic regression 
##### (3) ROC curve and 
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################

#######################
### Plot FOI vs proportion of suceptible
#######################
for (division_ in c(1,2,3)){
  if (division_ == 1){
    agelist = list(c(1:20),c(1:80))
    labelist = c("Age between 1 and 20","Age between 1 and 80")
  } else if (division_ == 2){
    agelist = list(c(1:4),c(5:9))
    labelist = c("Age between 1 and 4","Age between 5 and 9")
  } else if (division_ == 3){
    agelist = list(c(10:14),c(15:80))
    labelist = c("Age between 10 and 14","Age between 15 and 80")
  }
  
  
  pdf(paste(fig_save_place,"/FOI_vs_susceptible_with_PCC_", division_, ".pdf",sep=""),width = 15.*0.8, height = 4*0.8)
  par(mfrow=c(1,2),
      oma = c(4,3,0,0) + 0.1,
      mar = c(2,2,2,2) + 0.1)
  count_ <- 0
  for (age_list in agelist){
    count_ <- count_ + 1
    
    x0_50_mean <- c()
    mono1_50_mean <- c()
    mono2_50_mean <- c()
    mono3_50_mean <- c()
    mono4_50_mean <- c()
    for (i in 1:length(Surveilanceyear)){ # taking a weighted average (over age) of x0, ..., mono4.
      x0_50_mean[i] <- weighted.mean(x0_50[age_list,i],pop_dist_detail[age_list,i])
      mono1_50_mean[i] <- weighted.mean(mono1_50[age_list,i],pop_dist_detail[age_list,i])
      mono2_50_mean[i] <- weighted.mean(mono2_50[age_list,i],pop_dist_detail[age_list,i])
      mono3_50_mean[i] <- weighted.mean(mono3_50[age_list,i],pop_dist_detail[age_list,i])
      mono4_50_mean[i] <- weighted.mean(mono4_50[age_list,i],pop_dist_detail[age_list,i])
    }
    
    sus_ <- c() ## proportion of susceptible
    for (i in 1:length(case_sum_set)){
      if (serotype_number[i] == 1){
        sus_[i] <- mono4_50_mean[i]+mono3_50_mean[i]+mono2_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 2){
        sus_[i] <- mono4_50_mean[i]+mono3_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 3){
        sus_[i] <- mono4_50_mean[i]+mono2_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 4){
        sus_[i] <- mono3_50_mean[i]+mono2_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      }
    }
    
    FOI_all = FOI_1_50 + FOI_2_50 + FOI_3_50 + FOI_4_50
    dummy_set <- which(serotype_number == 1)
    plot(sus_[dummy_set],FOI_all[dummy_set],col="red",pch=19,xlim=c(0.,1),ylim=c(0,1)
         ,main=labelist[[count_]],cex=2,las=1)
    dummy_set <- which(serotype_number == 2)
    points(sus_[dummy_set],FOI_all[dummy_set],col="blue",pch=19,cex=2)
    dummy_set <- which(serotype_number == 3)
    points(sus_[dummy_set],FOI_all[dummy_set],col="green",pch=19,cex=2)
    dummy_set <- which(serotype_number == 4)
    points(sus_[dummy_set],FOI_all[dummy_set],col="gold",pch=19,cex=2)
    
    mtext("Force of infections",side=2,line=0.8,cex=1.2,outer=TRUE)
    mtext("Proportion of susceptible (including primary and secondary infections)",side=1,line=0.8,cex=1.2,outer=TRUE)
    if (count_ == 1){
      legend(x="bottomleft",
             c("DENV1","DENV2","DENV3","DENV4"), 
             pch=c(19,19,19,19),
             col=c("red","blue","green","gold"),
             ncol=2
             )
    }
    Number_threshold <- 0.
    thres_set <- which(FOI_all>Number_threshold)
    cor_dum <- cor.test(sus_[thres_set],FOI_all[thres_set],method = c("pearson"))
    print(labelist[count_])
    print(cor_dum)
    cor_expect <- cor_dum$estimate
    cor_2p5 <- cor_dum$conf.int[1]
    cor_97p5 <- cor_dum$conf.int[2]
    ### Plot regression line ########
    x=sus_[thres_set]
    y=FOI_all[thres_set]
    lm.out <- lm(y ~ x)
    newx = seq(min(x),max(x),by = 0.05)
    conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",level = 0.95)
    polygon(c(newx,rev(newx)),c(conf_interval[,2],rev(conf_interval[,3])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.3),border=NA)
    abline(lm.out, col="black")
    legend(x="topleft",
           c(paste0("PCC: ",format(cor_expect,digits = 2)), paste0("(95%CI: ",format(cor_2p5,digits = 2)," ~ ",format(cor_97p5,digits = 2),")")), #fill=c("red", "blue","green","gold"), 
           lty=c(1,1), 
           col=c("black","white"),
           cex=1.,bg='white')
  }
  dev.off()
  
  
  
  ###########################################
  ## logistic regression (epidemic or not)
  ###########################################
  Only_primary_susceptible <- 0 ## 1: yes, 0: no ("yes" means taking only susceptibles for primary infections. In the article, "no" is used.)
  if (Only_primary_susceptible==1){
    pdf(paste(fig_save_place,"/casenumber_vs_susceptible_with_logistic_onlyPrimaryInfections_", division_, ".pdf",sep=""),width = 11., height = 7)
  } else {
    pdf(paste(fig_save_place,"/casenumber_vs_susceptible_with_logistic_", division_, ".pdf",sep=""),width = 15.*0.8, height = 4*0.8)
  }
  par(mfrow=c(1,2),
      oma = c(4,3,0,0) + 0.1,
      mar = c(2,2,2,2) + 0.1)
  count_ <- 0
  for (age_list in agelist){
    count_ <- count_ + 1
    x0_50_mean <- c()
    mono1_50_mean <- c()
    mono2_50_mean <- c()
    mono3_50_mean <- c()
    mono4_50_mean <- c()
    for (i in 1:length(Surveilanceyear)){ # taking a weighted average (over age) of x0, ..., mono4.
      x0_50_mean[i] <- weighted.mean(x0_50[age_list,i],pop_dist_detail[age_list,i])
      mono1_50_mean[i] <- weighted.mean(mono1_50[age_list,i],pop_dist_detail[age_list,i])
      mono2_50_mean[i] <- weighted.mean(mono2_50[age_list,i],pop_dist_detail[age_list,i])
      mono3_50_mean[i] <- weighted.mean(mono3_50[age_list,i],pop_dist_detail[age_list,i])
      mono4_50_mean[i] <- weighted.mean(mono4_50[age_list,i],pop_dist_detail[age_list,i])
    }
    
    sus_ <- c() # proportion of susceptibles
    for (i in 1:length(case_sum_set)){
      if (serotype_number[i] == 1){
        sus_[i] <- mono4_50_mean[i]+mono3_50_mean[i]+mono2_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 2){
        sus_[i] <- mono4_50_mean[i]+mono3_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 3){
        sus_[i] <- mono4_50_mean[i]+mono2_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 4){
        sus_[i] <- mono3_50_mean[i]+mono2_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      }
    }
    
    if (Only_primary_susceptible==1){
      sus_ <- x0_50_mean[1:length(case_sum_set)]
    }
    
    Number_threshold <- 300 ## If the number of cases exceeds this value, we consider that year as epidemic (see SI of the article.)
    
    thres_set <- which(case_sum_set >= Number_threshold)
    thres_non_set <- which(case_sum_set < Number_threshold)
    x=sus_[thres_set]
    x_non=sus_[thres_non_set]
    plot(x,rep(1,length(x)),col="red",pch=19,xlim=c(0.,1),ylim=c(0,1),
         main=labelist[[count_]],cex=1.2,las=1)
    points(x_non,rep(0,length(x_non)),col="green",pch=19,cex=1.2)
    
    ### logistic regression #####  
    epione_nonepizero <- sus_
    epione_nonepizero[thres_set] <- 1
    epione_nonepizero[thres_non_set] <- 0
    v1 <-epione_nonepizero
    v2 <-sus_
    liste <-list(v1,v2)
    donnees<-data.frame(liste)
    colnames(donnees) <-c("epi","sus")
    mylogit <- glm(epi ~ sus, data = donnees, family = "binomial")
    xdum <- seq(0,1,0.001)
    predicted_ <- predict(mylogit, list(sus = seq(0,1,0.001)), type="link", se.fit=TRUE)
    se_high = inv.logit(predicted_$fit + (predicted_$se.fit*1.96))
    se_low = inv.logit(predicted_$fit - (predicted_$se.fit*1.96))
    expected = inv.logit(predicted_$fit)
    ydum_97p5 <- se_high
    ydum_2p5 <- se_low
    ydum <- expected
    lines(xdum,expected,lwd=2)
    polygon(c(xdum,rev(xdum)),c(ydum_2p5,rev(ydum_97p5)),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.3),border=NA)
    if (Only_primary_susceptible==1){
      mtext("Proportion of susceptible (including only primary infections)",side=1,line=0.8,cex=1.2,outer=TRUE)
    } else {
      mtext("Proportion of susceptible (including primary and secondary infections)",side=1,line=0.8,cex=1.2,outer=TRUE)
    }
    mtext("Probability of epidemics",side=2,line=0.8,cex=1.2,outer=TRUE)
    if (count_ == 1){
      legend(x="topleft",
             c("Epidemic year (1)",
               "Non-epidemic year (0)",
               "Probability of epidemic (logistic reg.)"), 
             lty=c(NA,NA,1,1), 
             pch=c(19,19,NA,NA),
             col=c("red","green","black"),
             lwd=c(NA,NA,1,2)
             )
    }
  }
  dev.off()
  
  
  
  
  #########################################################
  ##### ROC curve (epidepic or not)  ######################
  #########################################################
  Only_primary_susceptible <- 0 ## 1: yes, 0: no ("yes" means taking only susceptibles for primary infections. In the article, "no" is used.)
  color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') { # making a color bar
    scale = (length(lut)-1)/(max-min)
    points(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale + min
      rect(0.95,y,1,y+1/scale, col=lut[i], border=NA)
      if (i %% 200 == 1){
        text(0.9,y,format(round(i/length(lut), 1), nsmall = 1))
      }
    }
    i <- length(lut)
    y = (i-1)/scale + min
    text(0.9,y,format(round(i/length(lut), 1), nsmall = 1))
    text(0.9,y+0.15,"Threshold")
  }
  
  library(pROC)
  
  if (Only_primary_susceptible==1){
    pdf(paste(fig_save_place,"/ROC_curve_onlyPrimarySusceptible_", division_, ".pdf",sep=""),width = 15., height = 4)
  } else {
    pdf(paste(fig_save_place,"/ROC_curve_", division_, ".pdf",sep=""),width = 12.*0.8, height = 4*0.8)
  }
  par(mfrow=c(1,2),
      oma = c(4,3,0,0) + 0.1,
      mar = c(2,2,2,2) + 0.1)
  count_ <- 0
  for (age_list in agelist){
    count_ <- count_ + 1
    
    x0_50_mean <- c()
    mono1_50_mean <- c()
    mono2_50_mean <- c()
    mono3_50_mean <- c()
    mono4_50_mean <- c()
    for (i in 1:length(Surveilanceyear)){  # taking a weighted average (over age) of x0, ..., mono4.
      x0_50_mean[i] <- weighted.mean(x0_50[age_list,i],pop_dist_detail[age_list,i])
      mono1_50_mean[i] <- weighted.mean(mono1_50[age_list,i],pop_dist_detail[age_list,i])
      mono2_50_mean[i] <- weighted.mean(mono2_50[age_list,i],pop_dist_detail[age_list,i])
      mono3_50_mean[i] <- weighted.mean(mono3_50[age_list,i],pop_dist_detail[age_list,i])
      mono4_50_mean[i] <- weighted.mean(mono4_50[age_list,i],pop_dist_detail[age_list,i])
    }
    
    sus_ <- c()
    for (i in 1:length(case_sum_set)){
      if (serotype_number[i] == 1){
        sus_[i] <- mono4_50_mean[i]+mono3_50_mean[i]+mono2_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 2){
        sus_[i] <- mono4_50_mean[i]+mono3_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 3){
        sus_[i] <- mono4_50_mean[i]+mono2_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      } else if (serotype_number[i] == 4){
        sus_[i] <- mono3_50_mean[i]+mono2_50_mean[i]+mono1_50_mean[i]+x0_50_mean[i]
      }
    }
    
    if (Only_primary_susceptible==1){
      sus_ <- x0_50_mean[1:length(case_sum_set)]
    } 
    
    Number_threshold <- 300 ## If the number of cases exceeds this value, we consider that year as epidemic (see SI of the article.)
    
    thres_set <- which(case_sum_set >= Number_threshold)
    thres_non_set <- which(case_sum_set < Number_threshold)
    x=sus_[thres_set]
    x_non=sus_[thres_non_set]
    epione_nonepizero <- sus_
    epione_nonepizero[thres_set] <- 1
    epione_nonepizero[thres_non_set] <- 0
    v1 <-epione_nonepizero
    v2 <-sus_
    liste <-list(v1,v2)
    donnees<-data.frame(liste)
    colnames(donnees) <-c("epi","sus")
    
    ### ROC (error bars)
    roc_test <- roc(response=donnees$epi,predictor=donnees$sus)
    plot(1 - roc_test$specificities,roc_test$sensitivities,type="l",main=labelist[[count_]])
    range_xx <- seq(0, 1, .05)
    roc_confi <- ci.se(roc_test,specificities = range_xx,conf.level=0.95)
    polygon(c(1 - range_xx,rev(1 - range_xx)),c(roc_confi[,1],rev(roc_confi[,3])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.3),border=NA)
    
    ### ROC 
    Thresholdvalues <- seq(0.,1,0.001)
    datadummy <- donnees$sus
    data_result <- donnees$epi
    ROC_x <- Thresholdvalues
    ROC_y <- Thresholdvalues
    dummy_count <- 0
    for (Th in Thresholdvalues){
      dummy_count <- dummy_count + 1
      N_predict <- which(datadummy < Th)
      P_predict <- which(datadummy >= Th)
      TP <- length(which(data_result[P_predict]==1))
      FN <- length(which(data_result[N_predict]==1))
      FP <- length(which(data_result[P_predict]==0))
      TN <- length(which(data_result[N_predict]==0))
      ROC_x[dummy_count] <- FP/(FP + TN)
      ROC_y[dummy_count] <- TP/(TP + FN)
    }
    
    rbPal <- colorRampPalette(c("purple","blue","light blue","green","light green", "yellow", "orange", "red"))
      # This adds a column of color values based on the y values
    dat_Col <- rbPal(length(Thresholdvalues))[as.numeric(cut(Thresholdvalues,breaks = length(Thresholdvalues)))]
    points(ROC_x,ROC_y,pch=19,col=dat_Col)
    if (count_==1){
      plotrange <- Thresholdvalues
      color.bar(rbPal(length(plotrange)), min=0.02,max=0.72)   
      Roc_distance = ROC_x**2 + (1-ROC_y)**2
      Roc_x_save = ROC_x
      Roc_y_save = ROC_y
      Roc_threshold = Thresholdvalues
    }
    
    mtext("True positive rate",side=2,line=0.8,cex=1.2,outer=TRUE)
    mtext("False positive rate",side=1,line=0.8,cex=1.2,outer=TRUE)
  }
  ###############
  dev.off()
}






## to check the threshold values
min_set_ = which(Roc_distance == min(Roc_distance))
print(Roc_x_save[min_set_])
print(Roc_y_save[min_set_])
print(Roc_threshold[min_set_])

#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
##### Figure 6                     
##### Consistency check using synthetic data
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################

if (which_data == 'Synthetic'){ ## below is only for which_data == 'Synthetic'
  ## check how many percent of the parameters are within the confidence interval?

  ## these parameters are not used in the final version of the inference (so remove them by hand).
  Remove_sequence = c(79,# b_miss parameter we don't use
                      78,# secondary infection parameters we don't use
                      77,# secondary infection parameters we don't use
                      75,# secondary infection parameters we don't use
                      74,# secondary infection parameters we don't use
                      72,# secondary infection parameters we don't use
                      71,# secondary infection parameters we don't use
                      69,# secondary infection parameters we don't use
                      68,# secondary infection parameters we don't use
                      44,# 2015 FOI we don't use
                      45,# 2016 FOI we don't use
                      46,# 2017 FOI we don't use
                      47, # 2018 FOI we don't use
                      52 # Background phi (parameter to determine the shape of tanh function) we don't use (we fix it as a constant.)
  )
  remaining_sequence = c(1:(length(parameter50)))[-Remove_sequence]
  

  #####################
  ## check how many percent of the parameters are in the interval.
  #####################
  count_all = 0
  count_yes = 0
  for (i in remaining_sequence){
    count_all = count_all + 1
    if ((parameter2p5[i] < parameter_correct[i]) & (parameter97p5[i] > parameter_correct[i])){
      count_yes = count_yes + 1
    }
  }
  print(toString(count_yes/count_all))
  
  pdf(paste(fig_save_place,"/Parameter_howmany.pdf",sep=""),width = 5.5*0.8, height = 4*0.8)
  plot(c(1:length(parameter_correct[remaining_sequence])),parameter_correct[remaining_sequence],log="y",
       main=toString(count_yes/count_all),xlab='label of parameters',ylab='parameter values')
  xx = c(1:length(parameter_correct[remaining_sequence]))
  polygon(c(xx,rev(xx)),c(parameter2p5[remaining_sequence],rev(parameter97p5[remaining_sequence])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),border=NA)
  dev.off()
  
  #####################
  ## Reporting rate
  #####################
  # background
  pdf(paste(fig_save_place,"/PHI_background_compare.pdf",sep=""),width = 5.5*0.8, height = 4*0.8)
  par(mfrow=c(1,1))
  dum_ <- which(Surveilanceyear<2015)
  plot(Surveilanceyear[dum_],PHI_background_50[dum_]*PHI_Primary_50[1],type="l",ylim = c(0.005,0.07),xlim=c(1979,2014),xlab = "",ylab="")
  polygon(c(Surveilanceyear[dum_],rev(Surveilanceyear[dum_])),c(PHI_background_2p5[dum_]*PHI_Primary_50[1],rev(PHI_background_97p5[dum_]*PHI_Primary_50[1])),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),border=NA)
  points(Surveilanceyear[dum_],PHI_background_correct[dum_],ylim = c(0.005,0.07),xlim=c(1979,2014),xlab = "Year",ylab="",pch=4)#,ylab="Primary reporting rate of DENV1",main=Model_description[3])
  dev.off()
  
  # each serotype
  pdf(paste(fig_save_place,"/PHI_Primary_Secondary_par_eachDE1_compare.pdf",sep=""),width = 4., height = 4)
  # primary
  plot(primary_x-0.1,PHI_Primary_50/PHI_Primary_50[1],
       xlab = "",ylab="",type="p",xlim = c(0.9,4.1),ylim = c(0.01,5.4),xaxt = "n",col="black")
  points(primary_x-0.1,c(1,PHI_Primary_correct[2:4]),
         xlab = "",ylab="",type="p",xlim = c(0.9,4.1),ylim = c(0.01,6.),xaxt = "n",col="black",pch=4)
  # secondary
  points(c(1,2,3,4)+0.1,PHI_Secondary_v2_50[c(1,4,7,10)]/PHI_Secondary_v2_50[1],
         xlab = "",ylab="",type="p",ylim = c(0.01,6.),xaxt = "n",col="red")
  points(c(1,2,3,4)+0.1,PHI_Secondary_correct[c(1,4,7,10)]/PHI_Secondary_correct[1],
         xlab = "",ylab="",type="p",ylim = c(0.01,6.),xaxt = "n",col="red",pch=4)
  arrows(primary_x-0.1, PHI_Primary_2p5, primary_x-0.1, PHI_Primary_97p5, length=0.05, angle=90, code=3,col="black")
  arrows(c(1,2,3,4)+0.1, PHI_Secondary_v2_2p5[c(1,4,7,10)], c(1,2,3,4)+0.1, PHI_Secondary_v2_97p5[c(1,4,7,10)], length=0.05, angle=90, code=3,col="red")
  axis(1, at=1:4, labels=c("DENV-1","DENV-2","DENV-3","DENV-4"))  
  abline(h = 1, lty = 3)
  legend(x="topleft",
         c("Primary (prediction)","Primary (true)","Secondary (prediction)","Secondary (true)"), 
         pch=c(1,4,1,4),
         col=c("black","black","red","red"),
         cex=0.8,ncol=1,text.col = c('black','black','red','red'),bty = "n")
  dev.off()
  
  # secondary / primary
  pdf(paste(fig_save_place,"/PHI_Secondary_par_priDE1_compare.pdf",sep=""),width = 2.*0.9, height = 4*0.9)
  plot(c(1),PHI_Secondary_50[c(1)],type="p",
       xlab = "",ylab="",xaxt = "n",ylim = c(0.01,5.))
  points(c(1),PHI_Secondary_correct[c(1)],type="p",
         xlab = "",ylab="",xaxt = "n",ylim = c(0.01,5.),pch=4)
  arrows(c(1,2,3,4), PHI_Secondary_2p5[c(1)], c(1,2,3,4), PHI_Secondary_97p5[c(1)], length=0.05, angle=90, code=3)
  abline(h = 1, lty = 2)
  dev.off()

  # age factor
  pdf(paste(fig_save_place,"/Age_dependent_par_5_10_compare.pdf",sep=""),width = 5.*0.8, height = 4*0.8)
  age_fac_50 <- c(Age_k_50[2],1,Age_k_50[3],Age_k_50[4])
  age_fac_2p5 <- c(Age_k_2p5[2],1,Age_k_2p5[3],Age_k_2p5[4])
  age_fac_97p5 <- c(Age_k_97p5[2],1,Age_k_97p5[3],Age_k_97p5[4])
  plot(c(1:4),age_fac_50,xlab = "",type="l",ylab="",ylim = c(0.3,1.2),xaxt = "n",main="")
  age_fac_correct <- c(PHI_crossp_correct[2],1,PHI_crossp_correct[3],PHI_crossp_correct[4])
  points(c(1:4),age_fac_correct,xlab = "",type="p",ylab="",ylim = c(0.3,1.2),xaxt = "n",main="",pch=4)
  polygon(c(c(1:4),rev(c(1:4))),c(age_fac_2p5,rev(age_fac_97p5)),col=rgb(red=0.5,green=0.5,blue=0.5,alpha=0.5),border=NA)
  axis(1, at=1:4, labels=c("1-4","5-9","10-14","15-80"))
  dev.off()
  
  
  # k value in negative binomial distribution
  pdf(paste(fig_save_place,"/k_value_negbin_compare.pdf",sep=""),width = 2.*0.9, height = 4*0.9)
  plot(c(1),Age_k_50[1],main="",ylab="",xlab="",ylim=c(0.55,0.8),xaxt='n')
  points(c(1),PHI_crossp_correct[1],main="",ylab="",xlab="",ylim=c(0.55,0.8),xaxt='n',pch=4)
  arrows(c(1), Age_k_2p5[1], c(1), Age_k_97p5[1], 
         length=0.05, angle=90, code=3,col="black")
  dev.off()
  
  
  
  #####################
  ## Force of infection
  #####################
  colset <- c("red","blue","green","gold")
  pdf(paste(fig_save_place,"/FitFig1_FOI_Immunity_compare.pdf",sep=""),width = 10.*0.9, height = 4*0.9)
  par(mfrow=c(1,1),
      oma = c(4,3,0,0) + 0.1,
      mar = c(2,2,2,2) + 0.1)
  Ymaxvalue <- 1
  plot(Year_plot,FOI_1_50,type="l",
       xlim=c(1979,2014),ylim=c(0,Ymaxvalue),col="red",
       xlab = "",ylab="")
  polygon(c(Year_plot,rev(Year_plot)),c(FOI_1_2p5,rev(FOI_1_97p5)),col=rgb(red=1.,green=0.,blue=0.,alpha=0.3),border=NA)
  points(Year_plot[which(FOI_1_correct!=0)],FOI_1_correct[which(FOI_1_correct!=0)],col="red",pch=4)
  #
  lines(Year_plot,FOI_2_50,type="l",xlim=c(1945,2018),col="blue",
        ylim=c(0,Ymaxvalue),xlab = "",
        ylab="")
  polygon(c(Year_plot,rev(Year_plot)),c(FOI_2_2p5,rev(FOI_2_97p5)),col=rgb(red=0.,green=0.,blue=1.,alpha=0.3),border=NA)
  points(Year_plot[which(FOI_2_correct!=0)],FOI_2_correct[which(FOI_2_correct!=0)],col="blue",pch=4)
  lines(Year_plot,FOI_3_50,type="l",col="green",
        xlim=c(1945,2018),ylim=c(0,Ymaxvalue),xlab = "",
        ylab="")
  polygon(c(Year_plot,rev(Year_plot)),c(FOI_3_2p5,rev(FOI_3_97p5)),col=rgb(red=0.,green=1.,blue=0.,alpha=0.3),border=NA)
  points(Year_plot[which(FOI_3_correct!=0)],FOI_3_correct[which(FOI_3_correct!=0)],col="green",pch=4)
  lines(Year_plot,FOI_4_50,type="l",xlim=c(1945,2018),ylim=c(0,Ymaxvalue),col="gold",
        xlab = "",ylab="")
  polygon(c(Year_plot,rev(Year_plot)),c(FOI_4_2p5,rev(FOI_4_97p5)),col=rgb(red=1.,green=1.,blue=0.,alpha=0.5),border=NA)
  points(Year_plot[which(FOI_4_correct!=0)],FOI_4_correct[which(FOI_4_correct!=0)],col="gold",pch=4)
  legend( x="topleft", 
          legend=c("DENV-1 (prediction)","DENV-1 (true)",
                   "DENV-2 (prediction)","DENV-2 (true)",
                   "DENV-3 (prediction)","DENV-3 (true)",
                   "DENV-4 (prediction)","DENV-4 (true)"), 
          col=c("red","red",
                "blue","blue",
                "green","green",
                "gold","gold"), lwd=1, 
          lty=c(1,NA,
                1,NA,
                1,NA,
                1,NA), ncol=4,
          pch=c(NA,4,
                NA,4,
                NA,4,
                NA,4), merge=FALSE, cex=0.87,
          text.col = c("red","red",
                       "blue","blue",
                       "green","green",
                       "gold","gold") )
  dev.off()
}


