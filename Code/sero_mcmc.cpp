#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// (Year needs to have at longest one year interval: less than one year is fine.)
// Becareful, the matrix starts from 0
// [[Rcpp::export]]
NumericMatrix integrate_lambda_func_Integral(NumericVector lambda,
                                    NumericVector Year_,
                                    int TT, int L) {
  int i,j,jcount;
  double dummy_year, dummy_integ, dummy_year_save, integ_residue, year_residue;
  NumericMatrix integrate_lambda(L,TT);
  for (i = 1; i < TT+1; ++i) {
    jcount = 1;
    year_residue = 0.;
    integ_residue = 0.;
    for (j = 1; j < L+1; ++j) {
//      if (i == 31){
//        printf("i: %d, j: %d, jcount, %d \n", i,j, jcount);
//      }
      dummy_year = year_residue;
      dummy_integ = integ_residue;
      while (dummy_year < 1.){
        dummy_year_save = dummy_year;
        dummy_year += (Year_[i+L-jcount] - Year_[i+L-jcount-1]);
        if (dummy_year > 1.){
          dummy_integ += lambda[i+L-jcount-1] * (1. - dummy_year_save);
          year_residue = dummy_year - 1.;
          integ_residue = lambda[i+L-jcount-1] * year_residue;
          dummy_year = 1.;
          //          if (i == 31){
          //  printf("year_residue %f \n", year_residue);
          //  printf("integ_residue %f \n", integ_residue);
          //}
        } else {
          dummy_integ += lambda[i+L-jcount-1] * (Year_[i+L-jcount] - Year_[i+L-jcount-1]);
          year_residue = 0.;
          integ_residue = 0.;
          //          if (i == 31){
          //  printf("year_residue %f \n", year_residue);
          //  printf("integ_residue %f \n", integ_residue);
          //}
        }
        //        if (i == 31){
        //  printf("year0, year1 : %f, %f \n", Year_[i+L-jcount-1], Year_[i+L-jcount]);
        //  printf("lambda: %f, year %f \n", lambda[i+L-jcount-1], Year_[i+L-jcount-1]);
        //  printf("dummy_integ: %f, dummy_year %f:  \n", dummy_integ, dummy_year);
        //}
        jcount += 1;
      } 
      if (dummy_year != 1.) {printf("Error: %f \n", dummy_year);
        printf("i,j: %d, %d \n", i,j);
        }
      if (j == 1) {
        integrate_lambda(j-1,i-1) = dummy_integ;
      } else {
        integrate_lambda(j-1,i-1) = integrate_lambda(j-2,i-1) + dummy_integ;
      }
      //      if (i == 31){
      //  printf("integrate_lambda(j-1,i-1), %f \n", integrate_lambda(j-1,i-1));
      //  printf("\n");
      // }
    }
  }
  return (integrate_lambda);
}

// the same as integrate_lambda_func_Integral, but starting t is shifted. \int da lambda(t+DeltaT-a).
// Here Delta t is determined as shift_ratio * (duration of the time period that starts  at t)
// [[Rcpp::export]]
NumericMatrix integrate_lambda_func_Integral_SHIFTED(NumericVector lambda,
                                             NumericVector Year_,
                                             int TT, int L, double shift_ratio) {
  int i,j,jcount;
  double dummy_year, dummy_integ, dummy_year_save, integ_residue, year_residue;
  NumericMatrix integrate_lambda(L,TT);
  for (i = 1; i < TT+1; ++i) {
    jcount = 1;
    year_residue = shift_ratio*(Year_[L+i] - Year_[L+i-1]);
    integ_residue = shift_ratio* (Year_[L+i] - Year_[L+i-1])  *lambda[L+i-1];
    for (j = 1; j < L+1; ++j) {
      //      if (i == 31){
      //        printf("i: %d, j: %d, jcount, %d \n", i,j, jcount);
      //      }
      dummy_year = year_residue;
      dummy_integ = integ_residue;
      while (dummy_year < 1.){
        dummy_year_save = dummy_year;
        dummy_year += (Year_[i+L-jcount] - Year_[i+L-jcount-1]);
        if (dummy_year > 1.){
          dummy_integ += lambda[i+L-jcount-1] * (1. - dummy_year_save);
          year_residue = dummy_year - 1.;
          integ_residue = lambda[i+L-jcount-1] * year_residue;
          dummy_year = 1.;
          //          if (i == 31){
          //  printf("year_residue %f \n", year_residue);
          //  printf("integ_residue %f \n", integ_residue);
          //}
        } else {
          dummy_integ += lambda[i+L-jcount-1] * (Year_[i+L-jcount] - Year_[i+L-jcount-1]);
          year_residue = 0.;
          integ_residue = 0.;
          //          if (i == 31){
          //  printf("year_residue %f \n", year_residue);
          //  printf("integ_residue %f \n", integ_residue);
          //}
        }
        //        if (i == 31){
        //  printf("year0, year1 : %f, %f \n", Year_[i+L-jcount-1], Year_[i+L-jcount]);
        //  printf("lambda: %f, year %f \n", lambda[i+L-jcount-1], Year_[i+L-jcount-1]);
        //  printf("dummy_integ: %f, dummy_year %f:  \n", dummy_integ, dummy_year);
        //}
        jcount += 1;
      } 
      if (dummy_year != 1.) {
        if (i!= 47){
          printf("Error: %f \n", dummy_year);
          printf("i,j: %d, %d \n", i,j);
        }
      }
      if (j == 1) {
        integrate_lambda(j-1,i-1) = dummy_integ;
      } else {
        integrate_lambda(j-1,i-1) = integrate_lambda(j-2,i-1) + dummy_integ;
      }
      //      if (i == 31){
      //  printf("integrate_lambda(j-1,i-1), %f \n", integrate_lambda(j-1,i-1));
      //  printf("\n");
      // }
    }
  }
  return (integrate_lambda);
}



// [[Rcpp::export]]
NumericMatrix poi_rate_func(NumericVector lambda1,
                            NumericVector lambda2,
                            NumericVector lambda3,
                            NumericVector lambda4,
                            NumericMatrix x0,
                            NumericMatrix mono1,
                            NumericMatrix mono2,
                            NumericMatrix mono3,
                            NumericMatrix mono4,
                            NumericVector PHI_background,
                            NumericVector PHI_Primary,
                            NumericVector PHI_Secondary,
                            NumericVector serotype_number,
                            int TT, int L,
                            NumericVector Surveilanceyear_roof,
                            double ReAgeFac2,
                            double ReAgeFac3,
                            double ReAgeFac4) {
  int i,j;
  NumericMatrix poi_rate(L,TT);
  double delta_T;
  double age_fac;
  for (i = 1; i < TT+1; ++i) {
    delta_T = (Surveilanceyear_roof[i] - Surveilanceyear_roof[i-1]);
    if (serotype_number[i-1] == 1) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[1-1]*lambda1[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono2(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[1-1]*lambda1[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono3(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[2-1]*lambda1[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono4(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[3-1]*lambda1[i+L-1] * delta_T * age_fac;
      }
    } else if (serotype_number[i-1] == 2) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[2-1]*lambda2[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono1(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[4-1]*lambda2[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono3(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[5-1]*lambda2[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono4(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[6-1]*lambda2[i+L-1] * delta_T * age_fac;
      }
    } else if (serotype_number[i-1] == 3) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[3-1]*lambda3[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono1(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[7-1]*lambda3[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono2(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[8-1]*lambda3[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono4(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[9-1]*lambda3[i+L-1] * delta_T * age_fac;
      }
    } else if (serotype_number[i-1] == 4) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[4-1]*lambda4[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono1(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[10-1]*lambda4[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono2(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[11-1]*lambda4[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono3(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[12-1]*lambda4[i+L-1] * delta_T * age_fac;
      }
    }
  }
  return (poi_rate);
}


// [[Rcpp::export]]
NumericMatrix poi_rate_func_with_rate(NumericVector lambda1,
                            NumericVector lambda2,
                            NumericVector lambda3,
                            NumericVector lambda4,
                            NumericMatrix x0,
                            NumericMatrix mono1,
                            NumericMatrix mono2,
                            NumericMatrix mono3,
                            NumericMatrix mono4,
                            NumericVector PHI_background,
                            NumericVector PHI_Primary,
                            NumericVector PHI_Secondary,
                            NumericVector serotype_number,
                            int TT, int L,
                            NumericVector Surveilanceyear_roof,
                            double ReAgeFac2,
                            double ReAgeFac3,
                            double ReAgeFac4,
                            double rate_) {
  int i,j;
  NumericMatrix poi_rate(L,TT);
  double delta_T;
  double age_fac;
  for (i = 1; i < TT+1; ++i) {
    delta_T = (Surveilanceyear_roof[i] - Surveilanceyear_roof[i-1]) * rate_;
    if (serotype_number[i-1] == 1) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[1-1]*lambda1[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono2(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[1-1]*lambda1[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono3(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[2-1]*lambda1[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono4(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[3-1]*lambda1[i+L-1] * delta_T * age_fac;
      }
    } else if (serotype_number[i-1] == 2) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[2-1]*lambda2[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono1(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[4-1]*lambda2[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono3(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[5-1]*lambda2[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono4(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[6-1]*lambda2[i+L-1] * delta_T * age_fac;
      }
    } else if (serotype_number[i-1] == 3) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[3-1]*lambda3[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono1(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[7-1]*lambda3[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono2(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[8-1]*lambda3[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono4(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[9-1]*lambda3[i+L-1] * delta_T * age_fac;
      }
    } else if (serotype_number[i-1] == 4) {
      for (j = 1; j < L+1; ++j) {
        if (j >= 1 && j < 5){
          age_fac = ReAgeFac2;
        } else if (j >= 5 && j < 10){
          age_fac = 1.;
        } else if (j >= 10 && j < 15){
          age_fac = ReAgeFac3;
        } else {
          age_fac = ReAgeFac4;
        } 
        poi_rate(j-1,i-1)=x0(j-1,i-1)* PHI_background[i-1]*PHI_Primary[4-1]*lambda4[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono1(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[10-1]*lambda4[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono2(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[11-1]*lambda4[i+L-1] * delta_T * age_fac;
        poi_rate(j-1,i-1) += mono3(j-1,i-1)* PHI_background[i-1]*PHI_Secondary[12-1]*lambda4[i+L-1] * delta_T * age_fac;
      }
    }
  }
  return (poi_rate);
}




// [[Rcpp::export]]
NumericMatrix Substitution_lambda(NumericVector x,
                                  NumericVector serotype_number,
                                  NumericVector Year_,
                                  int lambda_length,
                                  int TT, int Year_length,
                                  int index1944,
                                  int index1978) {
  int i;
  NumericMatrix lambda_matrix(4,Year_length);
  for (i = 1; i < index1944+1; ++i) {
    lambda_matrix(4-1,i-1) = 0.;
    lambda_matrix(3-1,i-1) = 0.;
    lambda_matrix(2-1,i-1) = 0.;
    lambda_matrix(1-1,i-1) = 0.;
  }

  for (i = index1944+1; i < index1978+1; ++i) {
    lambda_matrix(4-1,i-1) = 0.;
    lambda_matrix(3-1,i-1) = 0.;
    lambda_matrix(2-1,i-1) = 0.;
    lambda_matrix(1-1,i-1) = 0.;
  }
  
  lambda_matrix(1-1,index1944-1) = x[lambda_length-1];
  lambda_matrix(3-1,index1944+20-1) = x[lambda_length-2];
  lambda_matrix(3-1,index1944+21-1) = x[lambda_length-2];
  lambda_matrix(3-1,index1944+22-1) = x[lambda_length-2];
  lambda_matrix(3-1,index1944+23-1) = x[lambda_length-2];
  lambda_matrix(3-1,index1944+24-1) = x[lambda_length-2];
  lambda_matrix(3-1,index1944+25-1) = x[lambda_length-2];
  lambda_matrix(2-1,index1944+27-1) = x[lambda_length-3];
  lambda_matrix(2-1,index1944+28-1) = x[lambda_length-3];
  lambda_matrix(2-1,index1944+29-1) = x[lambda_length-3];
  lambda_matrix(2-1,index1944+30-1) = x[lambda_length-3];
  lambda_matrix(2-1,index1944+31-1) = x[lambda_length-3];
  lambda_matrix(1-1,index1944+32-1) = x[lambda_length-4];
  lambda_matrix(1-1,index1944+33-1) = x[lambda_length-4];
  lambda_matrix(1-1,index1944+34-1) = x[lambda_length-4];
  

  for (i = 1; i < TT+1; ++i){
    if (serotype_number[i-1]== 1){
      lambda_matrix(1-1,index1978+i-1) = x[i-1];
      lambda_matrix(2-1,index1978+i-1) = 0.;
      lambda_matrix(3-1,index1978+i-1) = 0.;
      lambda_matrix(4-1,index1978+i-1) = 0.;
    } else if (serotype_number[i-1] == 2){
      lambda_matrix(2-1,index1978+i-1) = x[i-1];
      lambda_matrix(1-1,index1978+i-1) = 0.;
      lambda_matrix(3-1,index1978+i-1) = 0.;
      lambda_matrix(4-1,index1978+i-1) = 0.;
    } else if (serotype_number[i-1] == 3){
      lambda_matrix(3-1,index1978+i-1) = x[i-1];
      lambda_matrix(1-1,index1978+i-1) = 0.;
      lambda_matrix(2-1,index1978+i-1) = 0.;
      lambda_matrix(4-1,index1978+i-1) = 0.;
    } else if (serotype_number[i-1] == 4){
      lambda_matrix(4-1,index1978+i-1) = x[i-1];
      lambda_matrix(1-1,index1978+i-1) = 0.;
      lambda_matrix(2-1,index1978+i-1) = 0.;
      lambda_matrix(3-1,index1978+i-1) = 0.;
    } else if (serotype_number[i-1] == 0){
      lambda_matrix(4-1,index1978+i-1) = x[TT-1];
      lambda_matrix(3-1,index1978+i-1) = x[TT-2];
      lambda_matrix(2-1,index1978+i-1) = x[TT-3];
      lambda_matrix(1-1,index1978+i-1) = x[TT-4];
    } 
  }
  return (lambda_matrix);
}


// [[Rcpp::export]]
NumericMatrix poi_rate_cohort_func(NumericMatrix poi_rate,
                                   NumericVector roof_age,
                                   NumericMatrix C,
                                   int TT, int AA) {
  int i,j,k;
  double dummy;
  NumericMatrix poi_rate_cohort(AA,TT);
  for (i = 1; i < TT+1; ++i) {
    for (j = 1; j < AA+1; ++j) {
      if (j==AA){    	
        dummy = 0.;
        for (k = (roof_age[j-1]); k < (roof_age[j]); ++k) {
          dummy += poi_rate(k-1,i-1);
        }
        poi_rate_cohort(j-1,i-1) = dummy * C(j-1,i-1)/(roof_age[j]-roof_age[j-1]);
      } else{
        dummy = 0.;
        for (k = (roof_age[j-1]); k < (roof_age[j]); ++k) {
          dummy += poi_rate(k-1,i-1);
        }
        poi_rate_cohort(j-1,i-1) = dummy * C(j-1,i-1)/(roof_age[j]-roof_age[j-1]);
      }
    }
  }
  return (poi_rate_cohort);
}



// [[Rcpp::export]]
List mono_tilde_calculate(int L,int TT,NumericMatrix x0,
                          NumericMatrix mono1,
                          NumericMatrix mono2,
                          NumericMatrix mono3,
                          NumericMatrix mono4,
                          NumericVector which_Change_serotype_1,
                          NumericVector How_many_within_kkyears,
                          NumericMatrix dummy_sub_lambda,
                          NumericMatrix Time_interval_set,
                          NumericVector serotype_number_extended,
                          int AddBegining
) {
  
  NumericMatrix  mono1_tilde(L,TT);
  NumericMatrix  mono2_tilde(L,TT);
  NumericMatrix  mono3_tilde(L,TT);
  NumericMatrix  mono4_tilde(L,TT);
  NumericMatrix  x0_stepset(max(How_many_within_kkyears)+1,L);
  NumericVector  delx_stepset(L);
  NumericMatrix monotilde_seroset_dummy(4,L);
  List monotilde_seroset = List::create(mono1,mono2,mono3,mono4);
  double lambda_sum;
  int i,idum,i_year,j,jtime,serotype_dummy;
  int which_Change_serotype_1_num = which_Change_serotype_1.size();
  for (i = 1; i < L+1; ++i){
    for (j = 1; j < TT+1; ++j){
      mono1_tilde(i-1,j-1) = mono1(i-1,j-1);
      mono2_tilde(i-1,j-1) = mono2(i-1,j-1);
      mono3_tilde(i-1,j-1) = mono3(i-1,j-1);
      mono4_tilde(i-1,j-1) = mono4(i-1,j-1);
    }
  }
  for (idum = 1; idum < which_Change_serotype_1_num+1; ++idum){
    i = which_Change_serotype_1[idum-1];
    i_year = i + L;
    //    mono_seroset <- rbind(mono1[,i],mono2[,i],mono3[,i],mono4[,i])
    //    monotilde_seroset <- matrix(0, ncol=L, nrow=4)
    for (j = 1;j<L+1; ++j){
      x0_stepset(1-1,j-1) = x0(j-1,i-1);
    }
    for (jtime = 1; jtime < How_many_within_kkyears[i-1]+1; ++jtime){
      lambda_sum = dummy_sub_lambda(0,i_year-jtime-1) + dummy_sub_lambda(1,i_year-jtime-1) + dummy_sub_lambda(2,i_year-jtime-1) + dummy_sub_lambda(3,i_year-jtime-1);
      for (j = 1;j<L+1; ++j){
        x0_stepset(jtime+1-1,j-1) = x0_stepset(jtime-1,j-1) * exp(lambda_sum*Time_interval_set(i-1,jtime-1));
      }
      for (j = 1;j<L+1; ++j){
        if (x0_stepset(jtime+1-1,j-1)>1){
          x0_stepset(jtime+1-1,j-1) = 1;
        }
      }
    }
    for (serotype_dummy = 1;serotype_dummy<4+1; ++serotype_dummy){
      for (j = 1;j<L+1; ++j){
        delx_stepset[j-1] = 0;
      }
      for (jtime = 1;jtime<How_many_within_kkyears(i-1)+1; ++jtime){
        if (serotype_number_extended(i+AddBegining - jtime-1) == serotype_dummy){
          for (j = 1;j<L+1; ++j){
            delx_stepset[j-1] = delx_stepset[j-1] - (x0_stepset(jtime-1,j-1) - x0_stepset(jtime+1-1,j-1));
          }
        }
      }
      
      if (serotype_dummy==1){
        for (j = 1;j<L+1; ++j){
          mono1_tilde(j-1,i-1) = mono1(j-1,i-1) - delx_stepset(j-1);
          if (mono1_tilde(j-1,i-1)<0) {mono1_tilde(j-1,i-1)=0;}
        }
      } else if (serotype_dummy==2){
        for (j = 1;j<L+1; ++j){
          mono2_tilde(j-1,i-1) = mono2(j-1,i-1) - delx_stepset(j-1);
          if (mono2_tilde(j-1,i-1)<0) {mono2_tilde(j-1,i-1)=0;}
        }
      } else if (serotype_dummy==3){
        for (j = 1;j<L+1; ++j){
          mono3_tilde(j-1,i-1) = mono3(j-1,i-1) - delx_stepset(j-1);
          if (mono3_tilde(j-1,i-1)<0) {mono3_tilde(j-1,i-1)=0;}
        }
      } else if (serotype_dummy==4){
        for (j = 1;j<L+1; ++j){
          mono4_tilde(j-1,i-1) = mono4(j-1,i-1) - delx_stepset(j-1);
          if (mono4_tilde(j-1,i-1)<0) {mono4_tilde(j-1,i-1)=0;}
        }
      }
    }
  }
  monotilde_seroset[0] = mono1_tilde;
  monotilde_seroset[1] = mono2_tilde;
  monotilde_seroset[2] = mono3_tilde;
  monotilde_seroset[3] = mono4_tilde;
  return (monotilde_seroset);
}

