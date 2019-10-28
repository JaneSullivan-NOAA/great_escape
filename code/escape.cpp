// Escape ring model based on Haist et al. 2004 Appendix N
// jane.sullivan1@alaska.gov
// Last updated 2019-10-09

#include <TMB.hpp>
#include <numeric>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA SECTION ----
  
  DATA_INTEGER(slx_type)      // model switch (1 = symmetric selectivity, 2 = asymmetrical)
  DATA_INTEGER(nset)          // number of sets
  DATA_INTEGER(nlen)          // number of length bins in the model  
  DATA_VECTOR(len)            // length bins [nlen]
  
  DATA_VECTOR(npot_ctl)       // number of control pots fished by set [nset]
  DATA_VECTOR(npot_exp)       // experimental pots 
  
  DATA_MATRIX(ctl_dat)        // control pots: number of fish in each length bin by set [nlen, nset]
  DATA_MATRIX(exp_dat)        // experimental pots
  
  DATA_INTEGER(mu_log_s90)    // prior mu for log_s90
  DATA_INTEGER(mu_log_s10)    // prior mu for log_s10
  DATA_INTEGER(sig_log_s90)   // prior sigma for log_s90 
  DATA_INTEGER(sig_log_s10)   // prior sigma for log_s10
  
  // PARAMETER SECTION ----
  
  // Troubleshoot model
  PARAMETER(dummy);  
  
  // lengths where there is a 10, 50 and 90 percent probability that a fish will
  // be retained in the escape-ring traps
  PARAMETER(log_s50);
  PARAMETER(log_s90);
  PARAMETER(log_s10);

  // relative probability of entering a control pot (account for possible
  // difference in the degree to which fish are attracted to escape-ring and to
  // control trap)
  PARAMETER(log_delta);
  
  // Random effects for set [nset]
  PARAMETER_VECTOR(nu);
  
  Type s50 = exp(log_s50);
  Type s90 = exp(log_s90);
  Type s10 = exp(log_s10);

  Type delta = exp(log_delta);
  
  // MODEL -----
  
  // Prior on log_s90 and log_s10
  Type prior_log_s90 = 0;
  Type prior_log_s10 = 0;

  prior_log_s90 += Type(0.5) * square(log_s90 - mu_log_s90) / square(sig_log_s90);
  prior_log_s10 += Type(0.5) * square(log_s10 - mu_log_s10) / square(sig_log_s10);

  // Selectivity matrix (slx): probability that a fish of length i in set j that
  // is caught in an escape-ring pot will be retained in the pot
  matrix<Type> slx(nlen,nset);
  slx.setZero();
  
  // Fitted values
  vector<Type> fit_slx(nlen);
  fit_slx.setZero();

  Type delta_slx = 0;
  
  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
        
        switch(slx_type) {
        
        case 1 : // Logistic with 2 parameters
          
          if (len(i) <= s50) 
            delta_slx = s50 - (Type(2) * s50 - s90);
          else 
            delta_slx = s90 - s50;
          
          break;
          
        case 2 : // Logistic with 3 parameters
          
          if (len(i) <= s50) 
            delta_slx = s50 - s10;
          else 
            delta_slx = s90 - s50;
          
          break;
        }
      slx(i,j) = Type(1) / (Type(1) + exp( Type(-2) * log(Type(3)) * ((len(i) - s50 + nu(j)) / delta_slx)));
    }
    fit_slx(i) = Type(1) / (Type(1) + exp( Type(-2) * log(Type(3)) * ((len(i) - s50) / delta_slx)));
  }

    // std::cout << "len \n" << len;
  // std::cout << "slx \n" << slx;
  // 
  // Ratio (r) of the number of control pots to the number of escape-ring pots in set j:
  vector<Type> r(nset);
  r.setZero();

  for (int j = 0; j < nset; j++) {
    r(j) = npot_ctl(j) / npot_exp(j);
  }

  // Given that a fish is caught and retained in one of the set j pots, the
  // probability that it is an escape ring pot is phi:
  matrix<Type> phi(nlen,nset);
  phi.setZero();

  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      phi(i,j) = slx(i,j) / (slx(i,j) + delta * r(j));
    }
  }
  // std::cout << phi << "\n phi";
  
  // OBJECTIVE FUNCTION -----

  Type nll = 0;             // Negative log likelihood
  nll += prior_log_s90;     // Add prior on log_s90
  
  if (slx_type == 2) {
    nll += prior_log_s10;   // Only add log_s10 prior if estimating 3rd parameter 
  }
  
  // Negative log likelihood
  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      if(ctl_dat(i,j) + exp_dat(i,j) > 0)
      nll -= exp_dat(i,j) * log(phi(i,j)) + ctl_dat(i,j) * log(Type(1) - phi(i,j));
    }
  }

  //Adjust nll for random effect
  for (int j = 0; j < nset; j++) {
    nll -= dnorm(nu(j), Type(0), Type(1), true);
  }
  
  // Troubleshoot model
  // nll = dummy * dummy;
    
  // REPORT SECTION -----
  
  REPORT(fit_slx);
  REPORT(s50); 
  REPORT(s10);
  REPORT(s90);
  REPORT(slx);
  REPORT(phi);
  REPORT(prior_log_s90);
  REPORT(prior_log_s10);
  
  return(nll);
  
}
