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
  
  DATA_INTEGER(model)     // model switch
  DATA_INTEGER(nset)      // number of sets
  DATA_INTEGER(nlen)      // number of length bins in the model  
  DATA_VECTOR(len)        // length bins [nlen]
  
  DATA_VECTOR(npot_ctl)   // number of control pots fished by set [nset]
  DATA_VECTOR(npot_exp)   // experimental pots 
  
  DATA_MATRIX(ctl_dat)    // control pots: number of fish in each length bin by set [nlen, nset]
  DATA_MATRIX(exp_dat)    // experimental pots
  
  // PARAMETER SECTION ----
  
  // lengths where there is a 10, 50 and 90 percent probability that a fish will
  // be retained in the escape-ring traps
  PARAMETER(log_beta10);
  PARAMETER(log_beta50);
  PARAMETER(log_beta90);
  
  // relative probability of entering a control pot (account for possible
  // difference in the degree to which fish are attracted to escape-ring and to
  // control trap)
  PARAMETER(log_delta);
  
  // Random effects for set [nset]
  PARAMETER_VECTOR(nu);
  
  // Troubleshoot model
  PARAMETER(dummy);

  // MODEL -----
  
  Type beta10 = exp(log_beta10);
  Type beta50 = exp(log_beta50);
  Type beta90 = exp(log_beta90);
  Type delta = exp(log_delta);
  
  Type fixed_beta90 = Type(2) * beta50 - beta10;  // fixed beta90 for 2 parameter logistic curve
  Type alpha = beta50;  // use to define the assymetrical arms of the 3 parameter logistic curve
  
  // Selectivity matrix p = probability that a fish of length i in set j that is
  // caught in an escape-ring pot will be retained in the pot
  matrix<Type> slx(nlen, nset);
  slx.setZero();
  
  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {

      switch(model) {
      
      case 1 : // Logistic with 2 parameters
        if(len(i) <= alpha) {
          slx(i,j) = Type(1) / (Type(1) + exp( (Type(-2) * log(Type(3)) * len(i) - beta50 + nu(j)) / (beta50 - beta10) ));
        }  else {
          slx(i,j) = Type(1) / (Type(1) + exp( (Type(-2) * log(Type(3)) * len(i) - beta50 + nu(j)) / (fixed_beta90 - beta50) ));
        }
        break;

      case 2 : // Logistic with 3 parameters
        if(len(i) <= alpha) {
          slx(i,j) = Type(1) / (Type(1) + exp( (Type(-2) * log(Type(3)) * len(i) - beta50 + nu(j)) / (beta50 - beta10) ));
        }  else {
          slx(i,j) = Type(1) / (Type(1) + exp( (Type(-2) * log(Type(3)) * len(i) - beta50 + nu(j)) / (beta90 - beta50) ));
        }
      break;

      // case 3 : // Dome-shaped
      }
    }
  }
  
  // Ratio of the number of control pots to the number of escape-ring pots in set j
  vector<Type> r(nset);
  r.setZero();

  for (int j = 0; j < nset; j++) {
    r(j) = npot_ctl(j) / npot_exp(j);
  }

  // // Given that a fish is caught and retained in one of the set j pots, the
  // // probability that it is an escape ring pot is phi:
  matrix<Type> phi(nlen,nset);
  phi.setZero();

  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      phi(i,j) = slx(i,j) / (slx(i,j) + delta * r(j));
    }
  }

  // OBJECTIVE FUNCTION -----

  Type nll = 0;

  // Negative log likelihood
  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      nll += exp_dat(i,j) * log(phi(i,j)) + ctl_dat(i,j) * log(Type(1) - phi(i,j));
    }
  }

  //Adjust nll for random effect
  for (int j = 0; j < nset; j++) {
    nll -= dnorm(nu(j), Type(0), Type(1), true);
  }
  
  // Troubleshoot model
  // nll = dummy * dummy;
    
  // REPORT SECTION -----
  
  ADREPORT(beta10);
  ADREPORT(beta50);
  ADREPORT(beta90);
  ADREPORT(slx);
  ADREPORT(phi);
  
  return(nll);
  
}
