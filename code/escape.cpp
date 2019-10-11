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
  
  // Troubleshoot model
  PARAMETER(dummy);  
  
  // lengths where there is a 10, 50 and 90 percent probability that a fish will
  // be retained in the escape-ring traps
  PARAMETER(log_s50);
  PARAMETER(log_slxr); // selection range between 75 and 25 %
  // PARAMETER(log_uslx); // = exp(sel90 - sel50)

  // relative probability of entering a control pot (account for possible
  // difference in the degree to which fish are attracted to escape-ring and to
  // control trap)
  PARAMETER(log_delta);
  
  // Random effects for set [nset]
  PARAMETER_VECTOR(nu);
  
  Type s50 = exp(log_s50);
  Type slxr = exp(log_slxr);
  // Type s10 = s50 - exp(log_lslx);
  // Type s90 = s50 + exp(log_uslx);
  Type delta = exp(log_delta);
  
  // MODEL -----
  
  // Selectivity matrix (slx): probability that a fish of length i in set j that
  // is caught in an escape-ring pot will be retained in the pot
  vector<Type> slx(nlen);
  slx.setZero();

  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
        
        slx(i,j) = Type(1.0) / (Type(1.0) + exp( -Type(2.0) * log(Type(3.0)) * (len(i) - s50 + nu(j)) / slxr));
    }  
  }
// 
  std::cout << "len \n" << len;
  std::cout << "slx \n" << slx;
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
  
  REPORT(s50); 
  REPORT(slxr); 
  // REPORT(s10);
  // REPORT(s90);
  REPORT(slx);
  REPORT(phi);
  
  return(nll);
  
}
