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
  DATA_VECTOR(fitted_len)     // length vector for fitted values
  
  DATA_VECTOR(npot_ctl)       // number of control pots fished by set [nset]
  DATA_VECTOR(npot_exp)       // experimental pots 
  
  DATA_MATRIX(ctl_dat)        // control pots: number of fish in each length bin by set [nlen, nset]
  DATA_MATRIX(exp_dat)        // experimental pots
  
  // PARAMETER SECTION ----
  
  // Troubleshoot model
  PARAMETER(dummy);
  
  // Selectivity curve parameters
  PARAMETER(alpha);
  PARAMETER(beta);
  // PARAMETER(gamma);
  
  // relative probability of entering a control pot (account for possible
  // difference in the degree to which fish are attracted to escape-ring and to
  // control trap)
  PARAMETER(log_delta);
  
  // Random effects for set [nset]
  PARAMETER_VECTOR(nu);
  
  Type delta = exp(log_delta);
  
  // MODEL -----
  
  // Selectivity matrix (slx): probability that a fish of length i in set j that
  // is caught in an escape-ring pot will be retained in the pot
  matrix<Type> slx(nlen,nset);
  slx.setZero();
  
  Type delta_slx = 0;
  
  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      slx(i,j) = Type(1) / (Type(1) + exp(alpha + beta * len(i) + nu(j)));
    }
  }
  // std::cout << "len \n" << len;
  // std::cout << "slx \n" << slx;
  
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
  
  // Fitted values
  int n_fitted = fitted_len.size();
  vector<Type> fit_slx(n_fitted);
  fit_slx.setZero();
  
  for (int i = 0; i < n_fitted; i++) {
    fit_slx(i) = Type(1) / (Type(1) + exp(alpha + beta * fitted_len(i)));
  }
  
  vector<Type> fit_phi(n_fitted);
  fit_phi.setZero();
  
  for (int i = 0; i < n_fitted; i++) {
    fit_phi(i) = fit_slx(i) / (fit_slx(i) + delta);
  }
  
  // OBJECTIVE FUNCTION -----
  
  // Negative log likelihood
  Type nll = 0;
  
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
  REPORT(alpha);
  REPORT(beta);
  REPORT(slx);
  REPORT(phi);
  REPORT(fit_phi);
  REPORT(r);
  
  return(nll);
  
}
