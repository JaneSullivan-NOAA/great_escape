// Escape ring model inspired by Arana et al., Fisheries Research (2011)
// jane.sullivan1@alaska.gov
// Last updated 2019-11-15

#include <TMB.hpp>
#include <numeric>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA SECTION ----
  
  DATA_INTEGER(nset)          // number of sets [i]
  DATA_INTEGER(nlen)          // number of length bins in the model [j]
  DATA_INTEGER(ntrt)          // number of experimental treatments [k]
  DATA_VECTOR(len)            // length bins [j]
  DATA_VECTOR(fit_len)        // length vector for fitted values

  DATA_VECTOR(npot_ctl)       // number of control pots sampled by set [j]
  DATA_MATRIX(npot_trt)       // number of experimental treatment pots sampled by set [j,k]
  
  DATA_MATRIX(ctl_dat)        // control pots: number of fish in each length bin by set [i,j]
  DATA_ARRAY(trt_dat)         // experimental pots: number of fish in each length bin by set and treatment [i,j,k]
  
  DATA_INTEGER(prior_type)    // Priors based on theoretical curves to consrain selectivity. 0 = normal, 1 = beta
  DATA_IVECTOR(s0_index)      // Index of length where theoretical curves were 0 and 1 for each treatment [k]
  DATA_IVECTOR(s100_index)  
  DATA_SCALAR(sigma_s0)       // Sigma for normal prior
  DATA_SCALAR(sigma_s100)     
  DATA_SCALAR(s0_alpha)       // shape 1 parameter for beta prior
  DATA_SCALAR(s0_beta)        // shape 2 parameter for beta prior
  DATA_SCALAR(s100_alpha)      
  DATA_SCALAR(s100_beta)       
  
  // PARAMETER SECTION ----
  
  // Troubleshoot model
  PARAMETER(dummy);
  
  // Selectivity curve parameters for each experimental treatmnet [k]
  PARAMETER_VECTOR(s50);
  PARAMETER_VECTOR(slp);     
  
  // Relative probability of entering a control pot (account for possible
  // difference in the degree to which fish are attracted to escape-ring and to
  // control trap). This one parameter is shared between treatments.
  PARAMETER(log_delta);
  Type delta = exp(log_delta);  
  
  // Random effects for set [j], shared between treatments.
  PARAMETER_VECTOR(nu);
  
  // MODEL -----
  
  // i = length bin
  // j = set
  // k = treatment
  
  // Logistic selectivity model: probability that a fish of length i is caught
  // in escape-ring treatment pot k will be retained in the pot
  array<Type> slx(nlen,nset,ntrt);
  slx.setZero();

  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      for (int k = 0; k < ntrt; k++) {
        slx(i,j,k) = Type(1) / (Type(1) + exp(Type(-1.0) * slp(k) * (len(i) - s50(k)) + nu(j)));
      }
    }
  }
  // std::cout << "len \n" << len;
  // std::cout << "slx \n" << slx;
  
  // Efficiency ratio (r): ratio of control pots to escape-ring pots for
  // treatment k in set j times the delta, with a random effect at the set level
  // (nu):
  matrix<Type> r(nset,ntrt);
  r.setZero();

  for (int j = 0; j < nset; j++) {
    for (int k = 0; k < ntrt; k++) {
      r(j,k) = npot_ctl(j) / (npot_trt(j,k) * delta);
    }
  }

  // Given that a fish is caught and retained in one of the set j pots, the
  // probability that it is an escape ring pot is phi:
  array<Type> phi(nlen,nset,ntrt);
  phi.setZero();

  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      for (int k = 0; k < ntrt; k++) {
        phi(i,j,k) = slx(i,j,k) / (slx(i,j,k) + r(j,k));
      }
    }
  }
  // std::cout << phi << "\n phi";

  // FITTED VALUES

  // Fitted selectivity  
  matrix<Type> fit_slx(nlen,ntrt);
  fit_slx.setZero();
  
  for (int i = 0; i < nlen; i++) {
    for (int k = 0; k < ntrt; k++) {
      fit_slx(i,k) = Type(1) / (Type(1) + exp(Type(-1.0) * slp(k) * (len(i) - s50(k)))); 
    }
  }
  
  // Fit selectivity across a wider range of lengths for mcmc
  int nfit = fit_len.size();
  matrix<Type> full_slx(nfit,ntrt);
  full_slx.setZero();

  for (int i = 0; i < nfit; i++) {
    for (int k = 0; k < ntrt; k++) {
      full_slx(i,k) = Type(1) / (Type(1) + exp(Type(-1.0) * slp(k) * (fit_len(i) - s50(k)))); 
    }
  }
  
  // Fitted phi (probability of capture) 
  matrix<Type> fit_phi(nlen,ntrt);
  fit_phi.setZero();
  
  for (int i = 0; i < nlen; i++) {
    for (int k = 0; k < ntrt; k++) {
    fit_phi(i,k) = fit_slx(i,k) / (fit_slx(i,k) + delta);
    }
  }

  // OBJECTIVE FUNCTION -----

  // // Negative log likelihood
  Type nll = 0;
  
  for (int i = 0; i < nlen; i++) {
    for (int j = 0; j < nset; j++) {
      for (int k = 0; k < ntrt; k++) {
        if(ctl_dat(i,j) + trt_dat(i,j,k) > 0)
          nll -= trt_dat(i,j,k) * log(phi(i,j,k)) + ctl_dat(i,j) * log(Type(1) - phi(i,j,k));
      }
    }
  }

  //Adjust nll for random effect
  Type set_effect = 0;

  for (int j = 0; j < nset; j++) {
    set_effect -= dnorm(nu(j), Type(0), Type(1), true);
  }
  nll += set_effect;

  // Priors to constrain selectivity curve at 0 and 100
  
  // Initialize values
  Type s0 = 0;
  Type prior_s0 = 0;
  Type s100 = 0;
  Type prior_s100 = 0;
  
  for (int j = 0; j < nset; j++) {
    for (int k = 0; k < ntrt; k++) {
      
      // Extract selectivity proportions by prior lengths at 0 and 100% selected
      // by treatment
      s0 = fit_slx(s0_index(k),k);
      s100 = fit_slx(s100_index(k),k);
      
      if (prior_type == 0) { // Normal
        prior_s0 += Type(0.5) * square(s0 - Type(0)) / square(sigma_s0);
        prior_s100 += Type(0.5) * square(s100 - Type(1)) / square(sigma_s100);
      }
      
      if (prior_type == 1) { // Beta
        prior_s0 += dbeta(s0, s0_alpha, s0_beta, false);
        prior_s100 += dbeta(s100, s100_alpha,s100_beta, false);
      }
    }
  }
  nll += prior_s0;
  nll += prior_s100;
  
  // Troubleshoot model
  // nll = dummy * dummy;
  
  // REPORT SECTION -----
  
  REPORT(full_slx);  
  REPORT(fit_slx);
  REPORT(slx);
  REPORT(fit_phi);  
  REPORT(phi);
  REPORT(set_effect);
  REPORT(prior_s0);
  REPORT(prior_s100);
  REPORT(nll);
  
  return(nll);
  
}
