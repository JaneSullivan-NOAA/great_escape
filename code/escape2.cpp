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
  
  DATA_INTEGER(slx_type)      // model switch (1 = symmetric selectivity, 2 = asymmetrical)
  DATA_INTEGER(nset)          // number of sets [i]
  DATA_INTEGER(nlen)          // number of length bins in the model [j]
  DATA_INTEGER(ntrt)          // number of experimental treatments [k]
  DATA_VECTOR(len)            // length bins [j]
  DATA_VECTOR(fit_len)        // length vector for fitted values
  DATA_VECTOR(trt)            // escape ring treatment diameters [k]
  
  DATA_VECTOR(npot_ctl)       // number of control pots sampled by set [j]
  DATA_MATRIX(npot_trt)       // number of experimental treatment pots sampled by set [j,k]
  
  DATA_MATRIX(ctl_dat)        // control pots: number of fish in each length bin by set [i,j]
  DATA_ARRAY(trt_dat)         // experimental pots: number of fish in each length bin by set and treatment [i,j,k]
  
  DATA_VECTOR(theor_s50)      // theoretical s50's [k]
  DATA_VECTOR(theor_slp)      // theoretical slope of the logistic curve [k]
  DATA_SCALAR(wt_s50)        // weight for penalty on s50
  DATA_SCALAR(wt_slp)        // weight for penalty on slp
           
  DATA_SCALAR(sigma_s0)      // Priors to consrain the selectivity curves at slx = 0 and 100
  DATA_SCALAR(sigma_s100)    
  
  // PARAMETER SECTION ----
  
  // Troubleshoot model
  PARAMETER(dummy);
  
  // Selectivity curve parameters for each experimental treatmnet [k]
  PARAMETER_VECTOR(s50);
  PARAMETER_VECTOR(slp);     
  
  // relative probability of entering a control pot (account for possible
  // difference in the degree to which fish are attracted to escape-ring and to
  // control trap). This one parameter is shared between treatments.
  PARAMETER(log_delta);
  Type delta = exp(log_delta);  
  
  // Linear regression coefficients governing the relationship between length at
  // 50% selectivity (s50) and escape ring diameter (trt): s50 ~ trt
  PARAMETER(a1);
  PARAMETER(b1);  
  // Type(b1 = exp(log_b1));  // estimate in log space to assume positive slope? not sure this is necessary
  
  // Linear regression coefficients governing the relationship between the slope
  // of the logistic selectivity curve (slx_slp) and escape ring diameter (trt):
  // slx_slp ~ trt
  PARAMETER(a2);
  PARAMETER(b2);  
  
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
  int nfit = fit_len.size();
  matrix<Type> fit_slx(nfit,ntrt);
  fit_slx.setZero();

  for (int i = 0; i < nfit; i++) {
    for (int k = 0; k < ntrt; k++) {
      // fit_slx(i,k) = Type(1) / (Type(1) + exp(-(alpha(k) + beta(k) * fit_len(i))));
      fit_slx(i,k) = Type(1) / (Type(1) + exp(Type(-1.0) * slp(k) * (len(i) - s50(k)))); 
    }
  }

  // Fitted phi (probability of capture)
  matrix<Type> fit_phi(nfit,ntrt);
  fit_phi.setZero();

  for (int i = 0; i < nfit; i++) {
    for (int k = 0; k < ntrt; k++) {
    fit_phi(i,k) = fit_slx(i,k) / (fit_slx(i,k) + delta);
    }
  }

  // Length at 25%, 50%, and 75%
  // vector<Type> s25(ntrt);
  // vector<Type> s50(ntrt);
  // vector<Type> s75(ntrt);
  // 
  // for (int k = 0; k < ntrt; k++) {
  //   s25(k) = - (alpha(k) + log(Type(3))) / (beta(k));
  //   s50(k) = - alpha(k) / beta(k);
  //   s75(k) = - (alpha(k) - log(Type(3))) / (beta(k));
  // }

  // Selection range (slx_rng, difference between s75 and s25); selection factor (slx_fct,
  // s50 weighted by escape ring size, trt); slope of the logistic curve
  // (slx_slp)
  // vector<Type> slx_rng(ntrt);
  // vector<Type> slx_fct(ntrt);
  // vector<Type> slx_slp(ntrt);
  // 
  // for (int k = 0; k < ntrt; k++) {
  //   slx_rng(k) = s75(k) - s50(k);
  //   slx_fct(k) = s50(k) * trt(k);
  //   slx_slp(k) = beta(k) / Type(4);
  // }

  // Linear regression between length at 50% selectivity and escape ring diameter
  for (int k = 0; k < ntrt; k++) {
    s50(k) = a1 + b1 * trt(k);
  }

  // Linear regression between slope of logistic selectivity curve and escape
  // ring diameter
  for (int k = 0; k < ntrt; k++) {
    slp(k) = a2 + b2 * trt(k);
  }

  // OBJECTIVE FUNCTION -----

  // // Negative log likelihood
  Type nll = 0;
  
  // Type<matrix> sum_trt

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

  // What might the penalized likelihood be for a linear regression? One option
  // is to use the parameter estimates as the observed and the estimates from
  // the theoretical selectivity curves as the predicted, and use sum of squares
  // as the nll

  Type penl_s50 = 0;
  Type penl_slp = 0;

  for (int k = 0; k < ntrt; k++) {
    penl_s50 += square(theor_s50(k) - s50(k));
    penl_slp += square(theor_slp(k) - slp(k));
  }

  penl_s50 *= wt_s50;  // weight for the s50 penalty
  penl_slp *= wt_slp;  // weight for the slp penalty

  nll += penl_s50;
  nll += penl_slp;

  // Priors to constrain selectivity curve at 0 and 100
  Type s0 = 0;
  Type prior_s0 = 0;
  Type s100 = 0;
  Type prior_s100 = 0;
  
  for (int j = 0; j < nset; j++) {
    for (int k = 0; k < ntrt; k++) {
      s0 = fit_slx(0,k);
      s100 = fit_slx(nlen-1,k);
      prior_s0 += Type(0.5) * square(s0 - Type(0)) / square(sigma_s0);
      prior_s100 += Type(0.5) * square(s100 - Type(1)) / square(sigma_s100);
    }
  }
  nll += prior_s0;
  nll += prior_s100;
  
  // Troubleshoot model
  // nll = dummy * dummy;
  
  // REPORT SECTION -----
  
  REPORT(fit_slx);
  // REPORT(alpha);
  // REPORT(beta);
  REPORT(s50);
  REPORT(slp);
  REPORT(slx);
  REPORT(phi);
  REPORT(fit_phi);
  REPORT(r);
  // REPORT(s50);
  // REPORT(s25);
  // REPORT(s75);
  // REPORT(slx_rng);
  // REPORT(slx_fct);
  // REPORT(slx_slp);
  REPORT(set_effect);
  REPORT(penl_s50);
  REPORT(penl_slp);
  REPORT(nll);
  REPORT(prior_s0);
  REPORT(prior_s100);
  
  return(nll);
  
}
