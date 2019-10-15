//

// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data. 
data {
  // Define variables in data
  // Number of level-1 observations (an integer)
  int<lower=0> Ni; //number of units
  // Number of level-2 clusters
  int<lower=0> Nj; //number of conditions
  // Cluster IDs
  int<lower=1> leverl2[Ni];
  // RH from the data
  real RH;
  // Temp from the data
  int<lower=0> T;
  // Continuous outcome
  real y_ijk[Ni];//number of the outputs
  // Continuous predictor
  real t_ijk;
}

//Prepocessing of the data
transformed data{
  real x1;
  real x2;
  x1 = log(RH);
  x2 = 11605/(T+273.15);
}

// The parameters accepted by the model.
parameters {
  // Define parameters to estimate
  // Population intercept (a real number)
  real beta0;
  // Level-1
  real<lower=0> sigma;
  // Random effect
  matrix[Ni,2] mat;
  vector[2] mu_mat;
  matrix[2,2] sigma_mat;
  // Hyperparameters
  real<lower=0> A;
  real<lower=0> B;
  real<lower=0> delta_H;
}

//Parameters processing before the postier is computed
transformed parameters{
  // Population slope
  real beta1;
   beta1 = exp(log(A) + B*x1 + delta_H*x2); 
}


model {
  y_ijk ~ lognormal(beta0 + beta1*pow(t_ijk,mat[Nj,2]), sigma);
}


