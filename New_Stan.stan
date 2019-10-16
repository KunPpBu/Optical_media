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
  int<lower=0> N; //number of units
  // Number of level-2 clusters
  int<lower=0> K; //number of conditions
  // RH from the data
  real RH[N*K];
  // Temp from the data
  int<lower=0> T[N*K];
  // Continuous outcome
  real y_ijk[N*K];//number of the outputs
  // Continuous predictor
  real t_ijk[N*K];
}

//Prepocessing of the data
transformed data{
  real x1[N*K];
  real x2[N*K];
  x1 = log(RH);
  x2[N*K] = 11605/(T[N*K]+273.15);
}

// The parameters accepted by the model.
parameters {
  // Define parameters to estimate
  // Population intercept (a real number)
  real beta0[N];
  // Level-1
  real<lower=0> sigma;
  // Random effect
  matrix[N,2] mat;
  matrix[2,2] sigma_mat;
  // Hyperparameters
  vector[2] mu_mat;
  real<lower=0> A;
  real<lower=0> B;
  real<lower=0> delta_H;
  corr_matrix[2] sigma_for_prior;
}

//Parameters processing before the postier is computed
transformed parameters{
  // Population slope
  real beta1[(N-1)*5+K];
  for(k in 1:K){
    for(n in 1:N){
       beta1[(n-1)*5+k] = exp(log(A) + B*x1[(n-1)*5+k] + delta_H*x2[(n-1)*5+k]); 
    }
  }
   
   
}


model {
  sigma ~ gamma(1e-3,1e-3);
  mu_mat[1] ~ normal(1e-3,1e-3);
  mu_mat[2] ~ normal(1e-3,1e-3);
  A ~ normal(0,1e-3);
  B ~ normal(0, 1e-3);
  delta_H ~ normal(0,1e-3);
  mat[2] ~ multi_normal(mu_mat,sigma_mat);
  sigma_mat ~ wishart(3,sigma_for_prior);
  sigma_for_prior ~ lkj_corr_cholesky(1.5);
  for(n in 1:N){
    for(k in 1:K){
       y_ijk[(n-1)*5+k] ~ lognormal(beta0[n] + beta1[(n-1)*5+k]*pow(t_ijk[N*K],mat[N,2]), sigma);
    }
  }
 
}





