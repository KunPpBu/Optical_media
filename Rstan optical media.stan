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

// A short intro of STAN
//https://faculty.ai/blog/a-short-introduction-to-stan/


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; //the number of observations in each group
  int<lower=0> K; //the number of columns in the model matrix
  real<lower=0,upper=N*K> t; //time conditions
  real<lower=0,upper=N*K> RH;
  real<lower=0,upper=N*K> T;
  row_vector[N*K] y; //the response variable
  int<lower = 1> wishart_df; //degree of freedom in the wishart dist.
}

 transformed data{
  real<lower=0,upper=(N-1)*5+K> x1; //log(RH)
  real<lower=0,upper=(N-1)*5+K> x2;
  x2 = 11605/(T+273.15);
  x1 = log(RH);
}


// The parameters accepted by the model. 
parameters {
  row_vector[N*K] tau; //the standard deviation of the regression coefficients
  real<lower=0,upper=90> A;
  real<lower=0,upper=90> B;
  real<lower=0,upper=90> delta_H;
  matrix[N*K,2] mat;
  matrix[N*K,2] mu_mat;
  corr_matrix[2] sigma_mat;
  corr_matrix[2] sigma_for_prior;
}


transformed parameters{
    real<lower=0,upper=(N-1)*5+K> beta1;
    row_vector[N*K] mu;
    beta1 = exp(log(A) + B*x1 + delta_H*x2);
    mu[N*K] = mat[N*K,1] + beta1*(pow(t,mat[N*K,2]));
}

// The model to be estimated. We model the output
// 'y' to be log normally distributed with mean 'mu'
// and standard deviation 'sigma', (tau equals to inverse sigma)
model {
  vector[N*K] sigma; 
  matrix[N,N] R;
  matrix[2,2] sigma_1;


  //priors
  tau ~ gamma(1e-3,1e-3);
  mu_mat[1] ~ normal(1e-3,1e-3);
  mu_mat[2] ~ normal(1e-3,1e-3);
  A ~ normal(0,1e-3);
  B ~ normal(0, 1e-3);
  delta_H ~ normal(0,1e-3);
  mat[N*K,] ~ multi_normal(mu_mat[N*K,],sigma_mat);
  sigma_mat ~ wishart(wishart_df,sigma_for_prior);
  sigma_for_prior ~ lkj_corr_cholesky(1.5);

  for(n in 1:N){
    for(k in 1:K){
     y[(n-1)*5+k] ~ lognormal(mu[(n-1)*5+k], tau );  
    }}
}







