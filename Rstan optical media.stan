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
  vector[N*K] t; //time conditions
  vector[N*K] x1; //log(RH)
  vector[N*K] x2;
  vector[N*K] y; //the response variable
  // vector[N*K] RH;
  // vector[N*K] T;

  int<lower = 1> wishart_df; //degree of freedom in the wishart dist.
}

// transformed data{
//   vector[N*K] x1; //log(RH)
//   vector[N*K] x2;
//   x2 = 11605/(T+273.15);
//   x1 = log(RH);
//   
// }
// The parameters accepted by the model. 
parameters {
  vector<lower=0>[K] tau; //the standard deviation of the regression coefficients
  real<lower=0> logA;
  real<lower=0> B;
  real<lower=0> delta_H;
  matrix[(N-1)*5+K,(N-1)*5+K] mat;
  vector[N*K] mu_mat;
  corr_matrix[N*K] sigma_mat;
  corr_matrix[N*K] sigma_for_prior;
}


transformed parameters{
    vector[(N-1)*5+K] beta1;
    vector[(N-1)*5+K] mu;
    for (n in 1:N){
      for(k in 1:K){
     beta1[(n-1)*5+k] = exp(logA + B*x1[(n-1)*5+k] + delta_H*x2[(n-1)*5+k]);
     mu[(n-1)*5+k] = mat[(n-1)*5+k,1] + beta1[(n-1)*5+k]*(pow(t[(n-1)*5+k],mat[(n-1)*5+k,2]));
    }}
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
  logA ~ normal(0,1e-3);
  B ~ normal(0, 1e-3);
  delta_H ~ normal(0,1e-3);
  mat[,1] ~ multi_normal_cholesky_lpdf(mu_mat,sigma_mat);

  // sigma_mat ~ wishart(wishart_df,sigma_for_prior);
  // sigma_for_prior ~ lkj_corr_cholesky(1.5);

  
  for(n in 1:N){
    for(k in 1:K){
     y[(n-1)*5+k] ~ lognormal(mu[(n-1)*5+k], tau );  
    }}
}







