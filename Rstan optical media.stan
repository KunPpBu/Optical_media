//
//  edited on Friday
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

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; //the number of observations for each group
  int<lower=0> K; //the number of columns in the model matrix
  vector[N*K] t; //time conditions
  vector[N*K] x1; //log(RH)
  vector[N*K] x2; //(11605)/(T+273.15)
  vector[N*K] y; //the response variable
  // matrix[N,N] R;

}


// The parameters accepted by the model.
parameters {
  vector<lower=0>[K] tau; //the standard deviation of the regression coefficients
  cov_matrix[N] sigma_mat;
  // real<lower=0> mu_mat;
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] mu;
  vector[N] mu_mat;
  row_vector[N] mat;
  vector[N] sigma;
  vector[N] logA;
  vector[N] B;
  vector[N] delta_H;
  vector[N] beta1;
  matrix[N,N] R;
  matrix[2,2] sigma_1;


  //priors
  // sigma_mat ~ wishart(90,R);
  // logA ~ normal(0,1e-3);
  // B ~ normal(0, 1e-3);
  // delta_H ~ normal(0,1e-3);
  tau ~ gamma(1e-3,1e-3);
  // mu_mat[N,1] ~ normal(1e-3,1e-3);
  // mu_mat[N,2] ~ normal(1e-3,1e-3);


  for(n in 1:N){
   mat ~ multi_normal(mu_mat,sigma_mat);
   sigma_mat ~ wishart(90,R);
   mu_mat[1] ~ normal(1e-3,1e-3);
   mu_mat[2] ~ normal(1e-3,1e-3);
   beta1 = exp(logA + B*x1[n] + delta_H*x2[n]);
   mu[n] = mat[1] + (beta1[n]*pow(t[n],mat[2]));
   }
    y ~ lognormal(mu, tau );
}

// generated quantities {
// real<lower=0,upper=J> beta_j;
// beta_j = exp(logA+B*x1+delta_H*x2);
// }
