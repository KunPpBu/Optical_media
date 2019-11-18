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
  int<lower=1> N; //number of units
  // Number of level-2 clusters
  int<lower=1> K; //number of conditions
  // RH from the data
  real RH[(N-1)*5+K];
  // Temp from the data
  int<lower=1> T[(N-1)*5+K];
  // Continuous outcome
  vector[(N-1)*5+K] y_ijk;//number of the outputs
  // Continuous predictor
  real t_ijk[(N-1)*5+K];
  //zero mean for multi_normal 
  vector[2] mean_mat;
}

//Prepocessing of the data
transformed data{
  real x1[(N-1)*5+K];
  real x2[(N-1)*5+K];
  for(k in 1:K){
    for(n in 1:N){
      x1[(n-1)*5+k] = log(RH[(n-1)*5+k]);
      x2[(n-1)*5+k] = 11605/(T[(n-1)*5+k]+273.15);
      // print("x1:", x1);
    }
  }
  
}

// The parameters accepted by the model.
parameters {
  // Define parameters to estimate
  // Random effect
  // matrix[N,2] mat;
  // Level-1
  real<lower=0> sigma;
  // Hyperparameters
  vector[2] mu_mat;
  real A;
  real B;
  real delta_H;
  corr_matrix[2] sigma_mat;
  corr_matrix[2] sigma_for_prior;
}

//Parameters processing before the postier is computed
transformed parameters{
  // Random effect
  row_vector[N] beta0;
  row_vector[N] gamma;
  row_vector[(N-1)*5+K] beta1;
  row_vector[(N-1)*5+K] mu;
  // Population slope
  for(n in 1:N){
      beta0[n] = mu_mat[1];
      gamma[n] = mu_mat[2];
      // print("beta0:", beta0);
  }
  for(n in 1:N){
    for(k in 1:K){
      beta1[(n-1)*5+k] = exp(log(A) + B*x1[(n-1)*5+k] + delta_H*x2[(n-1)*5+k]); 
      mu[(n-1)*5+k] = beta0[n] + beta1[(n-1)*5+k] * pow(t_ijk[(n-1)*5+k],gamma[n]);
      // print("mu:", mu);
    }
  }
}


model {
  mu[(N-1)*5+K] ~ normal(0, sqrt(1000));
  sigma ~ gamma(1e-3,1e-3);
  mu_mat[1] ~ normal(0, sqrt(1000));
  mu_mat[2] ~ normal(0, sqrt(1000));
  A ~ normal(0, sqrt(1000));
  B ~ normal(0, sqrt(1000));
  delta_H ~ normal(0, sqrt(1000));
  mu_mat ~ multi_normal(mean_mat,sigma_mat);
  sigma_mat ~ wishart(2,sigma_for_prior);
  sigma_for_prior ~ lkj_corr(1);
  for (n in 1:N){
    for (k in 1:K){
      // target += lognormal_lpdf(y_ijk[(n-1)*5+k] | mu[(n-1)*5+k],sigma);
      y_ijk[(n-1)*5+k] ~ lognormal(mu[(n-1)*5+k], sigma);
    }
  }
  
}

