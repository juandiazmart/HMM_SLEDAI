
functions {
  void indlik_lp(int K, vector y,matrix X,vector alpha, vector[] theta, vector mu) {
    real acc[K];
    real gamma[rows(y), K];
      for (k in 1:K) 
      gamma[1,k] = poisson_lpmf(y[1] | log(mu[k]+ X[1]*alpha));
  for (t in 2:rows(y)) {
    for (k in 1:K) {
      for (j in 1:K)
        acc[j] = gamma[t-1, j] + log(theta[j, k]) + poisson_lpmf(y[t] | log(mu[k]+X[t]*alpha));
      gamma[t, k] = log_sum_exp(acc);
    }
  }
  target += log_sum_exp(gamma[rows(y)]);
  }
}


data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> I;
  int<lower=0> Npreds;
  vector[N] y;
  matrix[N,Npreds] X;
  int index[I+1];
}

parameters {
  simplex[K] theta[K];
  ordered[K] mu;
  vector[Npreds] alpha;
  
}

model {
  // priors
  target+= normal_lpdf(mu[1] | 0, 1);
  target+= normal_lpdf(mu[2] | 0, 1);
  // forward algorithm
  for (i in 1:I){
    indlik_lp(K, y[index[i]:(index[i+1]-1)],X[index[i]:(index[i+1]-1)],alpha, theta,mu);
  }
}

generated quantities {
  int<lower=1,upper=K> z_star[N];
  real log_p_z_star;
  vector[N] log_lik;
  vector[N] yrep;
  {
    int back_ptr[N, K];
    real best_logp[N, K];
    for (k in 1:K)
      best_logp[1, k] = poisson_lpmf(y[1] | log(mu[k]+ X[1]*alpha));
    for (t in 2:N) {
      for (k in 1:K) {
        best_logp[t, k] = negative_infinity();
        for (j in 1:K) {
          real logp;
          logp = best_logp[t-1, j] + log(theta[j, k]) + poisson_lpmf(y[1] | log(mu[k]+ X[t]*alpha));
          if (logp > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp;
          }
        }
      }
    }
    log_p_z_star = max(best_logp[N]);
    for (k in 1:K)
      if (best_logp[N, k] == log_p_z_star)
      z_star[N] = k;
    for (t in 1:(N - 1))
      z_star[N - t] = back_ptr[N - t + 1, z_star[N - t + 1]];
  }
}
