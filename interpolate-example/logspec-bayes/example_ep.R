#MODEL
code_ep = '
  data{
    int<lower=1> N;
    int<lower=1> B;
    real x[N];
    real y[N];
    int<lower=1> bin[N];
    vector[16] Mu_Cav;
    matrix[16,16] Sig_Cav;
  }

  transformed data{
    matrix[16,16] L;
    L <- cholesky_decompose(Sig_Cav);
  }

  parameters{
    vector[16] phi;
    matrix[8,B] eta;
  }

  transformed parameters{
    matrix[8,B] a;

    for(p in 1:8){
      for(b in 1:B){
        a[p,b] <- phi[p] + exp(phi[p + 8]) * eta[p,b];
      }
    }
  }

  model{
    vector[N] y_hat;
    vector[N] sigma_hat;

    phi ~ multi_normal_cholesky(Mu_Cav, L);

    for(p in 1:8){
      for(b in 1:B){
        eta[p,b] ~ normal(0,1);
      }
    }

    for(n in 1:N) {
      y_hat[n] <- 1 - inv_logit(a[5,bin[n]]) * exp( -0.5 * pow( ( x[n] - exp(a[6,bin[n]]) ) / exp(a[7,bin[n]]), 2 ) );
      y_hat[n] <- y_hat[n] * inv_logit( ( x[n] - exp(a[2,bin[n]]) ) / exp(a[3,bin[n]]) ) * exp(a[4,bin[n]]);
      y_hat[n] <- exp(a[1,bin[n]]) + y_hat[n];
      sigma_hat[n] <- exp(a[8, bin[n]]);
    }

    y ~ normal(y_hat, sigma_hat);
  }
'
