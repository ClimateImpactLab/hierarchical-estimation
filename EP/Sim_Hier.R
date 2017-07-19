# Hierarchical Model With a 2P-Dimensional Phi
# First P parameters correspond to hierarchical mean
# Second P parameters correspond to hierarchical standard deviation
# Pth parameter must correspond to local standard deviation
# func should be a function of x and the first P-1 parameters
sim_hier <- function(pri_mu = log(c(200, 2.5, 2.5)), pri_sig = c(0.5, 0.5, 0.5), J = 360, n_samp = 20, 
                     func = linear, raw_x = 0, raw_long = 0, period = 0){
  x <- y <- bin <- c();
  P <- length(pri_mu)/2;
  
  phi <- pri_mu + rnorm(P * 2, 0, 1) * sqrt(diag(pri_sig));
  mu_a <- phi[1:P];
  tau_a <- exp(phi[(P+1):(P * 2)]);
  a <- matrix(0, nrow = P, ncol = J);
  
  print(phi)
  
  for(j in 1:J){
    if(period == 0) {
      a[,j] <- mu_a[1:P] + rnorm(P, 0, 1) * tau_a[1:P];
    } else {
      a[,j] <- mu_a[1:P] + sin(period * j) * tau_a[1:P];
    }
    
    if(length(raw_x) < n_samp){
      x_j <- runif(n_samp, 0, 50);
    }else{
      x_j <- sample(raw_x[(raw_long < j) & (raw_long >= j-1)], n_samp, rep = FALSE);
    }
    y_j <- func(x_j, a[,j]) + rnorm(n_samp, 0, 1) * exp(a[P,j]);
    
    x <- c(x, x_j);
    y <- c(y, y_j);
    bin <- c(bin, rep(j, n_samp));
  }
  
  return(list(data = data.frame(x = x, y = y, bin = bin), phi = phi, a = a));
}

# Identity
identity <- function(x, a) {
  return(a[1]);
}

# Simple Linear Function
linear <- function(x, a) {
  return(exp(a[1]) + exp(a[2]) * x);
}

# Inverse logistic function
inv_logit <- function(x){ 
  return(1/(1 + exp(-x))) 
}

# Generalized inverse logistic
sigmoid <- function(x, a) {
  result <- exp(a[1]) + inv_logit((x - exp(a[2]))/exp(a[3])) * exp(a[4]);
  return(result)
}

# Generalized inverse logistic with Gaussian dip
sigmoid_dip <- function(x, a) {
  result <- 1 - inv_logit(a[5]) * exp(-0.5 * ((x - exp(a[6]))/exp(a[7]))^2);
  result <- result * inv_logit((x - exp(a[2]))/exp(a[3])) * exp(a[4]);
  result <- exp(a[1]) + result;
  return(result)
}

# Mixture Between Sigmoid and Simgoid with Dip
mixture_curve <- function(x, a) {
  u <- runif(length(x), 0, 1);
  cutoff <- ifelse(x > log(5), ifelse(x < log(17), inv_logit(a[8]), inv_logit(a[9])), 0);
  #cutoff <- ifelse(x<log(17), inv_logit(a[8]), inv_logit(a[9]));
  ind <- as.numeric(u < cutoff);
  return(ind * sigmoid_dip(x,a) + (1-ind) * sigmoid(x,a));
}