#function to run EP on the astronomy data
EP = function(fit = NULL, data = NULL, J = 360, K = 6, prior_Mu = log(c(200, 250, 5, .5, 50, 7, 1, 50, .5)), 
              prior_Sigma_inv = diag(9) * 6/10, S = 10, SMOOTH = 0.9, randomSites = FALSE, parallel = TRUE, mc_iter = 100){

  # Extract data
  x <- data$x;
  y <- data$y;
  bin <- data$bin;
  P <- length(prior_Mu);
  prior_Sigma_inv_Mu <- prior_Sigma_inv %*% prior_Mu;
  
  # Suffle Groups Among Sites If Needed
  if (randomSites) {
    bin_order <- bin_samp <- sample(unique(bin), replace = FALSE);
  } else {
    bin_order <- bin_samp <- 1:max(bin);
  }
  
  # Initialize global and local natural parameters
  Sigma_k_inv_Mu <- Sigma_k_inv <- Sigma_inv <- Sigma_inv_Mu <- Sigma_k_inv_tilt <- Post_Sigma <- list();
  Post_Mu <- matrix(0, ncol = P, nrow = S + 1);
  eta_j <- a_j <- tilt_fits <- list();
  for(s in 1:S) { eta_j[[s]] <- a_j[[s]] <- matrix(0, nrow = P/2, ncol = length(bin_samp)); } 
  
  for(k in 1:K) Sigma_k_inv_Mu[[k]] <- rep(0, P);
  for(k in 1:K) Sigma_k_inv[[k]] <- diag(P) * 0;
  Sigma_inv[[1]] <- prior_Sigma_inv;
  Sigma_inv_Mu[[1]] <- prior_Sigma_inv_Mu;
  
  Post_Sigma[[1]] <- solve(Sigma_inv[[1]]);
  Post_Mu[1,] <- prior_Mu;
  init_data <- list();

  # timers
  init_time <- proc.time();
  k_times <- matrix(0, nrow = K, ncol = S);

  # The actual algorithm...
  for(s in 1:S){
    #update natural parameters from previous iteration
    if (parallel) {
      Sigma_inv[[s+1]] <- prior_Sigma_inv;
      Sigma_inv_Mu[[s+1]] <- prior_Sigma_inv_Mu;
    } else if (s > 1) {
      Sigma_inv[[s]] <- Sigma_inv[[s-1]];
      Sigma_inv_Mu[[s]] <- Sigma_inv_Mu[[s-1]];
    }
    
    for(k in 1:K){
      k_times[k,s] <- proc.time()[3];
      cat(" Current Iteration Status: [", s, " out of ", S, "] \n");
      cat(" Current Partition Status: [", k, " out of ", K, "] \n");
      
      # 1. Update the Cavity Distribution...
      Sigma_inv_cav <- Sigma_inv[[s]] - Sigma_k_inv[[k]];
      Sigma_inv_Mu_cav <- Sigma_inv_Mu[[s]] - Sigma_k_inv_Mu[[k]];
      Sigma_cav <- solve(Sigma_inv_cav);
      Sigma_cav <- (Sigma_cav + t(Sigma_cav))/2; #preserve symmetry
      Mu_cav <- as.vector(Sigma_cav %*% Sigma_inv_Mu_cav);
    
      # 2. Find tilted distribution in Stan...
      # Extract current partition of the data...
      if (randomSites) { 
        bin_cur <- bin_samp[((k-1)*J/K+1) : (k*J/K)];
        subset <- which(bin %in% bin_cur);
        bin_k <- bin[subset];
        
        bin_order_k <- order(bin_cur);
        for(b in 1:length(bin_k)) {
          bin_k[b] <- bin_order_k[which(bin_cur == bin_k[b])];
        }
      } else { 
        subset <- which(bin <= k*J/K & bin > (k-1)*J/K);
        bin_k <- (ceiling(bin[subset] - 1))%%(J/K) + 1; 
        bin_order_k <- 1:(J/K);
      }
      y_k <- y[subset];
      x_k <- x[subset];
      B <- length(unique(bin_k));
      tilt_data <- list(N = length(y_k), M = P, B = B, x = x_k, y = y_k, bin = bin_k, Mu_Cav = Mu_cav, Sig_Cav = Sigma_cav);
      
      # Fit tilted distribution in Stan....
      for(i in 1:4) { init_data[[i]] <- list(eta = matrix(0, nrow = P/2, ncol = J/K), phi = Mu_cav);}
      tilt_fit <- sampling(fit, data = tilt_data, iter = mc_iter, chains = 4, init = init_data);
      tilt_fits[[k]] <- tilt_fit;
      
      # Extract local parameter means
      current_cols = ((J/K)*(k-1) + 1):(J/K*k)
      eta_k <- apply(extract(tilt_fit)$eta, c(2,3), mean);
      eta_j[[s]][, current_cols] <- eta_k[, bin_order_k];
      a_k <- apply(extract(tilt_fit)$a, c(2,3), mean);
      a_j[[s]][, current_cols] <- a_k[, bin_order_k];
      bin_order[current_cols] <- bin_order_k;
      
      # Extract global parameter mean and covariance matrix....
      Mu_tilt <- colMeans(extract(tilt_fit)$phi);
      Sigma_tilt <- matrix(0, nrow = P, ncol = P);
      for(i in 1:P)
        for(j in 1:P)
          Sigma_tilt[i,j] <- cov(extract(tilt_fit)$phi[,i], extract(tilt_fit)$phi[,j]);
      
      print(diag(Sigma_cav) - diag(Sigma_tilt));
      
      # Bias correction on covariance matrix....
      n <- length(extract(tilt_fit)$phi[,1]);
      Sigma_inv_tilt <- solve(Sigma_tilt) * (n-P-2)/(n-1);
      Sigma_inv_Mu_tilt <- Sigma_inv_tilt %*% Mu_tilt;
      Sigma_k_inv_tilt[[k]] <- Sigma_inv_tilt;
    
      delta = 1;
      repeat {
        # 3. Update the Site Distribution
        Sigma_k_inv[[k]] <- (Sigma_inv_tilt - Sigma_inv_cav) * delta; 
        Sigma_k_inv_Mu[[k]] <- (Sigma_inv_Mu_tilt - Sigma_inv_Mu_cav) * delta;
        
        # 4. Update g(phi) in Parallel/Serial Until Singular
        if (parallel) {
          Sigma_inv_proposal <- Sigma_inv[[s+1]] + Sigma_inv_cav * (1 - delta) + Sigma_k_inv[[k]];
          Sigma_inv_Mu_proposal <- Sigma_inv_Mu[[s+1]] + Sigma_inv_Mu_cav * (1 - delta) + Sigma_k_inv_Mu[[k]];
          
          n_negative <- length(which(diag(solve(Sigma_inv_proposal)) < 0));
          
          if((n_negative > 0) && (delta > 0.000001)){
            cat("Failure (Invalid Covariance): ...", delta, "... ");
            delta = delta * SMOOTH;
          } else if((n_negative > 0) && (delta < 0.000001)) {
            cat("\n\nSite was discarded!");
            break;
          } else {
            cat("\nSite was succesfully added!\n\n");
            Sigma_inv[[s+1]] <- Sigma_inv_proposal;
            Sigma_inv_Mu[[s+1]] <- Sigma_inv_Mu_proposal;
            break;
          }
        } else {
          Sigma_inv[[s]] <- Sigma_k_inv[[k]] + (1-delta) * Sigma_inv_cav;
          Sigma_inv_Mu[[s]] <- Sigma_k_inv_Mu[[k]] + (1-delta) * Sigma_inv_Mu_cav;
        }
      }
      
      k_times[k,s] <- proc.time()[3] - k_times[k,s];
    }
    
    eta_j[[s]] <- t(eta_j[[s]]);
    a_j[[s]] <- t(a_j[[s]]);
    
    # Convert Natural Parameters Back Into Usual Framework
    if (parallel) {
      Post_Sigma[[s+1]] = solve(Sigma_inv[[s+1]]);
      Post_Mu[s+1,] = as.vector(Post_Sigma[[s+1]] %*% Sigma_inv_Mu[[s+1]]);
    } else if (s > 1) {
      Post_Sigma[[s+1]] = solve(Sigma_inv[[s]]);
      Post_Mu[s+1,] = as.vector(Post_Sigma[[s+1]] %*% Sigma_inv_Mu[[s]]);
    }
  }

  #save the total time elapsed
  final_time <- proc.time();
  
  return(list(Post_Sigma = Post_Sigma, Post_Mu = Post_Mu, Sigma_k_inv = Sigma_k_inv, Sigma_k_inv_tilt = Sigma_k_inv_tilt, eta_j = eta_j, a_j = a_j, tilt_fits = tilt_fits,
              time = final_time - init_time, k_times = k_times, bin_samp = bin_samp, bin_order = bin_order));
}