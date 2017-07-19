require("rstan")
options(mc.cores = parallel::detectCores())

#  ..... Compile Stan Code .....
model_hier_sigmoid_dip <- stan_model(model_code = code_hier_sigmoid_dip, model_name = "Hierarchical Sigmoid with Dip");

# ..... Priors .....
J = 50;
pri_mu <- c(rep(0, 5), log(3.5), 0, 0, rep(-1, 8));
pri_sig <- diag(rep(4, 16));

init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu, eta = matrix(0, nrow = 8, ncol = J));}

# ...... Simulate the Data ......
N = 500;
actual_mu <- c(log(5.5), log(1.5), log(0.85), log(3.5), -0.471, log(3.37), log(0.46), log(0.14), rep(-1, 8));
actual_sig <- diag(rep(0.1, 16));

list_hier_sigmoid_dip <- hier_sim(actual_mu, actual_sig, J, N, sigmoid_dip, log(i100[i100>0]), long[i100>0]);
data_hier_sigmoid_dip <- list_hier_sigmoid_dip$data;

plot(data_hier_sigmoid_dip$x, data_hier_sigmoid_dip$y, xlab = "x", ylab = "y", main = "Hierarchical Sigmoid With Dip Model", bty = "l", col = data_hier_sigmoid_dip$bin);
plot(exp(data_hier_sigmoid_dip$x[data_hier_sigmoid_dip$bin==8]), exp(data_hier_sigmoid_dip$y[data_hier_sigmoid_dip$bin==8]), xlab = "x", ylab= "y", bty = 'l', col = 2);

# ...... EP Fit .......
fit_hier_sigmoid_EP_5 <- EP(fit = model_hier_sigmoid_dip, J = J, K = 5, data = data_hier_sigmoid_dip, S = 1, mc_iter = 1000, prior_Mu = pri_mu, prior_Sigma_inv = solve(pri_sig), randomSites = TRUE);
fit_hier_sigmoid_EP_10 <- EP(fit = model_hier_sigmoid_dip, J = J, K = 10, data = data_hier_sigmoid_dip, S = 1, mc_iter = 1000, prior_Mu = pri_mu, prior_Sigma_inv = solve(pri_sig), randomSites = TRUE);
fit_hier_sigmoid_EP_25 <- EP(fit = model_hier_sigmoid_dip, J = J, K = 25, data = data_hier_sigmoid_dip, S = 1, mc_iter = 1000, prior_Mu = pri_mu, prior_Sigma_inv = solve(pri_sig), randomSites = TRUE);

# ...... MC Fit ....... 
data <- list(N = N * J, B = J, x = data_hier_sigmoid_dip$x, y = data_hier_sigmoid_dip$y, bin = data_hier_sigmoid_dip$bin, Mu_Cav = pri_mu, Sig_Cav = pri_sig);
MC_time <- system.time(fit_hier_sigmoid_MC <- sampling(model_hier_sigmoid_dip, data = data, iter = 1000, chains = 4, init = init_data));
fit_hier_sigmoid_MC@par_dims$time <- MC_time[3];

# ...... Plot the Fits ......
plot_fits(data = data_hier_sigmoid_dip, fits = list(HMC = fit_hier_sigmoid_dip_MC, fit_hier_sigmoid_dip_EP_5, fit_hier_sigmoid_dip_EP_10, fit_hier_sigmoid_dip_EP_25), 
          names = c("HMC", "EP 5", "EP 10", "EP 25"), cols = c("black", "red", "blue", "darkgreen"), folder = "SigmoidDipFits", inc = 1, 
          funcs = list(sigmoid_dip), should_exp = T, phi_true = list_hier_sigmoid_dip$phi);
