### Generate some data for testing
# Y_{N×M} : indicator matrix such that y_{n,m} denote the existence of mutation m on n-th 
#           individual genome sequence (0=does not exist; 1=exist);
# X_{N×(D+1)}: the design matrix. For n-th patient, we have x_n = [1, x_{n,1}, x_{n,2}, ..., x_{n,D}]⊤
# β_k : a (D + 1) × 1 vector that stores the logistic regression intercept + coefficients in cluster k
#       for mutations in cluster k, i.e. βk = [βk,0, βk,1, βk,2, ..., βk,D]⊤
# π_{n,k}: indicates the probability of existence of mutation m on the n-th genome sequence given
#          that mutation m belongs to cluster k
# Z_{m,k}: whether or not (1 or 0) mutation m belongs to cluster k
# λ_{m,k} = P (Z_{m,k} = 1) :the probability that the mutation m belongs to cluster k.

rm(list=ls())

# ----------------Set the working directory to current path -------------------#
library("rstudioapi")
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

############################# Data Simulation ##################################
set.seed(888)
N <- 5000   # number of observations
M <- 20     # number of mutations
K <- 3      # number of clusters
D <- 10     # number of covariates
D_cont <- 5 # number of continuous covariates
D_cat <- D - D_cont  # number of categorical covariates


# Simulate predictors: 
# I. All continuous
# X <- cbind(rep(1, N), matrix(rnorm(N*D), nrow=N, ncol=D))

# II. Some continuous, some categorical
# Generate the continuous predictors
X_cont <- matrix(rnorm(N*D_cont), nrow=N, ncol=D_cont)

# Generate the categorical predictors
X_cat <- matrix(sample(0:1, N*D_cat, replace=TRUE, prob=c(0.9, 0.1)), nrow=N, ncol=D_cat)

# Bind the predictors together
X <- cbind(rep(1, N), X_cont, X_cat)

# Simulate coefficients
beta <- matrix(rnorm(K*(D+1)), nrow=K, ncol=D+1)

# # Define the logit link function
# logit <- function(x) {
#   exp(x) / (1 + exp(x))
# }

# Simulate cluster membership and observations
lambda <- matrix(NA, nrow=M, ncol=K)
Z <- matrix(NA, nrow=M, ncol=K)
pi <- matrix(NA, nrow=N, ncol=K)
Y <- matrix(NA, nrow=N, ncol=M)
for (m in 1:M) {
  
  # Sample lambda for each mutation m
  lambda[m,] <- runif(K)
  lambda[m,] <- lambda[m,] / sum(lambda[m,]) # normalize so that it sums to 1
  
  # Sample Z for mutation m
  Z[m,] <- rmultinom(n = 1, size = 1, prob = lambda[m,])
  k <- which(Z[m,] == 1)
  
  # delete loop
  pi[,k] <- plogis(beta[k,] %*% t(X))
  Y[,m] <- rbinom(N, size = 1, prob = pi[,k])
 
  # for (n in 1:N){
  #   pi[n,k] = logit(beta[k,] %*% X[n,])
  #   Y[n,m] <- rbinom(1, size = 1, prob = pi[n,k])
  # }
  
}

### Print true simulated data

# Create text file
file_name <- paste0("../Results/Simulation Results/TrueSimulation_K", K, "_N", N, "_M", M, "_D", D, ".txt")
file_conn <- file(file_name, open = "w")
clusterIndex <- rep(0, M)
MutationName <- rep("", M)
for(m in 1:M){
  MutationName[m] <- paste("mutation ", m)
  clusterIndex[m] <- which(Z[m,]==1)
}

for(k in 1:K){
  writeLines(paste("Cluster ", k), file_conn)
  writeLines(paste(MutationName[which(clusterIndex==k)]), file_conn)
  # TrueOR <- exp(TrueBeta[k,]) # print odds ratio
  # print(TrueOR)
  TrueBeta.k <- beta[k,]
  print(TrueBeta.k)
  write.table(TrueBeta.k, file_conn, row.names=TRUE, col.names=TRUE)
}

close(file_conn)

# Save simulated data
saveRDS(X, "../Data/Simulation Data/X.RData")
saveRDS(Y, "../Data/Simulation Data/Y.RData")
saveRDS(Z, "../Data/Simulation Data/Z.RData")
saveRDS(beta, "../Data/Simulation Data/beta.RData")
saveRDS(lambda, "../Data/Simulation Data/lambda.RData")
saveRDS(pi, "../Data/Simulation Data/pi.RData")



########################### Recall simulated data ########################

# Read in simulated data
TrueX <- readRDS("../Data/Simulation Data/X.RData")            # N x (D+1)
TrueY <- readRDS("../Data/Simulation Data/Y.RData")            # N x M
TrueZ <- readRDS("../Data/Simulation Data/Z.RData")            # M x K
TrueBeta <- readRDS("../Data/Simulation Data/beta.RData")      # K x (D+1)
TrueLambda <- readRDS("../Data/Simulation Data/lambda.RData")  # M x K
TruePi <- readRDS("../Data/Simulation Data/pi.RData")          # N x K

N <- nrow(TrueY)
M <- ncol(TrueY)
K <- ncol(TrueZ)
D <- ncol(TrueX)-1

############################## Run EM algorithm ###############################

# # Load required packages
# library(foreach)
# library(doParallel)
# 
# # Set up parallel processing
# num_cores <- detectCores()
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # Define function to run EM algorithm for a given K
# run_EM <- function(K) {
#   max_iter <- 20
#   ptm <- proc.time()
#   runEM_result <- EM_CLRM(TrueY, TrueX, K, max_iter, tol = 1e-6)
#   proc.time()-ptm
#   saveRDS(runEM_result, file = paste0("~/Dropbox/Bayesian Inference on CLRM/Results/Simulation Results/EMSimulation_K", K, "_N", N, "_M", M, "_D", D, ".RDS"))
# }
# 
# # Run EM algorithm in parallel for different values of K
# K.list <- c(2, 3, 4, 5)
# foreach(K = K.list) %dopar% {
#   run_EM(K)
# }
# 
# # Stop parallel processing
# stopCluster(cl)



########################### Summarize EM result #############################

# K.list <- c(2, 3, 4, 5)
# BIC.list <- c()
# 
# for (K in K.list) {
#   
#   runEM.result <- readRDS(paste0( "~/Dropbox/Bayesian Inference on CLRM/Results/Simulation Results/EMSimulation_K", K,"_N", N, "_M", M, "_D", D, ".RDS" ))
#   
#   # Extract beta, lambda, and Z from last iteration
#   beta <- runEM.result$beta.list[[length(runEM.result$beta.list)]]
#   lambda <- runEM.result$lambda.list[[length(runEM.result$lambda.list)]]
#   Z <- runEM.result$Z.list[[length(runEM.result$Z.list)]]
#   
#   # Plot log-likelihood
#   plot(
#     runEM.result$loglik[1:length(runEM.result$Z.list)],
#     type = "b",
#     main = paste("K =", K),
#     ylab = "log-likelihood"
#   )
#   
#   # Get the cluster membership summary for each mutation
#   M <- ncol(Y)
#   clusterIndex <- rep(0, M)
#   MutationName <- rep("", M)
#   for (m in 1:M) {
#     MutationName[m] <- paste("Mutation ", m)
#     clusterIndex[m] <- which(Z[m, ] == 1)
#   }
#   
# 
#   # Create text file
#   file_name <- paste0("Cluster_", K, "_N", N, "_M", M, "_EMresults.txt")
#   
#   # open a file connection for writing
#   file_conn <- file(file_name, open = "w")
#   
#   ### Write estimated results from EM algorithm to a file
#   for(k in 1:K){
#     writeLines(paste("Cluster ", k), file_conn)
#     writeLines(paste(MutationName[which(clusterIndex==k)]), file_conn)
#     OR <- as.matrix(exp(beta[k,])) # print odds ratio)
#     print(OR)
#     write.table(OR, file_conn, row.names=TRUE, col.names=FALSE)
#   }
#   
#   # close the file connection
#   close(file_conn)
#   
#   # Save BIC
#   BIC.list <- c(BIC.list, runEM.result$BIC)
# }
# 
# 
# plot(K.list, BIC.list, type="b")
# 

############################ Run Gibbs sampler #################################

# Load required packages
library(foreach)
library(doParallel)

# Set up parallel processing
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(4)

# # Define function to run Gibbs sampler for a given K
# source("gibbs_sampler.R")
# run_Gibbs <- function(K) {
#   num_iter <- 5500
#   const_var <- 1
#   scale_param <- 0.05
#   R <- 500   # Initial steps for updating beta_k with fixed proposal
#   p <- 500  # Frequency for updating beta_k with adaptive proposal
#   ptm <- proc.time()
#   Gibbs_result <- Gibbs_CLRM(TrueY, TrueX, K, num_iter, const_var, scale_param, R, p)
#   proc.time()-ptm
#   saveRDS(Gibbs_result, file = paste0("../Results/Simulation Results/2023-08-31 AdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_Iter", num_iter, ".RDS"))
# }

# Tune with mixture proposal
source("mixture_adaptive_gibbs_sampler.R") # With adaptive proposal on beta_k
run_Gibbs <- function(K) {
  num_iter <- 50000
  const_var <- 1
  scale_param <- 0.05  # Scales covariance in the fixed proposal
  R <- 100            # Initial steps for updating beta_k with fixed proposal
  p <- 1              # Frequency for updating beta_k with adaptive proposal
  psi <- 0.05
  chain_number <- 1
  ptm <- proc.time()
  Gibbs_result <- Gibbs_CLRM(Y, X, K, num_iter, const_var, scale_param, R, p, psi)
  proc.time() - ptm
  saveRDS(Gibbs_result, file = paste0("../Results/Simulation Results/MixtureAdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param,"_Chain", chain_number, ".RDS"))
}

run_Gibbs <- function(chain_number) {
  K <- 3
  num_iter <- 50000
  const_var <- 1
  scale_param <- 0.05  # Scales covariance in the fixed proposal
  R <- 100            # Initial steps for updating beta_k with fixed proposal
  p <- 1              # Frequency for updating beta_k with adaptive proposal
  psi <- 0.05
  ptm <- proc.time()
  Gibbs_result <- Gibbs_CLRM(TrueY, TrueX, K, num_iter, const_var, scale_param, R, p, psi)
  proc.time() - ptm
  saveRDS(Gibbs_result, file = paste0("../Results/Simulation Results/SimulationForAdaptation/MixtureAdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param,"_Chain", chain_number, ".RDS"))
}


# Run Gibbs sampler in parallel for different values of K
K.list <- c(2, 3, 4, 5)
foreach(K = K.list, .combine = 'c', .verbose = TRUE) %dopar% {
  run_Gibbs(K)
}

chain.number.list <- c(2,3)
foreach(chain.number = chain.number.list, .combine = 'c', .verbose = TRUE) %dopar% {
  run_Gibbs(chain.number)
}


# Stop parallel processing
stopCluster(cl)



########################### Summarize Gibbs result #############################

#### Find optimal K

K.list <- c(2, 3, 4, 5)
BIC.list <- c()
burnin <- 10000
num_iter <- 50000

for (K in K.list) {
  
  Gibbs_result_1 <- readRDS(paste0("../Results/Simulation Results/Simulation_1/MixtureAdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param,"_Chain1.RDS"))
  Gibbs_result_2 <- readRDS(paste0("../Results/Simulation Results/Simulation_1/MixtureAdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param,"_Chain2.RDS"))
  Gibbs_result_3 <- readRDS(paste0("../Results/Simulation Results/Simulation_1/MixtureAdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param,"_Chain3.RDS"))
  
  # Gibbs_result <- readRDS(paste0("../Results/Simulation Results/2023-08-31 AdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_Iter", num_iter, ".RDS"))
  # 
  # Gibbs_result <- readRDS(paste0("../Results/Simulation Results/2023-08-05 GibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, ".RDS"))
  # 
  selected.iterations <- burnin:num_iter
  
  # Draw trace plot of beta
  for (k in 1:K) {
    # png(paste0("../Figures/RealDataBeta_", k, "_TracePlot_K", K, "_N", N, "_M", M, "_D", D, ".png"), width=3500, height=3500, res=500)
    par(mfrow = c(4, 3))
    for (d in 1:(D+1)) {
      beta_values <- Gibbs_result_1$beta_list[k=1, d, selected.iterations]
      plot(
        beta_values,
        type = "l",
        xlab = "Iteration",
        ylab = bquote(beta[.(k=3) * "," * .(d-1)]),
        ylim = c(min(beta_values)-0.01, max(beta_values)+0.01)
      )
      lines(Gibbs_result_2$beta_list[k=1, d, selected.iterations], col="blue")
      lines(Gibbs_result_3$beta_list[k=3, d, selected.iterations], col="green")
      abline(h = TrueBeta[k=3, d], col = "red", lty=1)
      
      plot(density(Gibbs_result_1$beta_list[k=3, d, selected.iterations]),
           xlab="",
           ylab = bquote(beta[.(k=1) * "," * .(d-1)]),
           main = "")
      lines(density(Gibbs_result_2$beta_list[k=2, d, selected.iterations]), col="blue") # standard mcmc
      lines(density(Gibbs_result_3$beta_list[k=2, d, selected.iterations]), col="green")
      abline(v = quantile(Gibbs_result_1$beta_list[k=3, d, selected.iterations], c(0.025, 0.975)), col="black", lty=2)
      abline(v = quantile(Gibbs_result_2$beta_list[k=2, d, selected.iterations], c(0.025, 0.975)), col="blue", lty=2)
      abline(v = quantile(Gibbs_result_3$beta_list[k=2, d, selected.iterations], c(0.025, 0.975)), col="green", lty=2)
      abline(v = TrueBeta[k=2, d], col = "red", lty=1)
    }
    # dev.off()
  }

  
  # # Draw histogram of beta
  # par(mfrow = c(4, 4))
  # TrueBetaPlot <- rbind(TrueBeta[2,], TrueBeta[1,], TrueBeta[3,])
  # for (k in 1:K) {
  #   for (d in 1:(D+1)) {
  #     beta_values <- Gibbs_result$beta_list[k, d, ]
  #     hist(
  #       beta_values,
  #       xlab="",
  #       main = bquote(beta[.(k) * "," * .(d)]),
  #       freq = FALSE
  #     )
  #     abline(v = TrueBetaPlot[k,d], col="red")
  #     abline(v = mean(beta_values), col="blue")
  #     abline(v = quantile(beta_values, 0.025), col="blue", lty=2)
  #     abline(v = quantile(beta_values, 0.975), col="blue", lty=2)
  #   }
  # }
  
  # Look at acceptance rate of beta_k
  accept.matrix <- apply(Gibbs_result$accept_beta_list, 2, c) # convert array to matrix
  for (k in 1:K){
    print(paste0("beta_", k))
    accept.rate <- sum(accept.matrix[,k])/length(accept.matrix[,k])
    print(accept.rate)
  }
  
  
  
  # Draw trace plot of lambda
  posterior.lambda <- matrix(NA, nrow = M, ncol = K)
  par(mfrow = c(3, 2))
  for (m in 1:M) {
    for (k in 1:K) {
      lambda_values <- Gibbs_result_3$lambda_list[m, k, selected.iterations]
      plot(
        lambda_values,
        type = "l",
        xlab = "Iteration",
        ylab = bquote(lambda[.(m) * "," * .(k)])
      )
      posterior.lambda[m,k] <- mean(lambda_values)
    }
  }
  
  ### Print estimated results from Gibbs sampler
  
  # Get posterior Z and the cluster membership for each mutation
  clusterIndex <- rep(0, M)
  MutationName <- rep("", M)
  posterior.Z <- matrix(0, nrow = M, ncol = K)
  for (m in 1:M) {
    MutationName[m] <- paste("Mutation ", m)
    clusterIndex[m] <- which.max(rowSums(Gibbs_result_3$Z_list[m, , selected.iterations]))
    posterior.Z[m, clusterIndex[m]] <- 1
  }
  
  # Get posterior beta
  posterior.beta <- matrix(NA, nrow = K, ncol = D+1)
  for (k in 1:K) {
    for (d in 1:(D+1)) {
      posterior.beta[k, d] <- mean(Gibbs_result_3$beta_list[k, d, selected.iterations])
    }
  }
  
  
  # Create text file
  file_name <- paste0("../Results/Simulation Results/MixtureAdaptiveGibbsSamplerSimulation_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param,"_Chain3.txt")
  
  # open a file connection for writing
  file_conn <- file(file_name, open = "w")
  
  # Print cluster specific mutations, odds ratio, and 95% CI
  accept.matrix <- apply(Gibbs_result_3$accept_beta_list, 2, c) # convert array to matrix
  for (k in 1:K) {
    writeLines(paste("Cluster ", k), file_conn)
    writeLines(paste(MutationName[which(clusterIndex==k)]), file_conn)
    summaryTable = t(apply(Gibbs_result_3$beta_list[k, , selected.iterations], 1, function(x)
      c(mean(x), quantile(x, c(0.025, 0.975)))))
    # summaryTable = t(apply(Gibbs_result$beta_list[k, , ], 1, function(x)
    #   c(mean(x), quantile(x, c(0.025, 0.975)))))
    print(summaryTable)
    write.table(summaryTable, file_conn, row.names=TRUE, col.names=TRUE)
    
    # Write acceptance rate of beta_k
    accept.rate <- sum(accept.matrix[selected.iterations-1,k])/length(accept.matrix[selected.iterations-1,k])
    print(paste0("Acceptance rate of beta_", k, ":", accept.rate))
    writeLines(paste0("Acceptance rate of beta_", k, ":", accept.rate), file_conn)
  }
  

  # Compute BIC value
  bic <- BIC(TrueY, TrueX, posterior.beta, posterior.Z)
  writeLines(paste("BIC: ", bic), file_conn)
  

  # Store BIC value
  BIC.list <- c(BIC.list, bic)
  

  
  # close the file connection
  close(file_conn)
}

# Plot BIC
plot(K.list, BIC.list, ylab = "BIC", xlab="K", type="b", xaxt="n")
axis(1, at = K.list, labels = K.list)


