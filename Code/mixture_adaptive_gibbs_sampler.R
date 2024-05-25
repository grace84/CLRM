### Gibbs Sampler with Mixture Adaptive MH for CLRM 
# Y_{N×M} : indicator matrix such that y_{n,m} denote the existence of mutation m on n-th 
#           individual genome sequence (0=does not exist; 1=exist);
# X_{N×(D+1)}: the design matrix. For n-th patient, we have x_n = [1, x_{n,1}, x_{n,2}, ..., x_{n,D}]⊤
# β_k : a (D + 1) × 1 vector that stores the logistic regression intercept + coefficients in cluster k
#       for mutations in cluster k, i.e. βk = [βk,0, βk,1, βk,2, ..., βk,D]⊤
# π_{n,k}: indicates the probability of existence of mutation m on the n-th genome sequence given
#          that mutation m belongs to cluster k
# Z_{m,k}: whether or not (1 or 0) mutation m belongs to cluster k
# λ_{m,k} = P (Z_{m,k} = 1) :the probability that the mutation m belongs to cluster k.


### Define the response probability in a logistic regression model
pi <- function(X, beta){ 
  # exp.vals <- exp(X %*% t(beta)) # X: N x (D+1); beta: K x (D+1)
  # pi.matrix <- exp.vals / (1 + exp.vals)
  pi.matrix <- plogis(X %*% t(beta))
  return(pi.matrix) # N x K
}


### Define log-likelihood function p(Y|β, Z)
log_likelihood <- function(Y, X, beta, Z) {
  
  N <- nrow(Y)
  M <- ncol(Y)
  K <- ncol(Z)
  
  # Compute pi using current beta
  pi.matrix <- pi(X, beta)
  
  # Compute log-likelihood
  ptm <- proc.time()
  ll <- 0
  epsilon <- 1e-10  # Avoid -Inf in log
  # for (n in 1:N) {
    for (m in 1:M) {
      for (k in 1:K) {
        ll <- ll + sum(Z[m, k] * (Y[, m] * log(pi.matrix[, k]) + (1 - Y[, m]) * log((1-pi.matrix[, k]) + epsilon)))
       }
    }
  # }
  proc.time() - ptm
  
  return(ll)
}


# Compute proposal density of beta_k
log_proposal <- function(beta_star, beta_previous, Sigma_hat, Sigma_k, psi, scale_param, D) {
  log_prop <- log((1 - psi) * dmvnorm(beta_star, mean = beta_previous, sigma = (2.38^2 / (D + 1)) * Sigma_hat) +
                    psi * dmvnorm(beta_star, mean = beta_previous, sigma = scale_param * Sigma_k))
  return(log_prop)
}


# Compute posterior distribution of beta_k
log_posterior_beta_k <- function(Y, X, beta_k, mu_k, Sigma_k, Z, M_k, N, k){

  # Start profiling and specify the output file
  # Rprof("profile_output.out")

  # Likelihood
  ptm <- proc.time()
  logit_X_beta <- X %*% beta_k  #X: Nx(D+1); beta_k:(D+1)x1
  log_likelihood <- 0
  for (m in M_k) {
     # for (n in 1:N) {
      log_likelihood <- log_likelihood + sum(Z[m,k]*(Y[,m]*logit_X_beta-log(1+exp(logit_X_beta))))
     # }
  }
  proc.time()-ptm
  
  # Prior
  log_prior <- sum(dmvnorm(beta_k, mean = mu_k, sigma = Sigma_k, log = TRUE))
  
  
  # Stop profiling
  # Rprof(NULL)

  return(log_likelihood + log_prior)
}


# With adaptive step - update covariance every p iterations after initial R steps
update_beta_adaptive <- function(Y, X, beta, mu, Sigma, lambda, Z, 
                                 scale_param, iter, R, 
                                 beta_list, p, prev_emp_Sigma, psi) {

  
  library("mvtnorm")
  N <- nrow(X)
  D <- ncol(X) - 1
  M <- nrow(Z)
  K <- ncol(Z)
  
  # Keep track of acceptance
  acceptance_beta <- matrix(0, nrow=1, ncol=K)
  emp_Sigma <- array(0, c(D+1, D+1, K))
  
  # Loop over clusters k
  for (k in 1:K) {
    
    # Obtain the set of mutations assigned to cluster k
    M_k <- which(Z[, k] == 1)
    
    # If we're past the R steps, adjust the proposal distribution
    if (iter > R) {

      # Generate a random number to decide which Gaussian component to sample from
      if (runif(1) < psi) {
        beta_k_new <- as.vector(rmvnorm(1, mean = beta[k,], sigma = Sigma[,,k] * scale_param))
      } else {
        beta_k_new <- as.vector(rmvnorm(1, mean = beta[k,], sigma = 2.38^2 / (D+1) * prev_emp_Sigma[,,k]))
      }
      
      # Compute acceptance probability
      log_r <- log_proposal(beta[k,], beta_k_new, prev_emp_Sigma[,,k], Sigma[,,k], psi, scale_param, D) + 
        log_posterior_beta_k(Y, X, beta_k_new, mu[,,k], Sigma[,,k], Z, M_k, N, k) -
        log_proposal(beta_k_new, beta[k,], prev_emp_Sigma[,,k], Sigma[,,k], psi, scale_param, D) -
      log_posterior_beta_k(Y, X, beta[k,], mu[,,k], Sigma[,,k], Z, M_k, N, k)

    } else if (iter <= R) {
      
      beta_k_new <- as.vector(rmvnorm(1, mean = beta[k,], sigma = Sigma[,,k] * scale_param))
      
      # Compute acceptance probability
      log_r <- log_posterior_beta_k(Y, X, beta_k_new, mu[,,k], Sigma[,,k], Z, M_k, N, k) -
      log_posterior_beta_k(Y, X, beta[k,], mu[,,k], Sigma[,,k], Z, M_k, N, k)
      
      }
  

    # Accept or reject the proposed beta_k
    if (log(runif(1)) < log_r) {
      beta[k,] <- beta_k_new
      acceptance_beta[1,k] <- 1
    }else{
      beta[k,] <- beta[k,]
      acceptance_beta[1,k] <- 0
    }
    
    # Update empirical covariance matrix every iteration after R steps
    beta_list[k,,iter] <- beta[k,]
    if (iter == R){
      emp_Sigma[,,k] <- cov(t(beta_list[k,,]))
    }else if (iter > R & iter %% p == 0){    # iter between R and p (p>R) seems to produce error
      emp_Sigma[,,k] <- cov(t(beta_list[k,,]))
    }else{
      emp_Sigma[,,k] <- prev_emp_Sigma[,,k]
    }
  }
  
  
  return(list("new.beta" = beta,
              "Acceptance.beta" = acceptance_beta,
              "emp_Sigma" = emp_Sigma))
}



### Define function to update lambda and Z
update_lambda_Z <- function(Y, X, beta, lambda) {
  
  N <- nrow(Y)
  M <- ncol(Y)
  K <- nrow(beta)
  pi.matrix <- pi(X, beta)
  
  ### Update  p(Z_{m,k} = 1|y_m, θ(t))
  log.prob.Zmk <- matrix(NA, nrow=M, ncol=K)
  for (m in 1:M){
    for (k in 1:K){
      
      # p(y_m|Z_{m,k} = 1, θ(t))
      epsilon <- 1e-10
      # for (n in 1:N){
        log.prob.ym <- sum(Y[,m]*log(pi.matrix[,k]) + 
          (1-Y[,m])*log((1-pi.matrix[,k]) + epsilon))
      # } 
      # p(Z_{m,k} = 1, y_m|θ^(t))
      log.prob.Zmk[m,k] <- log(lambda[m,k]) + log.prob.ym
    }
  }
  
  # Estimated λ^(t+1)_{m,k} = p(Z_{m,k} = 1|y_m, θ(t))
  prob.Zmk.given.y <- matrix(NA, nrow=M, ncol=K)
  for (m in 1:M){
      scale.param <- 1000  # for a more spread out probability distribution
      max.log.prob.Zmk <- max(log.prob.Zmk[m,])
      prob.Zmk <- exp((log.prob.Zmk[m,] - max.log.prob.Zmk) / scale.param)
      prob.Zmk.given.y[m,] <- prob.Zmk/sum(prob.Zmk)
      # prob.Zmk.given.y[m,k] <- exp(log.prob.Zmk[m,k] - log(sum(exp(log.prob.Zmk[m,]))))
  }
  
  # Update Z
  Z <- matrix(0, nrow=M, ncol=K) # Initialize Z matrix
  for (m in 1:M) {
    k <- which(rmultinom(1, 1, prob.Zmk.given.y[m,])==1) 
    Z[m,k] <- 1
  }
  
  
  return(list("new.lambda" = prob.Zmk.given.y,
              "new.Z" = Z))
  
}


### K-means with constraints on group size
# run the k-means algorithm multiple times with different initial cluster assignments 
# and then choose the best result that satisfies your cluster size constraints.
custom_kmeans <- function(data, centers, min_cluster_size, nstart = 100, max_iter = 100) {
  
  best_fit <- NULL
  best_withinss <- Inf
  
  for (i in 1:nstart) {
    
    # Run k-means with random initial cluster assignments
    fit <- kmeans(data, centers, nstart = 1, iter.max = max_iter)
    
    # Check if the cluster size constraints are satisfied
    cluster_sizes <- table(fit$cluster)
    if (all(cluster_sizes >= min_cluster_size)) {
      
      # Calculate within-cluster sum of squares
      withinss <- fit$tot.withinss
      
      # Update the best fit if it has a lower withinss value
      if (withinss < best_withinss) {
        best_fit <- fit
        best_withinss <- withinss
      }
    }
  }
  
  if (is.null(best_fit)) {
    print("No feasible solution found.")
    best_fit <- fit
  }
  
  return(best_fit)
}


### Define function for Gibbs sampler
Gibbs_CLRM <- function(Y, X, K, num_iter, const_var, scale_param, R, p, psi) {
  
  library("DirichletReg")
  library("MASS")
  
  ### Initialize parameters 
  N <- nrow(Y)
  M <- ncol(Y)
  D <- ncol(X) - 1
  
  ### Store values of beta, lambda, Z at each iteration
  beta_list <- array(0, c(K, D+1, num_iter))
  Z_list <- array(0, c(M, K, num_iter))
  lambda_list <- array(0, c(M, K, num_iter))
  
  ### Store hyperparameters of beta
  mu_list <- array(0, c(1, D+1, K))
  Sigma_list <- array(0, c(D+1, D+1, K))
  emp_Sigma_list <- array(0, c(D+1, D+1, K, num_iter))
  
  ### Store acceptance of beta
  accept_beta_list <- array(0, c(1, K, num_iter-1))
  
  ### Initialize beta
  # Run GLM to obtain initial values of beta
  # beta_init = matrix(0, nrow = M, ncol = D+1)
  # for (m in 1:M){
  #   beta_init[m,] = suppressWarnings(coef(glm(Y[,m]~ X[,-1], family = "binomial", epsilon = 1e-6, maxit = 500)))
  # }
  
  beta_init <- t(apply(Y, 2, function(y) {
             suppressWarnings(coef(glm.fit(X, y, family = binomial())))
         }))
  
  # K-means clustering for mean of beta
  # set.seed(123)
  # (km.fit = kmeans(x = beta_init[,1:(D+1)], centers = K, iter.max = 100, nstart = 100, algorithm = "MacQueen"))  # nstart: how many random sets
  (km.fit = custom_kmeans(data=beta_init[,1:(D+1)], centers=K, min_cluster_size=2))
  
  # Generate initial values of each beta_k
  for (k in 1:K){
    
    mu_list[,,k] <- km.fit$centers[k,] 
    
    if (sum(km.fit$cluster==k) == 1){
      Sigma_list[,,k] <- const_var*diag(D+1)
    }else{
      within_cluster_var <- apply(beta_init[km.fit$cluster == k, ], 2, var)
      Sigma_list[,,k] <- diag(within_cluster_var)
      # within_cluster_sd <- apply(beta_init[km.fit$cluster == k, ], 2, sd)
      # Sigma_list[,,k] <- diag(within_cluster_sd)
    }
    
    beta_list[,,1][k,] <- mvrnorm(1, mu = mu_list[,,k], 
                                  Sigma = Sigma_list[,,k])
  }
  
  
  ### Initialize lambda
  lambda_list[,,1] <- rdirichlet(M, alpha=rep(1, K))
  
  ### Initialize Z matrix
  for (m in 1:M) {
    k <- which(rmultinom(1, 1, lambda_list[,,1][m,])==1)
    Z_list[,,1][m,k] <- 1
  }
  
  for (iter in 2:num_iter) {
    
    # Update lambda and Z
    update_lambda_Z_result <-  update_lambda_Z(Y, X, beta_list[,,iter-1], lambda_list[,,iter-1])
    lambda_list[,,iter] <- update_lambda_Z_result$new.lambda
    Z_list[,,iter] <- update_lambda_Z_result$new.Z
    
    # # Update beta
    # update_beta_result <- update_beta(Y, X, beta_list[,,iter-1], mu_list, Sigma_list, lambda_list[,,iter], Z_list[,,iter], scale_param) 
    # beta_list[,,iter] <- update_beta_result$new.beta
    # accept_beta_list[,,iter-1] <- update_beta_result$Acceptance.beta
    
    # Update beta using adaptive MCMC
    # update_beta_result <- update_beta_adaptive(Y, X, beta_list[,,iter-1], mu_list, Sigma_list, 
    #                                            lambda_list[,,iter], Z_list[,,iter], scale_param, 
    #                                            iter, burnin, beta_list[,,1:(iter-1)])
    # beta_list[,,iter] <- update_beta_result$new.beta
    # accept_beta_list[,,iter-1] <- update_beta_result$Acceptance.beta
    
    # Update beta using adaptive MCMC
    # ptm <- proc.time()
    update_beta_result <- update_beta_adaptive(Y, X, beta_list[,,iter-1], mu_list, Sigma_list,
                                               lambda_list[,,iter], Z_list[,,iter], scale_param,
                                               iter, R, beta_list[,,1:iter], p, emp_Sigma_list[,,,iter-1], psi)
    beta_list[,,iter] <- update_beta_result$new.beta
    accept_beta_list[,,iter-1] <- update_beta_result$Acceptance.beta
    emp_Sigma_list[,,,iter] <- update_beta_result$emp_Sigma
    # proc.time()-ptm
    
    # Print number every 10 iterations
    if (iter %% 10 == 0) {
      print(paste("Iteration:", iter))
    }
    
  }
  
  return(list("beta_list" = beta_list,
              "lambda_list" = lambda_list,
              "Z_list" = Z_list,
              "mu_list" = mu_list,
              "sigma_list" = Sigma_list,
              "accept_beta_list" = accept_beta_list,
              "emp_Sigma_list" = emp_Sigma_list))                 
}


### Compute AIC = -2 * log-likelihood + 2*number of parameters in the model
AIC <- function(Y, X, beta, Z){
  
  N <- nrow(Y)
  M <- ncol(Y)
  D <- ncol(X) - 1
  K <- ncol(Z)
  
  AIC.value <- -2 * log_likelihood(Y, X, beta, Z) + 2*(K*(D+1) + (K-1)*M)
  
  return(AIC.value)
}


### Compute BIC = -2 * log-likelihood + number of parameters in the model * log(number of observations)
BIC <- function(Y, X, beta, Z){
  
  N <- nrow(Y)
  M <- ncol(Y)
  D <- ncol(X) - 1
  K <- ncol(Z)
  # 
  # Count the number of empty clusters
  # empty_clusters <- sum(colSums(Z) == 0)

  BIC.value <- -2 * log_likelihood(Y, X, beta, Z) + (K*(D+1) + (K-1)*M)*log(N)
  # BIC.value <- -2 * log_likelihood(Y, X, beta, Z) + (K*(D+1) + K*M)*log(N) + empty_clusters*log(N)
  
  return(BIC.value)
}

# BIC.list.new <- matrix(NA, 1, 5)
# for (k in K.list){
#   BIC.list.new[1,k-1] <- BIC.list[k-1]-(k*(D+1) + k*M)*log(N) + (k*(D+1) + k*M)*log(N)*6
# }



# Define the complete data log-likelihood function
comploglik <- function(Y, X, Z, beta, lambda) {
  
  N <- nrow(X)
  M <- nrow(Z)
  K <- ncol(Z)
  
  loglik <- 0
  
  # log p(Z|λ)
  lambda[is.infinite(log(lambda))] <- 1e-308 # replace 'Inf' in log(lambda) with a very small number
  
  for (m in 1:M) {
    for (k in 1:K){
      loglik <- loglik + Z[m,k]*log(lambda[m,k])
    } 
  }
  
  # log p(y|Z, β)
  for (m in 1:M) {
    for (n in 1:N) {
      for (k in 1:K) {
        loglik <- loglik + Z[m,k]*(Y[n,m]*(t(X[n,]) %*% beta[k,])-log(1+exp(t(X[n,]) %*% beta[k,])))
      }
    }
  }
  
  return(loglik)
}


## Define the ICL criterion for CLRM
ICL <- function(Y, X, beta, Z, lambda) {
  
  N <- nrow(Y)
  M <- ncol(Y)
  D <- ncol(X) - 1
  K <- ncol(Z)
  
  
  ICL.value <- comploglik(Y, X, Z, beta, lambda) - (K*(D+1) + (K-1)*M)*log(N)/2  # Calculate the ICL criterion
  
  return(ICL.value)
}



