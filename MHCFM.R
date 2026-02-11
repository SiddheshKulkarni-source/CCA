
############################################################
### MHCFM script
### 
###
### High-level model idea:
###  - Two-view latent factor model:
###      X1 = mu1 + A1 Z + eps1
###      X2 = mu2 + A2 Z + eps2
###  - Z ~ N(0, I_d)
###  - eps blocks have covariance Phi1, Phi2 (can be diagonal IG updates)
###    and optionally Graphical Horseshoe (GHS) updates early on
###  - Multiplicative shrinkage on columns of A via eta (MGP-type)
###  - Posterior CCA is computed from implied grand covariance:
###      Sigma_grand = [[A1A1' + Phi1,  A1A2'],
###                     [A2A1',         A2A2' + Phi2]]
###  - Store thinned samples of parameters + CCA summaries
############################################################

##############################
### Libraries
##############################
library(MASS)     # mvrnorm
library(mvtnorm)  # dmvnorm
library(psych)    # fa() initialization

##############################
### Helper Functions
##############################

### Inverse-Gamma sampler (1 / Gamma)
rinvgamma <- function(n, shape, rate = 1, scale = 1/rate)
{
  if (missing(rate) && !missing(scale))
    rate <- 1/scale
  1/rgamma(n, shape, rate)
}

############################################################
### Build grand covariance from two views:
### sigma_grand =
###   [ A1 A1' + phi1      A1 A2' ]
###   [ A2 A1'             A2 A2' + phi2 ]
############################################################
sigma_grand_cal <- function(A_1, A_2, phi_1, phi_2, p1, p2)
{
  m_11 <- A_1 %*% t(A_1) + phi_1
  m_12 <- A_1 %*% t(A_2)
  m_21 <- A_2 %*% t(A_1)
  m_22 <- A_2 %*% t(A_2) + phi_2
  
  sigma_grand <- matrix(nrow = p1 + p2, ncol = p1 + p2)
  sigma_grand[1:nrow(m_11), 1:ncol(m_11)] <- m_11
  sigma_grand[1:nrow(m_12), ncol(m_11) + 1:ncol(m_12)] <- m_12
  sigma_grand[nrow(m_11) + 1:nrow(m_21), 1:ncol(m_21)] <- m_21
  sigma_grand[nrow(m_11) + 1:nrow(m_22), ncol(m_21) + 1:ncol(m_22)] <- m_22
  
  return(sigma_grand)
}

############################################################
### CCA computed from a grand covariance matrix
### Inputs:
###  - Covariance_matrix_Grand: covariance of concatenated (X1, X2)
###  - p1, p2: view dimensions
###  - CCA_select: number of canonical correlations to return
### Output:
###  - canonical correlations + canonical direction vectors
############################################################
CC_function_grand_2 <- function(Covariance_matrix_Grand, p1, p2, CCA_select)
{
  sigma_grand <- Covariance_matrix_Grand
  
  # Block submatrices
  m_11 <- sigma_grand[1:p1, 1:p1]           # Var(X1)
  m_12 <- sigma_grand[1:p1, p1+1:p2]        # Cov(X1, X2)  (NOTE: indexing as in your code)
  m_21 <- sigma_grand[p1+1:p2, 1:p1]        # Cov(X2, X1)
  m_22 <- sigma_grand[p1+1:p2, p1+1:p2]     # Var(X2)
  
  # Matrix inverse square roots via eigendecomposition
  eigen_x <- eigen(m_11)
  eigen_y <- eigen(m_22)
  
  sqrm_11 <- (eigen_x$vectors) %*% diag(eigen_x$values^-0.5) %*% solve(eigen_x$vectors)
  sqrm_22 <- (eigen_y$vectors) %*% diag(eigen_y$values^-0.5) %*% solve(eigen_y$vectors)
  
  # CCA eigen problems
  m1 <- sqrm_11 %*% m_12 %*% solve(m_22) %*% t(m_12) %*% sqrm_11
  m2 <- sqrm_22 %*% m_21 %*% solve(m_11) %*% t(m_21) %*% sqrm_22
  
  m1e <- (eigen(m1))
  m2e <- (eigen(m2))
  
  # Canonical correlations computed from eigenvectors
  Canonical_Correaltions <- c()
  for (j in 1:CCA_select) {
    Canonical_Correaltions[j] <-
      Re(m1e$vectors[, j]) %*% sqrm_11 %*% m_12 %*% sqrm_22 %*% Re(m2e$vectors[, j])
  }
  
  # Sign correction (keep CC positive; flip Y-direction accordingly)
  for (i in 1:length(Canonical_Correaltions)) {
    if (sign(Canonical_Correaltions[i]) == -1) {
      Canonical_Correaltions[i] <- abs(Canonical_Correaltions[i])
      m2e$vectors[, i] <- -1 * m2e$vectors[, i]
    }
  }
  
  return(list(
    CCA_grand = Canonical_Correaltions,
    Direction_CCA_Vec1_grand = Re(m1e$vectors[, 1:CCA_select]),
    Direction_CCA_Vec2_grand = Re(m2e$vectors[, 1:CCA_select])
  ))
}

############################################################
### Log-likelihood calculator for grand MVN
### Inputs:
###  - X_grand: observed concatenated data matrix
###  - mu_1, mu_2: mean vectors (view-wise)
###  - phi_grand: grand covariance matrix
### Output:
###  - vector of per-sample log-likelihoods
############################################################
lk_para_calculator2 <- function(X_grand, mu_1, mu_2, phi_grand)
{
  mu_grand <- rbind(mu_1, mu_2)         # concatenate means
  likelihood <- dmvnorm(X_grand, mean = mu_grand, sigma = phi_grand, log = TRUE)
  return(likelihood)
}

### simple norm helper (not used heavily in this snippet)
norm_vec <- function(x) sqrt(sum(x^2))

############################################################
### Graphical Horseshoe (GHS) update for precision matrices
### This is a node-wise Gibbs update scheme:
###  - Updates Omega (precision) and Sigma (covariance)
###  - Shrinkage controlled by lambda_sq (local) and tau_sq (global)
### NOTE: In MHCFM below, you use diagonal IG updates for Phi after iter>1,
###       and use this GHS routine only in the "else" branch (iter <= 1).
############################################################
GHS_modified_2 <- function(X, Lambda_sq, Nu, Omega, Sigma, tau_sq, xi) {
  
  n <- nrow(X)
  S <- t(X) %*% (X)      # scatter matrix
  p <- nrow(S)
  
  # Precompute indices excluding i for node-wise updates
  ind_all <- matrix(0, p - 1, p)
  for (i in 1:p) {
    if (i == 1) {
      ind = c(2:p)
    } else if (i == p) {
      ind = c(1:(p - 1))
    } else {
      ind = c(1:(i - 1), (i + 1):p)
    }
    ind_all[, i] <- ind
  }
  
  # Node-wise updates
  for (i in 1:p) {
    
    ind <- ind_all[, i]
    
    # Partition Sigma and S for node i
    sigma_11 <- Sigma[ind, ind]
    sigma_12 <- Sigma[ind, i]
    sigma_22 <- Sigma[i, i]
    s_21 <- S[ind, i]
    s_22 <- S[i, i]
    
    # Local shrinkage params for edges connected to i
    lambda_sq_12 <- Lambda_sq[ind, i]
    nu_12 <- Nu[ind, i]
    
    # Sample gamma (related to omega_22)
    gamma <- rgamma(1, shape = ((n/2) + 1), rate = s_22/2)
    
    # Conditional block inversion quantities
    inv_Omega_11 <- sigma_11 - ((sigma_12 %*% t(sigma_12)) / (sigma_22))
    
    # Posterior precision for beta = omega_12
    inv_C <- (s_22 * inv_Omega_11) + diag(1 / (lambda_sq_12 * tau_sq), p - 1)
    inv_C_chol <- chol(inv_C)
    
    # Posterior mean of beta
    mu_i <- -chol2inv(inv_C_chol) %*% s_21
    
    # Sample beta and build omega blocks
    beta <- mu_i + solve(inv_C_chol, rnorm(p - 1, 0, 1))
    omega_12 <- beta
    omega_22 <- gamma + t(beta) %*% inv_Omega_11 %*% beta
    
    # Update local shrinkage lambda_sq and nu
    rate <- (1/nu_12) + ((omega_12^2) / (2 * tau_sq))
    lambda_sq_12 <- rinvgamma(p - 1, shape = 1, rate = rate)
    nu_12 <- rinvgamma(p - 1, shape = 1, rate = (1 + (1 / lambda_sq_12)))
    
    # Write back Omega
    Omega[i, ind] <- omega_12
    Omega[ind, i] <- omega_12
    Omega[i, i] <- omega_22
    
    # Update Sigma from Omega blocks
    temp <- inv_Omega_11 %*% beta
    Sigma_11 <- inv_Omega_11 + temp %*% t(temp) / gamma
    sigma_12 <- -temp / gamma
    sigma_22 <- 1 / gamma
    
    Sigma[ind, ind] <- Sigma_11
    Sigma[i, i] <- sigma_22
    Sigma[i, ind] <- sigma_12
    Sigma[ind, i] <- sigma_12
    
    # Store shrinkage params symmetrically
    Lambda_sq[i, ind] <- lambda_sq_12
    Lambda_sq[ind, i] <- lambda_sq_12
    Nu[i, ind] <- nu_12
    Nu[ind, i] <- nu_12
  }
  
  # Global shrinkage update tau_sq
  omega_vector <- Omega[lower.tri(Omega, diag = FALSE)]
  lambda_sq_vector <- Lambda_sq[lower.tri(Lambda_sq, diag = FALSE)]
  rate1 <- 1/xi + sum((omega_vector^2) / (2 * lambda_sq_vector))
  
  tau_sq <- rinvgamma(1, shape = ((p * (p - 1) / 2) + 1) / 2, rate = rate1)
  xi <- rinvgamma(1, shape = 1, rate = (1 + 1/tau_sq))
  
  return(list(Omega = Omega, Lambda_sq = Lambda_sq, Nu = Nu,
              Tau_sq = tau_sq, Sigma = Sigma, xi = xi))
}

############################################################
###  MHCFM run function
### Inputs:
###   d          : latent factor dimension
###   thin       : thinning gap (store every 'thin' iters after burn)
###   burn_iter  : burn-in iterations
###   mcmc_iter  : sampling iterations after burn
###   zeta       : hyperparameter controlling eta prior (see Epost_iter update)
###   X_1, X_2   : data matrices (n x p1), (n x p2)
###   CCA_select : number of canonical correlations to compute
############################################################
MHCFM <- function(d,
                 thin,
                 burn_iter,
                 mcmc_iter,
                 zeta,
                 X_1,
                 X_2,
                 CCA_select)
{
  set.seed(666)
  
  ##########################################################
  ### Dimensions / iteration setup
  ##########################################################
  p1 <- ncol(X_1)
  p2 <- ncol(X_2)
  no_patients <- nrow(X_1)
  
  iter_max <- burn_iter + mcmc_iter
  
  ##########################################################
  ### Allocate full-chain storage (some are used, others legacy)
  ##########################################################
  A_1 <- array(dim = c(p1, d))    # loadings view1
  A_2 <- array(dim = c(p2, d))    # loadings view2
  mu_1 <- array(dim = c(1, p1))   # mean view1
  mu_2 <- array(dim = c(1, p2))   # mean view2
  
  phi_post_1 <- array(dim = c(p1, p1))   # covariance view1 (saved)
  phi_post_2 <- array(dim = c(p2, p2))   # covariance view2 (saved)
  
  Omega_post_1 <- array(dim = c(p1, p1)) # precision view1
  Omega_post_2 <- array(dim = c(p2, p2)) # precision view2
  
  Phi_grand_det <- array(dim = c(1, iter_max))
  log_det_data <- array(dim = c(1, iter_max))
  log_det_MCMC <- array(dim = c(1, iter_max))
  Log_likelihood_MCMC_Grand <- array(dim = c(no_patients, iter_max))
  
  tausqpost_1 <- array(dim = c(1, iter_max)) # global shrinkage A1
  tausqpost_2 <- array(dim = c(1, iter_max)) # global shrinkage A2
  
  ##########################################################
  ### Thinning schedule and thinned storage
  ##########################################################
  thin_seq <- seq(burn_iter + thin, iter_max, thin)
  
  A_1_thin <- array(dim = c(p1, d, length(thin_seq)))
  A_2_thin <- array(dim = c(p2, d, length(thin_seq)))
  mu_1_thin <- array(dim = c(1, p1, length(thin_seq)))
  mu_2_thin <- array(dim = c(1, p2, length(thin_seq)))
  
  phi_post_1_thin <- array(dim = c(p1, p1, length(thin_seq)))
  phi_post_2_thin <- array(dim = c(p2, p2, length(thin_seq)))
  
  Omega_post_1_thin <- array(dim = c(p1, p1, length(thin_seq)))
  Omega_post_2_thin <- array(dim = c(p2, p2, length(thin_seq)))
  
  log_det_MCMC_thin <- array(dim = c(1, length(thin_seq)))
  Log_likelihood_MCMC_Grand_thin <- array(dim = c(no_patients, length(thin_seq)))
  
  CCA_thin <- array(dim = c(1, CCA_select, length(thin_seq)))
  Direction_CCA_Vec1_thin <- array(dim = c(p1, CCA_select, length(thin_seq)))
  Direction_CCA_Vec2_thin <- array(dim = c(p2, CCA_select, length(thin_seq)))
  
  ##########################################################
  ### Shrinkage + latent factor bookkeeping
  ##########################################################
  eta <- array(dim = c(1, d, iter_max))
  eta_iter <- c()
  
  # Placeholders / initialization for GHS objects
  Omega_1_GHS <- Sigma_1_GHS <- Phi_1_iter <- array(dim = c(p1, p1))
  Omega_2_GHS <- Sigma_2_GHS <- Phi_2_iter <- array(dim = c(p2, p2))
  Phi_grand_iter <- array(dim = c(p1 + p2, p1 + p2))
  
  Lambda_sq_GHS_1 <- Nu_GHS_1 <- array(dim = c(p1, p1))
  Lambda_sq_GHS_2 <- Nu_GHS_2 <- array(dim = c(p2, p2))
  
  tau_sq_GHS_1 <- array(dim = c(1))
  tau_sq_GHS_2 <- array(dim = c(1))
  xi_GHS_1 <- array(dim = c(1))
  xi_GHS_2 <- array(dim = c(1))
  
  Z_post_iter <- array(dim = c(no_patients, d))
  
  # residual and conditional regression buffers
  X_delta_1 <- matrix(nrow = no_patients, ncol = p1)
  X_delta_2 <- matrix(nrow = no_patients, ncol = p2)
  X_tilda_1_iter <- array(dim = c(no_patients, p1))
  X_tilda_2_iter <- array(dim = c(no_patients, p2))
  
  A_1_iter <- array(dim = c(p1, d))
  A_2_iter <- array(dim = c(p2, d))
  
  Cpostmat_1_iter <- array(dim = c(p1, d))
  Cpostmat_2_iter <- array(dim = c(p2, d))
  
  lambdasqpostmat_1_iter <- array(dim = c(p1, d))
  lambdasqpostmat_2_iter <- array(dim = c(p2, d))
  
  Dpost_1_iter <- c()
  Dpost_2_iter <- c()
  tausqpost_1_iter <- c()
  tausqpost_2_iter <- c()
  Epost_iter <- c()
  
  ##########################################################
  ##### Initialization via factor analysis
  ##########################################################
  factor_analysis <- fa(cbind(X_1, X_2),
                        nfactors = d,
                        n.obs = no_patients,
                        rotate = "none",
                        scores = "regression")
  
  # Initialize latent factors (regression scores)
  Z_post_iter <- factor_analysis$scores
  
  # Initialize loadings split by view
  A_1_iter <- factor_analysis$loadings[1:p1, ]
  A_2_iter <- factor_analysis$loadings[p1+1:p2, ]   # (kept exactly as your code)
  
  ##########################################################
  ### Initialize means
  ##########################################################
  mu_1_iter <- apply(X_1, 2, mean)
  mu_2_iter <- apply(X_2, 2, mean)
  
  ##########################################################
  ### Initialize Phi_1 and Phi_2 (start as diagonal)
  ##########################################################
  Phi_1_iter <- 0.5 * diag(p1)
  Phi_2_iter <- 0.5 * diag(p2)
  
  # Grand covariance under current parameters
  Phi_grand_iter <- sigma_grand_cal(A_1_iter, A_2_iter, Phi_1_iter, Phi_2_iter, p1, p2)
  
  sigma <- 100
  
  # Initialize Omega as inverse of Phi (diagonal => inverse is easy)
  Omega_1_GHS <- solve(Phi_1_iter)
  Omega_2_GHS <- solve(Phi_2_iter)
  
  # Initialize log determinant chain
  log_det_MCMC[, 1] <- 0
  
  # Initialize eta and local shrinkage matrices
  eta_iter <- rep(1, d)
  
  tausqpost_1_iter <- (mean(apply(t(A_1_iter)^2, 2, mean))) / eta_iter[1]
  lambdasqpostmat_1_iter <- matrix(1, p1, d)
  
  tausqpost_2_iter <- (mean(apply(t(A_2_iter)^2, 2, mean))) / eta_iter[1]
  lambdasqpostmat_2_iter <- matrix(1, p2, d)
  
  # Initialize GHS hyperparameters
  Lambda_sq_GHS_1 <- Nu_GHS_1 <- matrix(1, nrow = p1, ncol = p1)
  Lambda_sq_GHS_2 <- Nu_GHS_2 <- matrix(1, nrow = p2, ncol = p2)
  
  tau_sq_GHS_1 <- 1; tau_sq_GHS_2 <- 1
  xi_GHS_1 <- 1;    xi_GHS_2 <- 1
  
  Phi_grand_det[, 1] <- 1
  
  ##########################################################
  ### Precompute grand data and its mean
  ##########################################################
  X_grand <- cbind(X_1, X_2)
  X_grand_mean <- apply(X_grand, 2, mean, na.rm = TRUE)
  
  ##########################################################
  ### Main Gibbs loop
  ##########################################################
  i <- 1
  for (iter in 2:iter_max) {
    
    print(iter)
    
    ########################################################
    ### Copy previous state into convenience aliases
    ########################################################
    Phi_1 <- Phi_1_iter
    Phi_2 <- Phi_2_iter
    A1 <- A_1_iter
    A2 <- A_2_iter
    eta_prod <- cumprod(eta_iter)
    
    ########################################################
    ### 1) Update grand mean mu | Phi_grand, X
    ########################################################
    E_mu_cal <- no_patients * (solve(Phi_grand_iter)) + diag(1/sigma, p1 + p2)
    variance_mu <- solve(E_mu_cal)
    mean_mu <- variance_mu %*% (solve(Phi_grand_iter)) %*% (X_grand_mean)
    mu_grand <- mvrnorm(1, mean_mu, variance_mu)
    
    mu_1_iter <- mu_grand[1:p1]
    mu_2_iter <- mu_grand[-c(1:p1)]
    
    ########################################################
    ### 2) Update latent factors Z_i | A, Omega, X
    ########################################################
    E <- diag(d) + t(A1) %*% Omega_1_GHS %*% A1 + t(A2) %*% Omega_2_GHS %*% A2
    Sigma_n <- solve(E)
    
    for (patient in 1:no_patients) {
      Mu_n <- Sigma_n %*% (
        t(A1) %*% Omega_1_GHS %*% (X_1[patient, ] - mu_1_iter) +
          t(A2) %*% Omega_2_GHS %*% (X_2[patient, ] - mu_2_iter)
      )
      Z_post_iter[patient, ] <- mvrnorm(1, Mu_n, Sigma_n)
    }
    
    ########################################################
    ### 3) Update loadings A_1 row-wise
    ########################################################
    for (j in 1:ncol(X_1)) {
      
      A_j <- A_1_iter[-j, ]
      mu_1_j <- mu_1_iter[-j]
      
      # Conditional regression coefficient
      phi_inversion_1 <- -Omega_1_GHS[-j, j] / Omega_1_GHS[j, j]
      
      # Construct tilde outcome (conditional residual)
      for (patient in 1:no_patients) {
        X_tilda_1_iter[patient, j] <-
          X_1[patient, j] - mu_1_iter[j] -
          phi_inversion_1 %*% (X_1[patient, -j] - mu_1_j - A_j %*% Z_post_iter[patient, ])
      }
      
      # conditional variance term
      phi_j_diag_1 <- 1 / Omega_1_GHS[j, j]
      
      # prior variances from global/local shrinkage
      lambda_1_inv <- c()
      for (k in 1:d) {
        lambda_1_inv[k] <- as.numeric(tausqpost_1_iter * (eta_prod[k]) * lambdasqpostmat_1_iter[j, k])
      }
      
      # guard tiny variances for numerical stability
      lambda_1_inv <- ifelse(lambda_1_inv < 1e-20, 1e-20, lambda_1_inv)
      
      # posterior precision and draw
      inv_C_1 <- phi_j_diag_1 * t(Z_post_iter) %*% (Z_post_iter) + diag(1 / lambda_1_inv)
      inv_C_chol_1 <- chol(inv_C_1)
      mu_i_1 <- phi_j_diag_1 * chol2inv(inv_C_chol_1) %*% t(Z_post_iter) %*% (X_tilda_1_iter[, j])
      A_1_iter[j, ] <- mu_i_1 + solve(inv_C_chol_1, rnorm(d))
    }
    
    ########################################################
    ### 4) Update loadings A_2 row-wise
    ########################################################
    for (j in 1:ncol(X_2)) {
      
      A_j <- A_2_iter[-j, ]
      mu_2_j <- mu_2_iter[-j]
      
      phi_inversion_2 <- -Omega_2_GHS[-j, j] / Omega_2_GHS[j, j]
      
      for (patient in 1:no_patients) {
        X_tilda_2_iter[patient, j] <-
          X_2[patient, j] - mu_2_iter[j] -
          phi_inversion_2 %*% (X_2[patient, -j] - mu_2_j - A_j %*% Z_post_iter[patient, ])
      }
      
      lambda_2_inv <- c()
      for (k in 1:d) {
        lambda_2_inv[k] <- as.numeric(tausqpost_2_iter * (eta_prod[k]) * lambdasqpostmat_2_iter[j, k])
      }
      
      lambda_2_inv <- ifelse(lambda_2_inv < 1e-20, 1e-20, lambda_2_inv)
      
      phi_j_diag_2 <- 1 / Omega_2_GHS[j, j]
      inv_C_2 <- phi_j_diag_2 * t(Z_post_iter) %*% (Z_post_iter) + diag(1 / lambda_2_inv)
      inv_C_chol_2 <- chol(inv_C_2)
      mu_i_2 <- phi_j_diag_2 * chol2inv(inv_C_chol_2) %*% t(Z_post_iter) %*% (X_tilda_2_iter[, j])
      A_2_iter[j, ] <- mu_i_2 + solve(inv_C_chol_2, rnorm(d))
    }
    
    ########################################################
    ### 5) Rescaling step:
    ### Normalize latent factor scales to reduce identifiability issues:
    ###   - RESCALE uses sd of Z columns
    ###   - multiply A by RESCALE and divide Z by RESCALE
    ########################################################
    RESCALE <- sqrt(diag(cov(Z_post_iter)))
    A_1_iter <- A_1_iter %*% diag(RESCALE)
    A_2_iter <- A_2_iter %*% diag(RESCALE)
    Z_post_iter <- Z_post_iter %*% diag(1/RESCALE)
    
    ########################################################
    ### 6) Residual calculation for Phi / Omega updates
    ########################################################
    for (patient in 1:no_patients) {
      X_delta_1[patient, ] <- X_1[patient, ] - mu_1_iter - A1 %*% Z_post_iter[patient, ]
      X_delta_2[patient, ] <- X_2[patient, ] - mu_2_iter - A2 %*% Z_post_iter[patient, ]
    }
    
    ########################################################
    ### 7) Residual covariance update:
    ### If iter > 1: update only diagonal Phi via IG (fast, stable)
    ### Else: do one-step full Graphical Horseshoe update
    ########################################################
    if (iter > 1) {
      
      for (j in 1:p1) {
        Phi_1_iter[j, j] <- rinvgamma(
          1, 0.1 + 0.5 * no_patients,
          0.1 + 0.5 * sum(X_delta_1[, j]^2)
        )
      }
      
      for (j in 1:p2) {
        Phi_2_iter[j, j] <- rinvgamma(
          1, 0.1 + 0.5 * no_patients,
          0.1 + 0.5 * sum(X_delta_2[, j]^2)
        )
      }
      
      # update precision diagonals
      diag(Omega_1_GHS) <- 1 / diag(Phi_1_iter)
      diag(Omega_2_GHS) <- 1 / diag(Phi_2_iter)
      
    } else {
      
      # Full GHS update for view 1
      GHS_Run_1 <- GHS_modified_2(
        X = X_delta_1,
        Lambda_sq = Lambda_sq_GHS_1,
        Nu = Nu_GHS_1,
        Omega = Omega_1_GHS,
        Sigma = Phi_1,
        tau_sq = tau_sq_GHS_1,
        xi = xi_GHS_1
      )
      Phi_1_iter <- GHS_Run_1$Sigma
      Omega_1_GHS <- GHS_Run_1$Omega
      Lambda_sq_GHS_1 <- GHS_Run_1$Lambda_sq
      Nu_GHS_1 <- GHS_Run_1$Nu
      tau_sq_GHS_1 <- GHS_Run_1$Tau_sq
      xi_GHS_1 <- GHS_Run_1$xi
      
      # Full GHS update for view 2
      GHS_Run_2 <- GHS_modified_2(
        X = X_delta_2,
        Lambda_sq = Lambda_sq_GHS_2,
        Nu = Nu_GHS_2,
        Omega = Omega_2_GHS,
        Sigma = Phi_2,
        tau_sq = tau_sq_GHS_2,
        xi = xi_GHS_2
      )
      Phi_2_iter <- GHS_Run_2$Sigma
      Omega_2_GHS <- GHS_Run_2$Omega
      Lambda_sq_GHS_2 <- GHS_Run_2$Lambda_sq
      Nu_GHS_2 <- GHS_Run_2$Nu
      tau_sq_GHS_2 <- GHS_Run_2$Tau_sq
      xi_GHS_2 <- GHS_Run_2$xi
    }
    
    ########################################################
    ### 8) Update eta (multiplicative shrinkage) with zeta prior
    ### Epost_iter provides an auxiliary IG update involving zeta
    ########################################################
    A1.sums <- apply(A_1_iter^2 / lambdasqpostmat_1_iter / tausqpost_1_iter, 2, sum) / 2
    A2.sums <- apply(A_2_iter^2 / lambdasqpostmat_2_iter / tausqpost_2_iter, 2, sum) / 2
    
    for (k in 2:d) {
      
      eta_prod_temp <- cumprod(eta_iter)
      eta_prod_temp[k:d] <- eta_prod_temp[k:d] / eta_iter[k]
      
      # auxiliary IG update using zeta
      Epost_iter[k] <- rinvgamma(1, 1, ((1/zeta^2) + (1/eta_iter[k])))
      
      # eta_k update
      eta_iter[k] <- rinvgamma(
        1,
        0.5 * (p1 + p2) * (d - k + 1),
        1/Epost_iter[k] +
          sum(A1.sums[k:d] / eta_prod_temp[k:d]) +
          sum(A2.sums[k:d] / eta_prod_temp[k:d])
      )
    }
    
    eta_prod <- cumprod(eta_iter)
    
    ########################################################
    ### 9) Update global shrinkage tau^2 for A_1 and A_2
    ########################################################
    A1.sums <- apply(A_1_iter^2 / lambdasqpostmat_1_iter, 2, sum) / 2
    A2.sums <- apply(A_2_iter^2 / lambdasqpostmat_2_iter, 2, sum) / 2
    A1.sums <- A1.sums / eta_prod
    A2.sums <- A2.sums / eta_prod
    
    Dpost_1_iter <- rinvgamma(1, 1, 1 + (1/tausqpost_1_iter))
    tausqpost_1_iter <- rinvgamma(1, (((p1 * d) + 1) / 2),
                                  (1/Dpost_1_iter) + sum(A1.sums))
    
    Dpost_2_iter <- rinvgamma(1, 1, 1 + (1/tausqpost_2_iter))
    tausqpost_2_iter <- rinvgamma(1, (((p2 * d) + 1) / 2),
                                  (1/Dpost_2_iter) + sum(A2.sums))
    
    # cap at 1 (your stabilization rule)
    tausqpost_1_iter <- min(tausqpost_1_iter, 1)
    tausqpost_2_iter <- min(tausqpost_2_iter, 1)
    
    ########################################################
    ### 10) Update local shrinkage lambda^2 for each loading entry
    ########################################################
    for (j in 1:p1) {
      for (k in 1:d) {
        Cpostmat_1_iter[j, k] <- rinvgamma(1, 1, 1 + (1/lambdasqpostmat_1_iter[j, k]))
        lambdasqpostmat_1_iter[j, k] <- rinvgamma(
          1, 1,
          (1/Cpostmat_1_iter[j, k]) +
            ((A_1_iter[j, k])^2 / (2 * (eta_prod[k]) * tausqpost_1_iter))
        )
      }
    }
    
    for (j in 1:p2) {
      for (k in 1:d) {
        Cpostmat_2_iter[j, k] <- rinvgamma(1, 1, 1 + (1/lambdasqpostmat_2_iter[j, k]))
        lambdasqpostmat_2_iter[j, k] <- rinvgamma(
          1, 1,
          (1/Cpostmat_2_iter[j, k]) +
            ((A_2_iter[j, k])^2 / (2 * (eta_prod[k]) * tausqpost_2_iter))
        )
      }
    }
    
    ########################################################
    ### 11) Compute CCA + diagnostics from implied grand covariance
    ########################################################
    Phi_grand_iter <- sigma_grand_cal(A_1_iter, A_2_iter, Phi_1_iter, Phi_2_iter, p1, p2)
    
    cc_cal <- CC_function_grand_2(Phi_grand_iter, p1 = p1, p2 = p2, CCA_select = CCA_select)
    
    log_det_MCMC[, iter] <- determinant(Phi_grand_iter, logarithm = TRUE)$modulus[1]
    
    Log_likelihood_MCMC_Grand[, iter] <-
      lk_para_calculator2(
        X_grand,
        mu_1 = as.matrix(mu_1_iter),
        mu_2 = as.matrix(mu_2_iter),
        phi_grand = Phi_grand_iter
      )
    
    # Store shrinkage chains
    eta[,, iter] <- eta_iter
    tausqpost_1[, iter] <- tausqpost_1_iter
    tausqpost_2[, iter] <- tausqpost_2_iter
    
    ########################################################
    ### 12) Save thinned draws
    ########################################################
    if (iter == thin_seq[i]) {
      
      A_1_thin[,, i] <- A_1_iter
      A_2_thin[,, i] <- A_2_iter
      
      phi_post_1_thin[,, i] <- Phi_1_iter
      phi_post_2_thin[,, i] <- Phi_2_iter
      
      Omega_post_1_thin[,, i] <- Omega_1_GHS
      Omega_post_2_thin[,, i] <- Omega_2_GHS
      
      mu_1_thin[,, i] <- mu_1_iter
      mu_2_thin[,, i] <- mu_2_iter
      
      # compute CCA for this thinned draw using thinned Phi
      Phi_grand_MCMC_thin <- sigma_grand_cal(
        A_1_thin[,, i], A_2_thin[,, i],
        phi_post_1_thin[,, i], phi_post_2_thin[,, i],
        p1, p2
      )
      
      CC_thin <- CC_function_grand_2(Phi_grand_MCMC_thin, p1, p2, CCA_select)
      
      CCA_thin[,, i] <- CC_thin$CCA_grand
      Direction_CCA_Vec1_thin[,, i] <- CC_thin$Direction_CCA_Vec1_grand
      Direction_CCA_Vec2_thin[,, i] <- CC_thin$Direction_CCA_Vec2_grand
      
      log_det_MCMC_thin[, i] <-
        determinant(Phi_grand_MCMC_thin, logarithm = TRUE)$modulus[1]
      
      Log_likelihood_MCMC_Grand_thin[, i] <-
        lk_para_calculator2(
          X_grand,
          mu_1 = as.matrix(mu_1_thin[,, i]),
          mu_2 = as.matrix(mu_2_thin[,, i]),
          phi_grand = Phi_grand_MCMC_thin
        )
      
      i <- i + 1
    }
    
    if (iter == iter_max) {
      print(paste0("Iteration Number", iter))
    }
  }
  
  ##########################################################
  ### Return thinned posterior draws + diagnostics
  ##########################################################
  return(list(
    A_1_MCMC = A_1_thin,
    A_2_MCMC = A_2_thin,
    Mu_1_MCMC = mu_1_thin,
    Mu_2_MCMC = mu_2_thin,
    Omega_1_MCMC = Omega_post_1_thin,
    Omega_2_MCMC = Omega_post_2_thin,
    log_det_MCMC = log_det_MCMC_thin,
    CCA_MCMC = CCA_thin,
    Direction_CCA_Vec1_MCMC = Direction_CCA_Vec1_thin,
    Direction_CCA_Vec2_MCMC = Direction_CCA_Vec2_thin,
    Log_likelihood_MCMC_Grand = Log_likelihood_MCMC_Grand_thin,
    tausqpost_1 = tausqpost_1,
    tausqpost_2 = tausqpost_2,
    eta = eta
  ))
}
############################################################
