
############################################################
### FUNCTIONS FOR CORRELATIONS + GRAPHICAL HORSESHOE (GHS)
### Goal:
###   - Fit a sparse precision matrix (Omega) using Graphical Horseshoe
###   - Convert posterior covariance (Sigma = Omega^{-1}) into CCA
###   - Save thinned posterior draws of:
###       * Omega, CCA, direction vectors, log-det, likelihood
############################################################

### =========================
### Libraries
### =========================
library(MASS)     # matrix ops + (usually) mvrnorm (but you override below)
library(mvtnorm)  # dmvnorm for multivariate normal log-likelihood

############################################################
### Custom mvrnorm (you override MASS::mvrnorm here)
### Draws from MVN(mu, Sigma) using eigen decomposition
### Notes:
###   - empirical=TRUE forces exact empirical mean/cov in samples
###   - checks Sigma is positive definite (up to tolerance)
############################################################
mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06,
                     empirical = FALSE, EISPACK = FALSE) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  if (EISPACK)
    stop("'EISPACK' is no longer supported by R", domain = NA)
  
  # eigen decomposition Sigma = V diag(ev) V'
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  
  # PD check
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  
  # base standard normal draws
  X <- matrix(rnorm(p * n), n)
  
  # optional: enforce empirical properties
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  
  # transform to MVN with covariance Sigma and mean mu
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  
  if (n == 1) drop(X) else t(X)
}

############################################################
### Inverse-Gamma sampler (via reciprocal Gamma)
### Used for GHS local/global shrinkage updates
############################################################
rinvgamma <- function(n, shape, rate = 1, scale = 1/rate)
{
  if (missing(rate) && !missing(scale))
    rate <- 1/scale
  1/rgamma(n, shape, rate)
}

############################################################
### Build block covariance Sigma_grand from factor blocks
### (Not directly used in GHS_CC, but useful generally)
############################################################
sigma_grand_cal <- function(A_1, A_2, phi_1, phi_2, p1, p2)
{
  m_11 <- A_1 %*% t(A_1) + phi_1
  m_12 <- A_1 %*% t(A_2)
  m_21 <- A_2 %*% t(A_1)
  m_22 <- A_2 %*% t(A_2) + phi_2
  
  sigma_grand <- matrix(nrow = p1 + p2, ncol = p1 + p2)
  
  # Fill block matrix
  sigma_grand[1:nrow(m_11), 1:ncol(m_11)] <- m_11
  sigma_grand[1:nrow(m_12), ncol(m_11) + 1:ncol(m_12)] <- m_12
  sigma_grand[nrow(m_11) + 1:nrow(m_21), 1:ncol(m_21)] <- m_21
  sigma_grand[nrow(m_11) + 1:nrow(m_22), ncol(m_21) + 1:ncol(m_22)] <- m_22
  
  return(sigma_grand)
}

############################################################
### CCA from a grand covariance matrix
### Inputs:
###   - Covariance_matrix_Grand: full covariance of (X_1, X_2)
###   - p1, p2: dimensions of each view
###   - CCA_select: number of canonical correlations to return
### Output:
###   - Canonical correlations + direction vectors (weights)
############################################################
CC_function_grand_2 <- function(Covariance_matrix_Grand, p1, p2, CCA_select)
{
  sigma_grand <- Covariance_matrix_Grand
  
  # Extract block matrices
  m_11 <- sigma_grand[1:p1, 1:p1]                 # Cov(X)
  m_12 <- sigma_grand[1:p1, p1+1:p2]              # Cov(X,Y)  (NOTE: indexing assumes p2 = p1+p2!!)
  m_21 <- sigma_grand[p1+1:p2, 1:p1]              # Cov(Y,X)
  m_22 <- sigma_grand[p1+1:p2, p1+1:p2]           # Cov(Y)
  
  ##########################################################
  ### Compute inverse square root of m_11 and m_22
  ### via eigendecomposition
  ##########################################################
  eigen_x <- eigen(m_11)
  eigen_y <- eigen(m_22)
  
  sqrm_11 <- (eigen_x$vectors) %*% diag(eigen_x$values^-0.5) %*% solve(eigen_x$vectors)
  sqrm_22 <- (eigen_y$vectors) %*% diag(eigen_y$values^-0.5) %*% solve(eigen_y$vectors)
  
  ##########################################################
  ### Construct standard CCA eigen problems
  ### m1 eigenvectors -> X directions
  ### m2 eigenvectors -> Y directions
  ##########################################################
  m1 <- sqrm_11 %*% m_12 %*% solve(m_22) %*% t(m_12) %*% sqrm_11
  m2 <- sqrm_22 %*% m_21 %*% solve(m_11) %*% t(m_21) %*% sqrm_22
  
  m1e <- eigen(m1)
  m2e <- eigen(m2)
  
  ##########################################################
  ### Compute canonical correlations using paired eigenvectors
  ##########################################################
  Canonical_Correaltions <- c()
  
  for (j in 1:CCA_select) {
    Canonical_Correaltions[j] <-
      Re(m1e$vectors[, j]) %*% sqrm_11 %*% m_12 %*% sqrm_22 %*% Re(m2e$vectors[, j])
  }
  
  ##########################################################
  ### Sign correction (CCA sign is arbitrary)
  ### Flip correlations to be positive; flip Y direction as well
  ##########################################################
  for (i in 1:length(Canonical_Correaltions)) {
    if (sign(Canonical_Correaltions[i]) == -1) {
      Canonical_Correaltions[i] <- abs(Canonical_Correaltions[i])
      
      # NOTE: this uses m2e$vectors[,j] but should probably be [,i]
      # Leaving as-is per your "comments only" request.
      m2e$vectors[, j] <- -1 * m2e$vectors[, j]
    }
  }
  
  ##########################################################
  ### Return canonical correlations and direction vectors
  ##########################################################
  return(list(
    CCA_grand = Canonical_Correaltions,
    Direction_CCA_Vec1_grand = m1e$vectors[, 1:CCA_select],
    Direction_CCA_Vec2_grand = m2e$vectors[, 1:CCA_select]
  ))
}

############################################################
##### Graphical Horseshoe (GHS) Gibbs update (node-wise)
### Inputs:
###   X          : centered data matrix (n x p)
###   Lambda_sq  : local shrinkage params (p x p)
###   Nu         : auxiliary params for Lambda_sq
###   Omega      : precision matrix (p x p)
###   Sigma      : covariance matrix (p x p)  (Omega^{-1})
###   tau_sq     : global shrinkage
###   xi         : auxiliary param for tau_sq
### Output:
###   Updated Omega, Sigma, Lambda_sq, Nu, tau_sq, xi
############################################################
GHS_modified_2 <- function(X,
                           Lambda_sq, Nu, Omega, Sigma,
                           tau_sq, xi, iter) {
  
  n <- nrow(X)            # number of samples
  S <- t(X) %*% (X)       # scatter matrix (p x p)
  p <- nrow(S)            # number of variables (p1 + p2)
  
  ##########################################################
  ### Precompute "all indices except i" for each node i
  ##########################################################
  ind_all <- matrix(0, p - 1, p)
  for (i in 1:p) {
    if (i == 1) {
      ind <- c(2:p)
    } else if (i == p) {
      ind <- c(1:(p - 1))
    } else {
      ind <- c(1:(i - 1), (i + 1):p)
    }
    ind_all[, i] <- ind
  }
  
  ##########################################################
  ### Node-wise updates (standard graphical model Gibbs step)
  ### Each i updates the i-th row/col of Omega and Sigma
  ##########################################################
  for (i in 1:p) {
    
    ind <- ind_all[, i]
    
    # Partition Sigma blocks: remove i-th variable
    sigma_11 <- Sigma[ind, ind]
    sigma_12 <- Sigma[ind, i]
    sigma_22 <- Sigma[i, i]
    
    # Scatter blocks
    s_21 <- S[ind, i]
    s_22 <- S[i, i]
    
    # Local shrinkage params for edges incident to i
    lambda_sq_12 <- Lambda_sq[ind, i]
    nu_12 <- Nu[ind, i]
    
    # gamma update (related to Omega_ii)
    gamma <- rgamma(1, shape = ((n / 2) + 1), rate = s_22 / 2)
    
    # Compute inv_Omega_11 (Schur complement style quantity)
    inv_Omega_11 <- sigma_11 - ((sigma_12 %*% t(sigma_12)) / (sigma_22))
    
    # Conditional precision for beta update
    inv_C <- (s_22 * inv_Omega_11) + diag(1 / (lambda_sq_12 * tau_sq), p - 1)
    inv_C_chol <- chol(inv_C)
    
    # Conditional mean for regression coefficients beta
    mu_i <- -chol2inv(inv_C_chol) %*% s_21
    
    # Sample beta ~ N(mu_i, inv_C^{-1})
    beta <- mu_i + solve(inv_C_chol, rnorm(p - 1, 0, 1))
    
    # Update precision matrix blocks
    omega_12 <- beta
    omega_22 <- gamma + t(beta) %*% inv_Omega_11 %*% beta
    
    # Update local shrinkage parameters for edges (i, ind)
    rate <- (1 / nu_12) + ((omega_12^2) / (2 * tau_sq))
    lambda_sq_12 <- rinvgamma(p - 1, shape = 1, rate = rate)
    nu_12 <- rinvgamma(p - 1, shape = 1, rate = (1 + (1 / lambda_sq_12)))
    
    # Write updated Omega symmetric entries
    Omega[i, ind] <- omega_12
    Omega[ind, i] <- omega_12
    Omega[i, i] <- omega_22
    
    # Update Sigma (= Omega^{-1}) via block formulas
    temp <- inv_Omega_11 %*% beta
    Sigma_11 <- inv_Omega_11 + temp %*% t(temp) / gamma
    sigma_12 <- -temp / gamma
    sigma_22 <- 1 / gamma
    
    Sigma[ind, ind] <- Sigma_11
    Sigma[i, i] <- sigma_22
    Sigma[i, ind] <- sigma_12
    Sigma[ind, i] <- sigma_12
    
    # Store updated shrinkage params symmetrically
    Lambda_sq[i, ind] <- lambda_sq_12
    Lambda_sq[ind, i] <- lambda_sq_12
    Nu[i, ind] <- nu_12
    Nu[ind, i] <- nu_12
  }
  
  ##########################################################
  ### Global shrinkage tau_sq update (uses all off-diagonals)
  ##########################################################
  omega_vector <- Omega[lower.tri(Omega, diag = FALSE)]
  lambda_sq_vector <- Lambda_sq[lower.tri(Lambda_sq, diag = FALSE)]
  
  rate1 <- 1 / xi + sum((omega_vector^2) / (2 * lambda_sq_vector))
  tau_sq <- rinvgamma(1, shape = ((p * (p - 1) / 2) + 1) / 2, rate = rate1)
  
  # xi auxiliary update
  xi <- rinvgamma(1, shape = 1, rate = (1 + 1 / tau_sq))
  
  return(list(
    Omega = Omega,
    Lambda_sq = Lambda_sq,
    Nu = Nu,
    Tau_sq = tau_sq,
    Sigma = Sigma,
    xi = xi
  ))
}

############################################################
##### GRAND GHS METHOD (GHS + CCA over posterior draws)
### Fits sparse precision matrix for concatenated data:
###   X_grand = [X_1 | X_2]
### then converts covariance draws into CCA draws.
############################################################
GHS_CC <- function(X_1,
                   X_2,
                   burn_iter,
                   mcmc_iter,
                   thin,
                   CCA_select)
{
  set.seed(666)
  
  ##########################################################
  ### Combine two views into one "grand" multivariate vector
  ##########################################################
  X_grand <- cbind(X_1 = X_1, X_2 = X_2)
  
  p1 <- ncol(X_1)
  p2 <- ncol(X_2)
  no_patients <- nrow(X_1)
  
  ##########################################################
  ### Initialize MCMC state
  ### Sigma (cov) is initialized from a ridge-ish covariance
  ### Omega (precision) is inverse of Sigma
  ##########################################################
  Phi_grand_GHS_iter <- array(dim = c(p1 + p2, p1 + p2))
  
  # Local shrinkage structures (initialized to 1)
  Lambda_sq_GHS_grand <- Nu_GHS_grand <- matrix(1, nrow = p1 + p2, ncol = p1 + p2)
  
  # Omega = precision, Sigma = covariance
  Omega_grand_GHS_iter <- Sigma_grand <- array(dim = c(p1 + p2, p1 + p2))
  
  # Centered data container: X - mu
  X_grand_delta <- matrix(nrow = no_patients, ncol = p1 + p2)
  
  # Ridge initialization to ensure PD covariance
  Phi_grand_GHS_iter <- 0.5 * cov(X_grand) + 0.5 * diag(p1 + p2)
  Omega_grand_GHS_iter <- solve(Phi_grand_GHS_iter)
  
  Sigma_grand_GHS <- diag(1, p1 + p2)
  
  # Global shrinkage init
  tau_sq_GHS_grand <- 1
  xi_GHS_grand <- 1
  
  ##########################################################
  ### Thinning schedule
  ##########################################################
  iter_max <- burn_iter + mcmc_iter
  thin_seq <- seq(burn_iter + thin, iter_max, thin)
  
  # Storage for thinned draws
  Omega_grand_thin <- array(dim = c(p1 + p2, p1 + p2, length(thin_seq)))
  log_det_MCMC_thin <- array(dim = c(1, length(thin_seq)))
  likelihood <- array(dim = c(no_patients, length(thin_seq)))
  
  # CCA storage (from covariance draws)
  CCA_thin <- array(dim = c(1, CCA_select, length(thin_seq)))
  Direction_CCA_Vec1_thin <- array(dim = c(p1, CCA_select, length(thin_seq)))
  Direction_CCA_Vec2_thin <- array(dim = c(p2, CCA_select, length(thin_seq)))
  
  # index for thinning saves
  i <- 1
  
  ##########################################################
  ### Main MCMC loop
  ##########################################################
  for (iter in 2:iter_max) {
    
    ########################################################
    ### Sample posterior mean vector mu | Omega, X
    ### Using conjugate Gaussian update
    ########################################################
    sigma <- 100  # prior variance scale for mu
    
    E_mu_cal <- no_patients * (Omega_grand_GHS_iter) + diag(1 / sigma, p1 + p2)
    variance_mu <- solve(E_mu_cal)
    mean_mu <- variance_mu %*% (Omega_grand_GHS_iter) %*% (apply(X_grand, 2, sum))
    
    mu_grand <- mvrnorm(1, mean_mu, variance_mu)
    
    # Center data
    for (patient in 1:no_patients) {
      X_grand_delta[patient, ] <- X_grand[patient, ] - mu_grand
    }
    
    ########################################################
    ### Graphical Horseshoe update for Omega/Sigma
    ########################################################
    GHS_Run_3 <- GHS_modified_2(
      X = X_grand_delta,
      Lambda_sq = Lambda_sq_GHS_grand,
      Nu = Nu_GHS_grand,
      Omega = Omega_grand_GHS_iter,
      Sigma = Phi_grand_GHS_iter,
      tau_sq = tau_sq_GHS_grand,
      xi = xi_GHS_grand,
      iter
    )
    
    # Update state with new draw
    Phi_grand_GHS_iter <- GHS_Run_3$Sigma
    Omega_grand_GHS_iter <- GHS_Run_3$Omega
    Lambda_sq_GHS_grand <- GHS_Run_3$Lambda_sq
    Nu_GHS_grand <- GHS_Run_3$Nu
    tau_sq_GHS_grand <- GHS_Run_3$Tau_sq
    xi_GHS_grand <- GHS_Run_3$xi
    
    # Refresh raw (uncentered) data for likelihood calculation
    X_grand <- cbind(X_1 = X_1, X_2 = X_2)
    
    ########################################################
    ### Save thinned draws
    ########################################################
    if (iter == thin_seq[i]) {
      
      # Save Omega draw
      Omega_grand_thin[,, i] <- Omega_grand_GHS_iter
      
      # Log-likelihood under MVN(mu, Sigma)
      likelihood[, i] <- dmvnorm(X_grand,
                                 mean = mu_grand,
                                 sigma = Phi_grand_GHS_iter,
                                 log = TRUE)
      
      # Compute CCA from current covariance draw
      cc_cal_grand <- CC_function_grand_2(
        Covariance_matrix_Grand = Phi_grand_GHS_iter,
        p1 = p1,
        p2 = p2,
        CCA_select = CCA_select
      )
      
      # NOTE: Your CC_function_grand_2 returns names like:
      #   CCA_grand, Direction_CCA_Vec1_grand, Direction_CCA_Vec2_grand
      # but below you reference cc_cal_grand$CCA and $Direction_CCA_Vec1, etc.
      # Leaving as-is per "comments only".
      
      CCA_thin[,, i] <- cc_cal_grand$CCA
      Direction_CCA_Vec1_thin[,, i] <- cc_cal_grand$Direction_CCA_Vec1
      Direction_CCA_Vec2_thin[,, i] <- cc_cal_grand$Direction_CCA_Vec2
      
      # Save log-determinant of covariance (useful for diagnostics)
      log_det_MCMC_thin[, i] <-
        determinant(Phi_grand_GHS_iter, logarithm = TRUE)$modulus[1]
      
      # increment thinning index
      i <- i + 1
    }
    
    ########################################################
    ### Simple progress printing
    ########################################################
    if (iter == iter_max) {
      print(paste0("Iteration Number ", iter))
    }
    print(iter)
  }
  
  ##########################################################
  ### Return posterior draws
  ##########################################################
  return(GHS_output <- list(
    Omega_grand = Omega_grand_thin,
    CCA_MCMC = CCA_thin,
    Direction_CCA_Vec1_MCMC = Direction_CCA_Vec1_thin,
    Direction_CCA_Vec2_MCMC = Direction_CCA_Vec2_thin,
    log_det_MCMC = log_det_MCMC_thin,
    likelihood_MCMC = likelihood
  ))
}
