############################################################
##### Parameters (example usage — commented out)
############################################################

# load("Breast_Cancer_Cleaned.Rdata")
# X_1 <- Data_Set_1$X_1_dna     # DNA copy number matrix
# X_2 <- Data_Set_1$X_2_rna     # RNA expression matrix
# p1  <- ncol(X_1)              # number of DNA features
# p2  <- ncol(X_2)              # number of RNA features
# no_patients <- nrow(X_1)

############################################################
### =========================
### Libraries
### =========================
############################################################

library(GFA)     # Group Factor Analysis (Bayesian multi-view factor model)
library(MASS)    # Linear algebra utilities
library(coda)    # MCMC diagnostics (ESS, etc.)

############################################################
### =========================
### Helper Functions
### =========================
############################################################

### Euclidean norm (used for direction normalization)
norm_vec <- function(x) sqrt(sum(x^2))

############################################################
### Inverse-Gamma sampler
### Used for Bayesian variance sampling if needed
############################################################
rinvgamma <- function(n, shape, rate = 1, scale = 1 / rate) {
  if (missing(rate) && !missing(scale)) rate <- 1 / scale
  1 / rgamma(n, shape, rate)
}

############################################################
### Build grand covariance matrix from GFA components
### Constructs full block covariance:
###
### [ A1 A1' + phi1      A1 A2' ]
### [ A2 A1'              A2 A2' + phi2 ]
############################################################
sigma_grand_cal <- function(A_1, A_2, phi_1, phi_2, p1, p2) {
  
  # Within-view covariances
  m_11 <- A_1 %*% t(A_1) + phi_1
  m_22 <- A_2 %*% t(A_2) + phi_2
  
  # Cross-view covariances
  m_12 <- A_1 %*% t(A_2)
  m_21 <- A_2 %*% t(A_1)
  
  # Assemble full covariance
  sigma_grand <- matrix(0, nrow = p1 + p2, ncol = p1 + p2)
  
  sigma_grand[1:p1, 1:p1] <- m_11
  sigma_grand[1:p1, (p1 + 1):(p1 + p2)] <- m_12
  sigma_grand[(p1 + 1):(p1 + p2), 1:p1] <- m_21
  sigma_grand[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- m_22
  
  sigma_grand
}

############################################################
### =========================
### CCA from Model-Implied Covariance
### =========================
###
### Computes canonical correlations and directions from
### the covariance implied by GFA posterior draws.
###
### Returns:
### - First CCA_select canonical correlations
### - Corresponding canonical direction vectors
############################################################
CC_function_GFA <- function(X_1, X_2,
                            A_1, A_2,
                            phi_1, phi_2,
                            CCA_select) {
  
  ##########################################################
  ### Construct model-implied covariance blocks
  ##########################################################
  m_11 <- A_1 %*% t(A_1) + phi_1
  m_12 <- A_1 %*% t(A_2)
  m_21 <- A_2 %*% t(A_1)
  m_22 <- A_2 %*% t(A_2) + phi_2
  
  ##########################################################
  ### Compute inverse square roots of covariance matrices
  ### via eigendecomposition
  ##########################################################
  eigen_x <- eigen(m_11)
  eigen_y <- eigen(m_22)
  
  sqrm_11 <- eigen_x$vectors %*%
    diag(eigen_x$values^(-0.5)) %*%
    solve(eigen_x$vectors)
  
  sqrm_22 <- eigen_y$vectors %*%
    diag(eigen_y$values^(-0.5)) %*%
    solve(eigen_y$vectors)
  
  ##########################################################
  ### Construct CCA eigenvalue problems
  ##########################################################
  m1 <- sqrm_11 %*% m_12 %*% solve(m_22) %*%
    t(m_12) %*% sqrm_11
  
  m2 <- sqrm_22 %*% m_21 %*% solve(m_11) %*%
    t(m_21) %*% sqrm_22
  
  m1e <- eigen(m1)
  m2e <- eigen(m2)
  
  ##########################################################
  ### Compute canonical correlations
  ##########################################################
  J <- min(nrow(m1e$vectors),
           nrow(m2e$vectors),
           CCA_select)
  
  Canonical_Correaltions <- rep(NA_real_, CCA_select)
  
  for (j in 1:J) {
    Canonical_Correaltions[j] <-
      Re(m1e$vectors[, j]) %*%
      sqrm_11 %*% m_12 %*%
      sqrm_22 %*%
      Re(m2e$vectors[, j])
  }
  
  ##########################################################
  ### Sign correction
  ### CCA directions are sign-invariant,
  ### so flip signs to ensure positive correlations
  ##########################################################
  for (j in 1:J) {
    if (!is.na(Canonical_Correaltions[j]) &&
        sign(Canonical_Correaltions[j]) == -1) {
      
      Canonical_Correaltions[j] <- abs(Canonical_Correaltions[j])
      m2e$vectors[, j] <- -1 * m2e$vectors[, j]
    }
  }
  
  return(list(
    CCA_grand = Canonical_Correaltions[1:CCA_select],
    Direction_CCA_Vec1_grand =
      m1e$vectors[, 1:CCA_select, drop = FALSE],
    Direction_CCA_Vec2_grand =
      m2e$vectors[, 1:CCA_select, drop = FALSE]
  ))
}

############################################################
### =========================
### GFA Wrapper for Posterior CCA
### =========================
###
### Fits GFA model and extracts posterior draws of:
### - Canonical correlations
### - Canonical directions
###
### Returns full MCMC samples
############################################################
GFA_function <- function(d,
                         X_1,
                         X_2,
                         burn_iter,
                         mcmc_iter,
                         p1 = NULL,
                         p2 = NULL,
                         iter_saved,
                         CCA_select) {
  
  set.seed(666)
  
  ##########################################################
  ### Ensure inputs are matrices
  ##########################################################
  if (!is.matrix(X_1)) X_1 <- as.matrix(X_1)
  if (!is.matrix(X_2)) X_2 <- as.matrix(X_2)
  
  p1 <- ncol(X_1)
  p2 <- ncol(X_2)
  
  ##########################################################
  ### Allocate storage for posterior samples
  ##########################################################
  CCA_GFA <- array(NA_real_,
                   dim = c(1, CCA_select, iter_saved))
  
  Direction_CCA_vec1_GFA <-
    array(NA_real_,
          dim = c(p1, CCA_select, iter_saved))
  
  Direction_CCA_vec2_GFA <-
    array(NA_real_,
          dim = c(p2, CCA_select, iter_saved))
  
  ##########################################################
  ### Configure GFA MCMC options
  ##########################################################
  default_values <- getDefaultOpts()
  default_values$iter.max    <- mcmc_iter
  default_values$iter.saved  <- iter_saved
  default_values$iter.burnin <- burn_iter
  
  ##########################################################
  ### Fit Bayesian Group Factor Analysis model
  ##########################################################
  GFA_analysis <- gfa(Y = list(X_1, X_2),
                      K = 50,
                      opts = default_values)
  
  K_fit <- GFA_analysis$K
  
  ##########################################################
  ### Extract posterior draws and compute CCA
  ##########################################################
  for (iter in 1:iter_saved) {
    
    # Extract loading matrices for both views
    A1 <- GFA_analysis$posterior$W[iter, 1:p1, , drop = FALSE]
    A2 <- GFA_analysis$posterior$W[iter,
                                   (p1 + 1):(p1 + p2),
                                   , drop = FALSE]
    
    A1 <- matrix(A1, nrow = p1)
    A2 <- matrix(A2, nrow = p2)
    
    # Residual precision parameters
    tau1 <- GFA_analysis$posterior$tau[iter, 1:p1]
    tau2 <- GFA_analysis$posterior$tau[iter, (p1 + 1):(p1 + p2)]
    
    phi1 <- diag(1 / tau1, p1)
    phi2 <- diag(1 / tau2, p2)
    
    # Compute CCA from model-implied covariance
    CCA_GFA_2 <- CC_function_GFA(
      X_1 = X_1,
      X_2 = X_2,
      A_1 = A1,
      A_2 = A2,
      phi_1 = phi1,
      phi_2 = phi2,
      CCA_select = CCA_select
    )
    
    # Store canonical correlations
    cca_vec <- as.numeric(CCA_GFA_2$CCA_grand)
    CCA_GFA[1, , iter] <- abs(cca_vec[1:CCA_select])
    
    # Store direction vectors
    Direction_CCA_vec1_GFA[, , iter] <-
      Re(CCA_GFA_2$Direction_CCA_Vec1_grand)
    
    Direction_CCA_vec2_GFA[, , iter] <-
      Re(CCA_GFA_2$Direction_CCA_Vec2_grand)
  }
  
  ##########################################################
  ### Return posterior samples
  ##########################################################
  return(list(
    CCA_MCMC = CCA_GFA,
    Direction_CCA_Vec1_MCMC = Direction_CCA_vec1_GFA,
    Direction_CCA_Vec2_MCMC = Direction_CCA_vec2_GFA,
    K_fit = K_fit
  ))
}





























