##### Parameters 

# load("Breast_Cancer_Cleaned.Rdata")

# ex<-1    
# 
# X_1<-Data_Set_1$X_1_dna
# 
# X_2<-Data_Set_1$X_2_rna
# 
# p1<-ncol(X_1)
# 
# p2<-ncol(X_2)

#no_experiment<-1

#burn_iter<-60000

#mcmc_iter<-540000

# no_patients<-nrow(X_1) 


###
###   LOAD AN R LIBRARY (SIMULATES MULTIVARIATE NORMAL)
###


### =========================
### Libraries
### =========================
library(GFA)
library(MASS)
library(coda)

### =========================
### Helpers
### =========================

# Norm
norm_vec <- function(x) sqrt(sum(x^2))

# Inverse-gamma sampler
rinvgamma <- function(n, shape, rate = 1, scale = 1 / rate) {
  if (missing(rate) && !missing(scale)) rate <- 1 / scale
  1 / rgamma(n, shape, rate)
}

# Build grand covariance from blocks (if you need it elsewhere)
sigma_grand_cal <- function(A_1, A_2, phi_1, phi_2, p1, p2) {
  m_11 <- A_1 %*% t(A_1) + phi_1
  m_12 <- A_1 %*% t(A_2)
  m_21 <- A_2 %*% t(A_1)
  m_22 <- A_2 %*% t(A_2) + phi_2
  
  sigma_grand <- matrix(0, nrow = p1 + p2, ncol = p1 + p2)
  sigma_grand[1:p1, 1:p1] <- m_11
  sigma_grand[1:p1, (p1 + 1):(p1 + p2)] <- m_12
  sigma_grand[(p1 + 1):(p1 + p2), 1:p1] <- m_21
  sigma_grand[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- m_22
  sigma_grand
}

### =========================
### CCA from model-implied covariance (A_1, A_2, phi_1, phi_2)
### Returns exactly CCA_select values + directions
### =========================
CC_function_GFA <- function(X_1, 
                        X_2,
                        A_1, 
                        A_2, 
                        phi_1, 
                        phi_2,
                        CCA_select) {

  # model-implied cov blocks
  m_11 <- A_1 %*% t(A_1) + phi_1
  m_12 <- A_1 %*% t(A_2)
  m_21 <- A_2 %*% t(A_1)
  m_22 <- A_2 %*% t(A_2) + phi_2

  # eigendecompositions for inverse square roots
  eigen_x <- eigen(m_11)
  eigen_y <- eigen(m_22)

  sqrm_11 <- (eigen_x$vectors) %*% diag(eigen_x$values^(-0.5)) %*% solve(eigen_x$vectors)
  sqrm_22 <- (eigen_y$vectors) %*% diag(eigen_y$values^(-0.5)) %*% solve(eigen_y$vectors)

  # CCA matrices
  m1 <- sqrm_11 %*% m_12 %*% solve(m_22) %*% t(m_12) %*% sqrm_11
  m2 <- sqrm_22 %*% m_21 %*% solve(m_11) %*% t(m_21) %*% sqrm_22

  m1e <- eigen(m1)
  m2e <- eigen(m2)

  # how many CCs can we compute?
  J <- min(nrow(m1e$vectors), nrow(m2e$vectors), CCA_select)

  Canonical_Correaltions <- rep(NA_real_, CCA_select)

  for (j in 1:J) {
    Canonical_Correaltions[j] <-
      Re(m1e$vectors[, j]) %*% sqrm_11 %*% m_12 %*% sqrm_22 %*% Re(m2e$vectors[, j])
  }

  # sign correction + flip corresponding y-direction if needed
  for (j in 1:J) {
    if (!is.na(Canonical_Correaltions[j]) && sign(Canonical_Correaltions[j]) == -1) {
      Canonical_Correaltions[j] <- abs(Canonical_Correaltions[j])
      m2e$vectors[, j] <- -1 * m2e$vectors[, j]
    }
  }

  return(list(
    CCA_grand = Canonical_Correaltions[1:CCA_select],
    Direction_CCA_Vec1_grand = m1e$vectors[, 1:CCA_select, drop = FALSE],
    Direction_CCA_Vec2_grand = m2e$vectors[, 1:CCA_select, drop = FALSE]
  ))
}

### =========================
### GFA wrapper that returns CCA posterior draws
### FIXED: correct indexing + correct list element names
### =========================
GFA_function <- function(d, X_1, X_2, burn_iter, mcmc_iter, p1 = NULL, p2 = NULL,
                         iter_saved, CCA_select) {
  
  set.seed(666)
  
  # enforce matrices
  if (!is.matrix(X_1)) X_1 <- as.matrix(X_1)
  if (!is.matrix(X_2)) X_2 <- as.matrix(X_2)
  
  p1 <- ncol(X_1)
  p2 <- ncol(X_2)
  
  # storage
  CCA_GFA <- array(NA_real_, dim = c(1, CCA_select, iter_saved))
  Direction_CCA_vec1_GFA <- array(NA_real_, dim = c(p1, CCA_select, iter_saved))
  Direction_CCA_vec2_GFA <- array(NA_real_, dim = c(p2, CCA_select, iter_saved))
  
  # GFA options
  default_values <- getDefaultOpts()
  default_values$iter.max <- mcmc_iter
  default_values$iter.saved <- iter_saved
  default_values$iter.burnin <- burn_iter
  
  # Fit GFA (K is fixed at 50 in your original code; keep it)
  GFA_analysis <- gfa(Y = list(X_1, X_2), K = 50, opts = default_values)
  
  K_fit <- GFA_analysis$K
  
  # posterior arrays: W[iter, feature, k], tau[iter, feature]
  for (iter in 1:iter_saved) {
    
    # Correct slicing:
    # - features are 1:(p1+p2)
    # - view1 is 1:p1
    # - view2 is (p1+1):(p1+p2)
    A1 <- GFA_analysis$posterior$W[iter, 1:p1, , drop = FALSE]
    A2 <- GFA_analysis$posterior$W[iter, (p1 + 1):(p1 + p2), , drop = FALSE]
    
    # drop iter dimension and ensure matrices
    A1 <- matrix(A1, nrow = p1)
    A2 <- matrix(A2, nrow = p2)
    
    tau1 <- GFA_analysis$posterior$tau[iter, 1:p1]
    tau2 <- GFA_analysis$posterior$tau[iter, (p1 + 1):(p1 + p2)]
    
    phi1 <- diag(1 / tau1, p1)
    phi2 <- diag(1 / tau2, p2)
    
    CCA_GFA_2 <- CC_function_GFA(
      X_1 = X_1,
      X_2 = X_2,
      A_1 = A1,
      A_2 = A2,
      phi_1 = phi1,
      phi_2 = phi2,
      CCA_select = CCA_select
    )
    
    # Ensure exactly CCA_select values (pad with NA if needed)
    cca_vec <- as.numeric(CCA_GFA_2$CCA_grand)
    if (length(cca_vec) < CCA_select) {
      cca_vec <- c(cca_vec, rep(NA_real_, CCA_select - length(cca_vec)))
    } else if (length(cca_vec) > CCA_select) {
      cca_vec <- cca_vec[1:CCA_select]
    }
    
    CCA_GFA[1, , iter] <- abs(cca_vec)
    
    # Directions (take real part; force correct dims)
    Direction_CCA_vec1_GFA[, , iter] <- Re(CCA_GFA_2$Direction_CCA_Vec1_grand)[, 1:CCA_select, drop = FALSE]
    Direction_CCA_vec2_GFA[, , iter] <- Re(CCA_GFA_2$Direction_CCA_Vec2_grand)[, 1:CCA_select, drop = FALSE]
  }
  
  return(list(
    CCA_MCMC = CCA_GFA,
    Direction_CCA_Vec1_MCMC = Direction_CCA_vec1_GFA,
    Direction_CCA_Vec2_MCMC = Direction_CCA_vec2_GFA,
    K_fit = K_fit
  ))
}




#K_fit_ex<-array(dim=c(1,no_experiment))

#############################################################################

# ## CCA by GFA Function#####
# 
# GFA_Output<-GFA_function(X_1=X_1,X_2=X_2,p1,p2,d_actual,iter_saved=iter_saved,CCA_select)
# 
# #### saving the files ####
# 
# filename_GFA_output<-paste0(20,"GFA_BC_Output",sep="")
# 
# save(GFA_Output,file=paste(filename_GFA_output,".Rdata",sep=""))































