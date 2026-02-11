

##### Required Library #####
# coda is used for MCMC diagnostics (effective sample size)
library(coda)

##### Utility Functions #####

# Compute Euclidean norm of a vector
# Used to normalize canonical direction vectors
norm_vec <- function(x) sqrt(sum(x^2))


############################################################
# Function: CC_function_grand_2
# Purpose:
#   Computes canonical correlations and direction vectors
#   from a full covariance matrix (grand covariance).
# Inputs:
#   Covariance_matrix_Grand : (p1+p2) x (p1+p2) covariance matrix
#   p1, p2                  : dimensions of view 1 and view 2
#   CCA_select              : number of canonical correlations to compute
# Output:
#   List containing canonical correlations and direction vectors
############################################################

CC_function_grand_2 <- function(Covariance_matrix_Grand, p1, p2, CCA_select)
{
  sigma_grand <- Covariance_matrix_Grand
  
  # Partition covariance matrix into block components
  m_11 <- sigma_grand[1:p1, 1:p1]          # Variance of View 1
  m_12 <- sigma_grand[1:p1, p1+1:p2]       # Cross-covariance (1 -> 2)
  m_21 <- sigma_grand[p1+1:p2, 1:p1]       # Cross-covariance (2 -> 1)
  m_22 <- sigma_grand[p1+1:p2, p1+1:p2]    # Variance of View 2
  
  # Compute inverse square roots using eigen-decomposition
  eigen_x <- eigen(m_11)
  eigen_y <- eigen(m_22)
  
  sqrm_11 <- eigen_x$vectors %*%
    diag(eigen_x$values^(-0.5)) %*%
    solve(eigen_x$vectors)
  
  sqrm_22 <- eigen_y$vectors %*%
    diag(eigen_y$values^(-0.5)) %*%
    solve(eigen_y$vectors)
  
  # Construct CCA eigenvalue problems
  m1 <- sqrm_11 %*% m_12 %*% solve(m_22) %*% t(m_12) %*% sqrm_11
  m2 <- sqrm_22 %*% m_21 %*% solve(m_11) %*% t(m_21) %*% sqrm_22
  
  m1e <- eigen(m1)
  m2e <- eigen(m2)
  
  # Compute canonical correlations
  Canonical_Correaltions <- numeric(CCA_select)
  
  for(j in 1:CCA_select) {
    Canonical_Correaltions[j] <-
      Re(m1e$vectors[, j]) %*%
      sqrm_11 %*%
      m_12 %*%
      sqrm_22 %*%
      Re(m2e$vectors[, j])
  }
  
  # Ensure canonical correlations are positive (sign identifiability)
  for(i in 1:length(Canonical_Correaltions)) {
    if(sign(Canonical_Correaltions[i]) == -1) {
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
# Function: sigma_grand_cal
# Purpose:
#   Constructs full covariance matrix from factor loadings
#   and view-specific residual covariance matrices.
############################################################

sigma_grand_cal <- function(A_1, A_2, phi_1, phi_2, p1, p2)
{
  # Block components of covariance
  m_11 <- A_1 %*% t(A_1) + phi_1
  m_12 <- A_1 %*% t(A_2)
  m_21 <- A_2 %*% t(A_1)
  m_22 <- A_2 %*% t(A_2) + phi_2
  
  # Assemble grand covariance matrix
  sigma_grand <- matrix(0, nrow = p1 + p2, ncol = p1 + p2)
  
  sigma_grand[1:p1, 1:p1] <- m_11
  sigma_grand[1:p1, (p1+1):(p1+p2)] <- m_12
  sigma_grand[(p1+1):(p1+p2), 1:p1] <- m_21
  sigma_grand[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- m_22
  
  return(sigma_grand)
}


############################################################
# Function: HDI_function
# Purpose:
#   Computes highest density interval (HDI approximation via quantiles)
#   and checks whether a true value lies inside the interval.
############################################################

HDI_function <- function(col_data, true_val, alpha)
{
  col_data <- Re(col_data)
  true_val <- Re(true_val)
  
  # Equal-tailed credible interval
  ci_hdi <- quantile(col_data, c(alpha/2, 1 - alpha/2))
  
  # Indicator: does interval contain true_val?
  Ind_CCA <- ifelse(ci_hdi[[1]] < true_val & true_val < ci_hdi[[2]], 1, 0)
  
  return(list(HDI = ci_hdi, Count = Ind_CCA))
}


############################################################
# Function: Analysis_Summary
# Purpose:
#   Post-processes MHCFM MCMC output.
#   Computes posterior summaries, credible intervals,
#   effective sample sizes, and shrinkage diagnostics.
############################################################

Analysis_Summary <- function(Output,
                             Data,
                             shrinkage_threshold = 0.2,
                             plot_indicator,
                             plot_components = 10,
                             alpha.CC,
                             alpha.dir,
                             shrinkage_check = TRUE)
{
  
  # Extract MCMC dimensions
  thin_seq <- 1:length(Output$CCA_MCMC[,1,])
  
  p1 <- nrow(Output$Direction_CCA_Vec1_MCMC[,1,])
  p2 <- nrow(Output$Direction_CCA_Vec2_MCMC[,1,])
  CCA_select <- ncol(Output$CCA_MCMC)
  
  iter_saved <- length(Output$CCA_MCMC[,1,])
  
  ##########################################################
  # Store posterior samples
  ##########################################################
  
  CCA_thin <- Output$CCA_MCMC
  Direction_CCA_Vec1_thin <- Output$Direction_CCA_Vec1_MCMC
  Direction_CCA_Vec2_thin <- Output$Direction_CCA_Vec2_MCMC
  
  ##########################################################
  # Sign correction for identifiability
  ##########################################################
  
  Direction_vector1_SC <- array(dim = c(p1, CCA_select, iter_saved))
  Direction_vector2_SC <- array(dim = c(p2, CCA_select, iter_saved))
  
  dir_vec <- 1  # Only first canonical vector adjusted
  
  for(iter in 1:iter_saved)
  {
    vec_1 <- Direction_CCA_Vec1_thin[, dir_vec, iter]
    
    if(sign(vec_1[which.max(abs(vec_1))]) == -1)
    {
      Direction_vector1_SC[, dir_vec, iter] <-
        -Direction_CCA_Vec1_thin[, dir_vec, iter]
      
      Direction_vector2_SC[, dir_vec, iter] <-
        -Direction_CCA_Vec2_thin[, dir_vec, iter]
    }
    else
    {
      Direction_vector1_SC[, dir_vec, iter] <-
        Direction_CCA_Vec1_thin[, dir_vec, iter]
      
      Direction_vector2_SC[, dir_vec, iter] <-
        Direction_CCA_Vec2_thin[, dir_vec, iter]
    }
  }
  
  ##########################################################
  # Posterior Means
  ##########################################################
  
  Direction_CCA_Vec1_mean <-
    apply(Direction_vector1_SC, c(1,2), mean)
  
  Direction_CCA_Vec2_mean <-
    apply(Direction_vector2_SC, c(1,2), mean)
  
  ##########################################################
  # Credible intervals for canonical correlations
  ##########################################################
  
  HDI_CCA <- matrix(NA, nrow = CCA_select, ncol = 2)
  
  for(j in 1:CCA_select)
  {
    HDI_cal <- HDI_function(
      col_data = CCA_thin[, j, ],
      true_val = 0,
      alpha = alpha.CC
    )
    
    HDI_CCA[j, ] <- HDI_cal$HDI
  }
  
  ##########################################################
  # Effective sample sizes
  ##########################################################
  
  effective_size_first_CCA <-
    effectiveSize(CCA_thin[, 1, ])
  
  effective_size_log_likelihood <-
    effectiveSize(apply(Output$Log_likelihood_MCMC_Grand, 2, sum))
  
  ##########################################################
  # Shrinkage diagnostic
  ##########################################################
  
  if(shrinkage_check)
  {
    Shrinkage_calculator <-
      mean(CCA_thin[,1,] < shrinkage_threshold)
    
    if(Shrinkage_calculator > 0.5)
      message("Run DFSM Model; Overshrinkage suspected")
    else
      message("Continue with NDFSM; No evidence of overshrinkage")
  }
  else
  {
    Shrinkage_calculator <- NA
  }
  
  ##########################################################
  # Final Output
  ##########################################################
  
  return(list(
    Shrinkage_calculator = Shrinkage_calculator,
    effective_size_first_CCA = effective_size_first_CCA,
    effective_size_log_likelihood = effective_size_log_likelihood,
    HDI_CCA = HDI_CCA
  ))
}





















