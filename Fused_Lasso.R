############################################################
### LOAD REQUIRED LIBRARY
### PMA package is used for Penalized CCA (including fused lasso)
############################################################
library(PMA)

############################################################
### MAIN FUNCTION
### Performs Fused Lasso CCA analysis between DNA and RNA views
############################################################
Fused_lasso_analysis <- function(Data_Set_1 = Data_Set_1,
                                 CCA_select = 10,
                                 gene_chr_250 = gene_chr_250)
{
  
  ##########################################################
  ### Extract the two data views
  ##########################################################
  X_1 <- Data_Set_1$X_1_dna   # DNA copy number data
  X_2 <- Data_Set_1$X_2_rna   # RNA expression data
  
  no_patients   <- nrow(X_1)  # Number of samples
  no_experiment <- 1          # Single experiment structure
  
  p1 <- ncol(X_1)  # Dimension of view 1 (DNA)
  p2 <- ncol(X_2)  # Dimension of view 2 (RNA)
  
  gene_chr <- gene_chr_250    # Chromosome indicator vector
  
  ##########################################################
  ### Helper: Euclidean norm
  ##########################################################
  norm_vec <- function(x) sqrt(sum(x^2))
  
  ##########################################################
  ### Function: Perform CCA using Fused Lasso penalty
  ### Computes both ordered and standard versions
  ##########################################################
  CCA_by_fused_lasso <- function(X_1, X_2, K)
  {
    ### Ordered penalty (fused lasso structure)
    CCA_fused_lasso_cal_O <- CCA(X_1, X_2,
                                 typex = "ordered",
                                 typez = "ordered",
                                 K = K)
    
    ### Standard penalty (lasso without ordering)
    CCA_fused_lasso_cal_std <- CCA(X_1, X_2,
                                   typex = "standard",
                                   typez = "standard",
                                   K = K)
    
    return(list(
      CCA_fused_lasso_cal_O_CCA = CCA_fused_lasso_cal_O$cors,
      CCA_Direction_X_1_fused_lasso_O = CCA_fused_lasso_cal_O$u,
      CCA_Direction_X_2_fused_lasso_O = CCA_fused_lasso_cal_O$v,
      CCA_fused_lasso_cal_std_CCA = CCA_fused_lasso_cal_std$cors,
      CCA_Direction_X_1_fused_lasso_std = CCA_fused_lasso_cal_std$u,
      CCA_Direction_X_2_fused_lasso_std = CCA_fused_lasso_cal_std$v
    ))
  }
  
  ##########################################################
  ### Storage arrays for canonical correlations and directions
  ##########################################################
  
  CCA_FL_O   <- array(dim = c(1, CCA_select, no_experiment))
  CCA_FL_std <- array(dim = c(1, CCA_select, no_experiment))
  
  Direction_CCA_Vec1_FL_O   <- array(dim = c(p1, CCA_select, no_experiment))
  Direction_CCA_Vec2_FL_O   <- array(dim = c(p2, CCA_select, no_experiment))
  Direction_CCA_Vec1_FL_std <- array(dim = c(p1, CCA_select, no_experiment))
  Direction_CCA_Vec2_FL_std <- array(dim = c(p2, CCA_select, no_experiment))
  
  ##########################################################
  ### RUN FUSED LASSO CCA
  ##########################################################
  
  ex <- 1
  
  FL <- CCA_by_fused_lasso(X_1 = X_1,
                           X_2 = X_2,
                           K   = CCA_select)
  
  ### Store ordered results
  Direction_CCA_Vec1_FL_O[,,ex] <- FL$CCA_Direction_X_1_fused_lasso_O
  Direction_CCA_Vec2_FL_O[,,ex] <- FL$CCA_Direction_X_2_fused_lasso_O
  CCA_FL_O[,,ex] <- sort(FL$CCA_fused_lasso_cal_O_CCA,
                         decreasing = TRUE)
  
  ### Store standard results
  Direction_CCA_Vec1_FL_std[,,ex] <- FL$CCA_Direction_X_1_fused_lasso_std
  Direction_CCA_Vec2_FL_std[,,ex] <- FL$CCA_Direction_X_2_fused_lasso_std
  CCA_FL_std[,,ex] <- sort(FL$CCA_fused_lasso_cal_std_CCA,
                           decreasing = TRUE)
  
  ##########################################################
  ### SIGN ALIGNMENT (CCA directions are sign-invariant)
  ### We align signs for interpretability
  ##########################################################
  
  dir_vec <- 1  # Analyze first canonical pair
  
  ### ----- STANDARD VERSION SIGN CORRECTION -----
  vec_1 <- FL$CCA_Direction_X_1_fused_lasso_std[, dir_vec]
  
  if(sign(max(vec_1)) == -1)
  {
    Direction_vector1_SC_std <- -FL$CCA_Direction_X_1_fused_lasso_std[, dir_vec]
    Direction_vector2_SC_std <- -FL$CCA_Direction_X_2_fused_lasso_std[, dir_vec]
  } else {
    Direction_vector1_SC_std <- FL$CCA_Direction_X_1_fused_lasso_std[, dir_vec]
    Direction_vector2_SC_std <- FL$CCA_Direction_X_2_fused_lasso_std[, dir_vec]
  }
  
  ##########################################################
  ### Chromosome 1 analysis (standard case)
  ##########################################################
  
  Direction_CCA_Vec2_fused_lasso_std_mean <- 
    cbind(Direction_vector2_SC_std, gene_chr)
  
  Direction_CCA_Vec2_fused_lasso_std_mean_chr <-
    subset(Direction_CCA_Vec2_fused_lasso_std_mean,
           gene_chr == 1)
  
  Direction_CCA_Vec2_fused_lasso_std_mean_chr_1_normed <-
    Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1] /
    norm_vec(Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1])
  
  ##########################################################
  ### SUMMARY STATISTICS
  ##########################################################
  
  Significant_Copy_Number_Std <-
    length(Direction_vector1_SC_std[Direction_vector1_SC_std != 0])
  
  Significant_Gene_Count_Std <-
    length(Direction_vector2_SC_std[Direction_vector2_SC_std != 0])
  
  Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_Std <-
    round(sum(Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1]^2), 4) * 100
  
  ##########################################################
  ### RETURN RESULTS
  ##########################################################
  
  Summary_Analysis <- list(
    CCA_FL_std = CCA_FL_std,
    CCA_FL_O = CCA_FL_O,
    Significant_Copy_Number_Std = Significant_Copy_Number_Std,
    Significant_Gene_Count_Std = Significant_Gene_Count_Std,
    Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_Std =
      Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_Std
  )
  
  return(Summary_Analysis)
}

    
    
    
















