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


library(GFA)

library(MASS)

library(coda)


#### Norm 1

norm_vec <- function(x) sqrt(sum(x^2))


####### rinvgamma#####

rinvgamma<-function(n, shape, rate = 1, scale = 1/rate) 
{
  if (missing(rate) && !missing(scale)) 
    rate <- 1/scale
  1/rgamma(n, shape, rate)
  
  
}


### likelihood Calculation###

sigma_grand_cal<-function(A_1,A_2,phi_1,phi_2,p1,p2)
{
  m_11<-A_1%*%t(A_1)+phi_1
  m_12<-A_1%*%t(A_2)
  m_21<-A_2%*%t(A_1)
  m_22<-A_2%*%t(A_2)+phi_2
  sigma_grand<-matrix(nrow=p1+p2,ncol=p1+p2)
  sigma_grand[1:nrow(m_11),1:ncol(m_11)]<-m_11
  sigma_grand[1:nrow(m_12),ncol(m_11)+1:ncol(m_12)]<-m_12
  sigma_grand[nrow(m_11)+1:nrow(m_21),1:ncol(m_21)]<-m_21
  sigma_grand[nrow(m_11)+1:nrow(m_22),ncol(m_21)+1:ncol(m_22)]<-m_22
  return(sigma_grand)
  
}

### CCA function Calculation#######

CC_function<-function(X_1,X_2,A_1,A_2,phi_1,phi_2,CCA_select)
{
  m_11<-A_1%*%t(A_1)+phi_1   #Variance of X
  m_12<-A_1%*%t(A_2)         #Covariance between X and Y
  m_21<-A_2%*%t(A_1)        #Covariance between Y and X
  m_22<-A_2%*%t(A_2)+phi_2   #Variance of Y
  
  #### finding the square  root of the matrices using diagonalization
  
  eigen_x<-eigen(m_11)
  
  eigen_y<-eigen(m_22)
  
  sqrm_11<-(eigen_x$vectors)%*%diag(eigen_x$values^-0.5)%*%solve(eigen_x$vectors)
  
  sqrm_22<-(eigen_y$vectors)%*%diag(eigen_y$values^-0.5)%*%solve(eigen_y$vectors)
  
  ########## Finding the weights /directions of CCA########## 
  
  m1<-sqrm_11%*%m_12%*%solve(m_22)%*%t(m_12)%*%sqrm_11
  
  m2<-sqrm_22%*%m_21%*%solve(m_11)%*%t(m_21)%*%sqrm_22
  
  m1e<-eigen(m1) ; m2e<-eigen(m2)
  
  ##### Canonical Correlations#######
  Canonical_Correaltions<-c()
  
  Canonical_Correaltions_corrected<-c()
  
  for(j in 1:min(nrow(m1e$vectors),nrow(m2e$vectors)))
  {
    
    Canonical_Correaltions[j]<-Re(m1e$vectors[,j])%*%sqrm_11%*%m_12%*%sqrm_22%*%Re(m2e$vectors[,j])
    
    
    
  }
  
  
  for(i in 1:length(Canonical_Correaltions))
  {
    if(sign(Canonical_Correaltions[i])==-1)
    {   
      
      Canonical_Correaltions[i]<-abs(Canonical_Correaltions[i])
      
      m2e$vectors[,i]<--1*m2e$vectors[,i]
      
      
    }
    
    
  } 
  
  return(list(CCA_grand=Canonical_Correaltions[1:CCA_select],Direction_CCA_Vec1_grand=m1e$vectors[,1:CCA_select],Direction_CCA_Vec2_grand=m2e$vectors[,1:CCA_select]))
  
}

##################################################


CC_function_grand_2<-function(Covariance_matrix_Grand,p1,p2,CCA_select)
{
  
  sigma_grand<-Covariance_matrix_Grand
  m_11<-sigma_grand[1:p1,1:p1]                     				          #Variance of X
  m_12<-sigma_grand[1:p1,p1+1:p2]                                   #Covariance between X and Y
  m_21<-sigma_grand[p1+1:p2,1:p1]                                		#Covariance between Y and X
  m_22<-sigma_grand[p1+1:p2,p1+1:p2]                        				#Variance of Y
  
  
  #### finding the square  root of the matrices using diagonalization
  
  eigen_x<-eigen(m_11)
  
  eigen_y<-eigen(m_22)
  
  sqrm_11<-(eigen_x$vectors)%*%diag(eigen_x$values^-0.5)%*%solve(eigen_x$vectors)
  
  sqrm_22<-(eigen_y$vectors)%*%diag(eigen_y$values^-0.5)%*%solve(eigen_y$vectors)
  
  ########## Finding the weights /directions of CCA########## 
  
  m1<-sqrm_11%*%m_12%*%solve(m_22)%*%t(m_12)%*%sqrm_11
  
  m2<-sqrm_22%*%m_21%*%solve(m_11)%*%t(m_21)%*%sqrm_22
  
  m1e<-eigen(m1) ; m2e<-eigen(m2)
  
  ##### Canonical Correlations#######
  Canonical_Correaltions<-c()
  
  for(j in 1:CCA_select)
  {
    Canonical_Correaltions[j]<-Re(m1e$vectors[,j])%*%sqrm_11%*%m_12%*%sqrm_22%*%Re(m2e$vectors[,j])
    
  }
  
  
  
  for(i in 1:length(Canonical_Correaltions))
  {
    if(sign(Canonical_Correaltions[i])==-1)
    {   
      
      Canonical_Correaltions[i]<-abs(Canonical_Correaltions[i])
      
      m2e$vectors[,i]<--1*m2e$vectors[,i]
      
      
    }
    
    
  }
  
  
  ### return the direction vectors###
  
  return(list(CCA_grand=Canonical_Correaltions,Direction_CCA_Vec1_grand=Re(m1e$vectors[,1:CCA_select]),Direction_CCA_Vec2_grand=Re(m2e$vectors[,1:CCA_select])))
  
}



##### GFA function### 


GFA_function<-function(d,X_1,X_2,burn_iter,mcmc_iter,p1,p2,iter_saved,CCA_select)
{
  
  set.seed(666)  
  
  p1<-ncol(X_1); p2<-ncol(X_2)
  
  CCA_GFA<-array(dim=c(1,CCA_select,iter_saved))
  
  Direction_CCA_vec1_GFA<-array(dim=c(p1,CCA_select,iter_saved))
  
  Direction_CCA_vec2_GFA<-array(dim=c(p2,CCA_select,iter_saved))
  
  default_values<-getDefaultOpts()
  
  default_values$iter.max<-mcmc_iter
  
  default_values$iter.saved<-iter_saved
  
  default_values$iter.burnin<-burn_iter
  
  GFA_analysis<-gfa(Y=list(X_1,X_2),K=50,opts=default_values)
  
  K_fit<-GFA_analysis$K
  
  for(iter in 1:iter_saved)
  {
    
    CCA_GFA_2<-CC_function(X_1=X_1,X_2=X_2,A_1=GFA_analysis$posterior$W[iter,1:p1,],A_2=GFA_analysis$posterior$W[iter,p1+1:p2,],phi_1=diag(1/GFA_analysis$posterior$tau[iter,1:p1],p1),phi_2=diag(1/GFA_analysis$posterior$tau[iter,p1+1:p2],p2),CCA_select)
    
    CCA_GFA[1,,iter]<-abs(CCA_GFA_2$CCA)
    
    Direction_CCA_vec1_GFA[,,iter]<-Re(CCA_GFA_2$Direction_CCA_Vec1)
    
    Direction_CCA_vec2_GFA[,,iter]<-Re(CCA_GFA_2$Direction_CCA_Vec2)
    
    
  }
  
  
  return(list(CCA_MCMC=CCA_GFA,Direction_CCA_Vec1_MCMC=Direction_CCA_vec1_GFA,Direction_CCA_Vec2_MCMC=Direction_CCA_vec2_GFA,K_fit=K_fit))
  
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































