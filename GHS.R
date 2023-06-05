
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

# no_patients<-nrow(X_1) #No of Experiments


###
###   LOAD AN R LIBRARY (SIMULATES MULTIVARIATE NORMAL)
###

library(MASS)

library(mvtnorm)

###
###   FUNCTIONS TO CALCULATE CORRELATIONS
###



###########mvrnorm####

mvrnorm<-function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (EISPACK) 
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
  
  
}


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

######################################################################################

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
  {Canonical_Correaltions[j]<-Re(m1e$vectors[,j])%*%sqrm_11%*%m_12%*%sqrm_22%*%Re(m2e$vectors[,j])
  
  }
  
  
 
  for(i in 1:length(Canonical_Correaltions))
  {
    if(sign(Canonical_Correaltions[i])==-1)
    {   
      
      Canonical_Correaltions[i]<-abs(Canonical_Correaltions[i])
    
    m2e$vectors[,j]<--1*m2e$vectors[,j]
    
    
    }
    
    
  }
  
  
  ### return the direction vectors###
  
  return(list(CCA_grand=Canonical_Correaltions,Direction_CCA_Vec1_grand=m1e$vectors[,1:CCA_select],Direction_CCA_Vec2_grand=m2e$vectors[,1:CCA_select]))
  
}



##### Graphical Horseshoe Calculation######

GHS_modified_2<-function(X,Lambda_sq,Nu,Omega,Sigma,tau_sq,xi,iter){
  
  n<-nrow(X)
  
  S<-t(X)%*%(X)
  p<-nrow(S)
  omega_save <- array(0,c(p,p))
  ind_all <-matrix(0,p-1,p)
  for(i in 1:p){
    if(i==1){ind = c(2:p)
    }else if(i==p){ind = c(1:(p-1))
    }else{ind = c(1:(i-1),(i+1):p)}
    ind_all[,i] <- ind
  }
  
  for (i in 1:p){
    ind <- ind_all[,i]
    sigma_11 <- Sigma[ind,ind]
    sigma_12 <- Sigma[ind,i]
    sigma_22 <- Sigma[i,i]
    s_21 <- S[ind,i]
    s_22 <- S[i,i]
    lambda_sq_12 <- Lambda_sq[ind,i]
    nu_12 <- Nu[ind,i]
    gamma <- rgamma(1,shape=((n/2)+1),rate=s_22/2)
    inv_Omega_11 <- sigma_11 - ((sigma_12%*%t(sigma_12))/(sigma_22))
    inv_C <- (s_22*inv_Omega_11) + diag(1/(lambda_sq_12*tau_sq),p-1)
    inv_C_chol <- chol(inv_C)
    mu_i <- -chol2inv(inv_C_chol)%*%s_21
    beta <- mu_i+ solve(inv_C_chol,rnorm(p-1,0,1))
    omega_12 <- beta
    omega_22 <- gamma + t(beta)%*%inv_Omega_11%*%beta
    rate <- (1/nu_12)+((omega_12^2)/(2*tau_sq))
    lambda_sq_12 <- rinvgamma(p-1,shape=1,rate=rate)
    nu_12 <- rinvgamma(p-1,shape=1,rate=(1+(1/lambda_sq_12)))
    Omega[i,ind] <- omega_12
    Omega[ind,i] <- omega_12
    Omega[i,i] <- omega_22
    temp <- inv_Omega_11%*%beta
    Sigma_11 <- inv_Omega_11 + temp%*%t(temp)/gamma
    sigma_12 <- -temp/gamma
    sigma_22 <- 1/gamma
    Sigma[ind,ind] <- Sigma_11
    Sigma[i,i] <- sigma_22
    Sigma[i,ind] <- sigma_12
    Sigma[ind,i] <- sigma_12
    Lambda_sq[i,ind] <- lambda_sq_12
    Lambda_sq[ind,i] <- lambda_sq_12
    Nu[i,ind] <- nu_12
    Nu[ind,i] <- nu_12
  }
  
  
  omega_vector <- Omega[lower.tri(Omega,diag = FALSE)]
  lambda_sq_vector <- Lambda_sq[lower.tri(Lambda_sq,diag = FALSE)]
  rate1 <- 1/xi + sum((omega_vector^2)/(2*lambda_sq_vector))
  tau_sq <- rinvgamma(1,shape=((p*(p-1)/2)+1)/2, rate=rate1)
  xi <- rinvgamma(1,shape=1,rate=(1+ 1/tau_sq))
  return(list(Omega=Omega,Lambda_sq=Lambda_sq,Nu=Nu,Tau_sq= tau_sq, Sigma=Sigma,xi=xi))
  

}


##### Grand GHS method####### 

GHS_CC<-function(X_1,X_2,burn_iter,mcmc_iter,thin,CCA_select)
{
  
  X_grand <- cbind(X_1=X_1,X_2=X_2)
  
  p1<-ncol(X_1) ;p2<-ncol(X_2) 
  
  ### MCMC initialization ###
  
  Phi_grand_GHS_iter<-array(dim=c(p1+p2,p1+p2))
  
  Lambda_sq_GHS_grand<-Nu_GHS_grand<-matrix(1,nrow=p1+p2,ncol=p1+p2)
  
  Omega_grand_GHS_iter<-Sigma_grand<-array(dim=c(p1+p2,p1+p2))
  
  X_grand_delta<-matrix(nrow=no_patients,ncol=p1+p2)
  
  Phi_grand_GHS_iter<-0.5*cov(X_grand)+0.5*diag(p1+p2)
  
  Omega_grand_GHS_iter<-solve(Phi_grand_GHS_iter)
  
  Sigma_grand_GHS<-diag(1,p1+p2)
  
  tau_sq_GHS_grand<-1
  
  xi_GHS_grand<-1
  
  Dpost_grand<-1
  
  sigma<-100
  
  #############  Thinning Storage #########
  
  iter_max<-burn_iter+mcmc_iter
  
  thin_seq<-seq(burn_iter+thin,iter_max,thin)
  
  Omega_grand_thin<-array(dim=c(p1+p2,p1+p2,length(thin_seq)))
  
  Phi_grand_det_thin<-array(dim=c(1,length(thin_seq)))
 
  log_det_MCMC_thin<-array(dim=c(1,length(thin_seq)))
  
  Log_likelihood_MCMC_Grand_thin<-array(dim=c(no_patients,length(thin_seq)))
  
  CCA_thin<-array(dim=c(1,CCA_select,length(thin_seq)))
  
  Direction_CCA_Vec1_thin<-array(dim=c(p1,CCA_select,length(thin_seq)))
  
  Direction_CCA_Vec2_thin<-array(dim=c(p2,CCA_select,length(thin_seq)))
  
  likelihood<-array(dim=c(no_patients,length(thin_seq)))
  
  i<-1
  
  iter<-2
  
  for(iter in 2:iter_max)
  {
    
    ########## Mu Calculation####
    
    E_mu_cal<-no_patients*(Omega_grand_GHS_iter)+diag(1/sigma,p1+p2)    
    
    variance_mu<-solve(E_mu_cal)  
    
    mean_mu<-variance_mu%*%(Omega_grand_GHS_iter)%*%(apply(X_grand,2,sum))
    
    mu_grand<-mvrnorm(1,mean_mu,variance_mu)    
    
    for(patient in 1:no_patients)
    {
      X_grand_delta[patient,]<-X_grand[patient,]-mu_grand
      
      
    }
    
    GHS_Run_3<-GHS_modified_2(X=X_grand_delta,Lambda_sq=Lambda_sq_GHS_grand,Nu=Nu_GHS_grand,Omega=Omega_grand_GHS_iter,Sigma=Phi_grand_GHS_iter,tau_sq =tau_sq_GHS_grand,xi=xi_GHS_grand,iter)
    
    Phi_grand_GHS_iter<-GHS_Run_3$Sigma 
    
    Omega_grand_GHS_iter<-GHS_Run_3$Omega
    
    Lambda_sq_GHS_grand<-GHS_Run_3$Lambda_sq; Nu_GHS_grand<-GHS_Run_3$Nu; tau_sq_GHS_grand<-GHS_Run_3$Tau_sq ;xi_GHS_grand<-GHS_Run_3$xi
    
    X_grand=cbind(X_1=X_1,X_2=X_2)
    
    if(iter==thin_seq[i])
    {
    
    Omega_grand_thin[,,i]<-Omega_grand_GHS_iter
    
    likelihood[,i]<-dmvnorm(X_grand,mean=mu_grand,sigma=Phi_grand_GHS_iter,log=TRUE)
      
    cc_cal_grand<-CC_function_grand_2(Covariance_matrix_Grand=Phi_grand_GHS_iter,p1,p2,CCA_select)
    
    CCA_thin[,,i]<-cc_cal_grand$CCA
    
    Direction_CCA_Vec1_thin[,,i]<-cc_cal_grand$Direction_CCA_Vec1
    
    Direction_CCA_Vec2_thin[,,i]<-cc_cal_grand$Direction_CCA_Vec2
    
    log_det_MCMC_thin[,i]<-determinant(Phi_grand_GHS_iter,logarithm = TRUE)$modulus[1]
    
    i=i+1
    
    }
    
    
   ##############################################################################
    if(iter==iter_max)
   { print(paste0("Iteration Number",iter))}
    
  }
  
  
 return(list(Omega_grand=Omega_grand_thin,CCA_MCMC=CCA_thin,Direction_CCA_Vec1_MCMC=Direction_CCA_Vec1_thin,Direction_CCA_Vec2_MCMC=Direction_CCA_Vec2_thin,log_det_MCMC=log_det_MCMC_thin,likelihood_MCMC=likelihood))
  
}


######## Grand  GHS#################
# 
# GHS_output<-GHS_grand_CCA_2(X_grand=cbind(X_1=X_1,X_2=X_2),p1,p2,burn_iter,mcmc_iter,thin=7,CCA_select=10)
# 
# #### Storing CCA outputs####
# 
# filename_Grand_GHS<-paste(20,"GHS_BC_Output",sep="_")
# 
# save(GHS_output,file=paste(filename_Grand_GHS,".Rdata",sep=""))
# 
# 


