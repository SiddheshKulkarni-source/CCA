###
###   RECEIVE ARGUMENTS FROM CLUSTER
###


# args<-commandArgs()
# 
# ex<-as.numeric(args)[6]   
# 
# cat(ex,sep="\n")

ex <-1
###
###   LOAD AN R LIBRARY (SIMULATES MULTIVARIATE NORMAL)
###


library(MASS)


###
###   FUNCTIONS TO Generate the data and 
###

###### AR Structure###########

#### Toeplitz Structure########

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  rho^exponent
}

############  Data Generation Function #####


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

Data_genaration1<-function(no_patients,ex,rho=TRUE,scaling=1)
{
  d <- 1
  p1 <- 100
  p2 <- 50
  rho1 <- .4*rho
  rho2 <- .2*rho

  set.seed(1234*ex)
  ## Latent variable generation
  
  Z<-mvrnorm(no_patients,rep(0,d),diag(d))
  
  ### Generation of PSi_1 and Psi_2 as initial data geneartion#####
  psi_1<-ar1_cor(p1,rho1)
  
  psi_2<-ar1_cor(p2,rho2)
  
  ##### Projection matrix generation#############
  
  A_1_data<- diag(0,nrow=p1,ncol=d)  #Fix pattern generation
  A_1_data[c(1,11,21)] <- 1
  A_1_data <- scaling * A_1_data
  
  A_2_data<- diag(0,nrow=p2,ncol=d)  #Fix pattern generation
  A_2_data[c(1,11)] <- c(1,-1)
  A_2_data <- scaling * A_2_data 
  
  ##### Drawing Mu_1 and Mu_2############
  
  Mu_1_data<-rnorm(p1)*0
  
  Mu_2_data<-rnorm(p2)*0
  
  ### The Simulated Data Set##########
  X_1<-matrix(nrow=no_patients,ncol=p1)
  X_2<-matrix(nrow=no_patients,ncol=p2)
  
  for(i in 1:no_patients)
  {
    X_1[i,]<-mvrnorm(1,Mu_1_data+A_1_data%*%Z[i,],psi_1)
    X_2[i,]<-mvrnorm(1,Mu_2_data+A_2_data%*%Z[i,],psi_2)
    
  }
  
  X_1_central<-matrix(nrow=no_patients,ncol=p1)
  
  X_2_central<-matrix(nrow=no_patients,ncol=p2)
  
  for(i in 1:no_patients)
  {
    X_1_central[i,]<-X_1[i,]-Mu_1_data+A_1_data%*%Z[i,]
    
    X_2_central[i,]<-X_2[i,]-Mu_2_data+A_2_data%*%Z[i,]
    
  }  
  
 return(list(Z=Z,A_1_data=A_1_data,A_2_data=A_2_data,Mu_1_data=Mu_1_data,Mu_2_data=Mu_2_data,X_1=X_1, X_1_central= X_1_central,X_2=X_2, X_2_central= X_2_central,X_grand=cbind(X_1,X_2),psi_1=psi_1,psi_2=psi_2,p1=p1,p2=p2,rho1=rho1,rho2=rho2,d=d))
  
  
}



Data_genaration2<-function(no_patients,ex,rho=TRUE)
{
  d <- 10
  p1 <- 100
  p2 <- 50
  rho1 <- .4*rho
  rho2 <- .2*rho

  ##### Projection matrix generation#############
  
 set.seed(112358)  
  A_1_data<- diag(0,nrow=p1,ncol=d)  #Fix pattern generation
  A_1_data[c(1,11,21),1] <- 1
  A_1_data[,-1] <- rbinom( (d-1)*p1, 1, .05) * rnorm((d-1)*p1)
  

	#(2*rbinom((d-1)*p1, 1, .5)-1)
  A_1_data <- A_1_data %*% diag( c(1,seq(.45,.05,length=d-1) ))

  
  A_2_data<- diag(0,nrow=p2,ncol=d)  #Fix pattern generation
  A_2_data[c(1,11),1] <- c(1,-1)
  A_2_data[,-1] <- rbinom( (d-1)*p2, 1, .05) * rnorm((d-1)*p2)
	#(2*rbinom((d-1)*p2, 1, .5)-1)
  A_2_data <- A_2_data %*% diag( c(1,seq(.45,.05,length=d-1) ))


  ##### Drawing Mu_1 and Mu_2############
  
  Mu_1_data<-rnorm(p1)*0
  
  Mu_2_data<-rnorm(p2)*0
  
  set.seed(1234*ex)
  ## Latent variable generation
  
  Z<-mvrnorm(no_patients,rep(0,d),diag(d))
  
  ### Generation of PSi_1 and Psi_2 as initial data geneartion#####
  psi_1<-ar1_cor(p1,rho1)
  
  psi_2<-ar1_cor(p2,rho2)
  
  ### The Simulated Data Set##########
  X_1<-matrix(nrow=no_patients,ncol=p1)
  X_2<-matrix(nrow=no_patients,ncol=p2)
  
  for(i in 1:no_patients)
  {
    X_1[i,]<-mvrnorm(1,Mu_1_data+A_1_data%*%Z[i,],psi_1)
    X_2[i,]<-mvrnorm(1,Mu_2_data+A_2_data%*%Z[i,],psi_2)
    
  }
  
  X_1_central<-matrix(nrow=no_patients,ncol=p1)
  X_2_central<-matrix(nrow=no_patients,ncol=p2)
  
    for(i in 1:no_patients)
  {
    X_1_central[i,]<-X_1[i,]-Mu_1_data+A_1_data%*%Z[i,]
    
    X_2_central[i,]<-X_2[i,]-Mu_2_data+A_2_data%*%Z[i,]
    
  }  
  
  return(list(Z=Z,A_1_data=A_1_data,A_2_data=A_2_data,Mu_1_data=Mu_1_data,Mu_2_data=Mu_2_data,X_1=X_1, X_1_central= X_1_central,X_2=X_2, X_2_central= X_2_central,X_grand=cbind(X_1,X_2),psi_1=psi_1,psi_2=psi_2,p1=p1,p2=p2,d=d,rho1=rho1,rho2=rho2))
  
  
}



CC_function<-function(X_1,X_2,A_1,A_2,phi_1,phi_2,p1,p2)
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
  {Canonical_Correaltions[j]<-Re(m1e$vectors[,j])%*%sqrm_11%*%m_12%*%sqrm_22%*%Re(m2e$vectors[,j])
 
  }
  
 
  return(list(CCA=Canonical_Correaltions,Direction_CCA_Vec1=round(m1e$vectors,3),Direction_CCA_Vec2=round(m2e$vectors,3)))
  
}


#### Actual Calculations#### 

### Option 1:
DG1<- Data_genaration1(no_patients=300,ex)
CC_cal_Data1<- CC_function(X_1=DG1$X_1,X_2 = DG1$X_2,A_1 =DG1$A_1_data,A_2 = DG1$A_2_data, phi_1 = DG1$psi_1, phi_2 = DG1$psi_2,p1=DG1$p1,p2=DG$p2)
Data_Set_1<-list(DG1,CC_cal_Data1)

### Option 2:
DG2<- Data_genaration1(no_patients=50,ex)
CC_cal_Data2<- CC_function(X_1=DG2$X_1,X_2 = DG2$X_2,A_1=DG2$A_1_data,A_2 = DG2$A_2_data, phi_1 = DG2$psi_1, phi_2 = DG2$psi_2,p1=DG2$p1,p2=DG$p2)
Data_Set_2<-list(DG2,CC_cal_Data2)

### Option 3:
DG3<- Data_genaration1(no_patients=300,ex, rho=FALSE)
CC_cal_Data3<- CC_function(X_1=DG3$X_1,X_2 = DG3$X_2,A_1 =DG3$A_1_data,A_2 = DG3$A_2_data, phi_1 = DG3$psi_1, phi_2 = DG3$psi_2,p1=DG3$p1,p2=DG3$p2)
Data_Set_3<-list(DG3,CC_cal_Data3)

### Option 4:
DG4<- Data_genaration1(no_patients=50,ex, rho=FALSE)
CC_cal_Data4<- CC_function(X_1=DG4$X_1,X_2 = DG4$X_2,A_1=DG4$A_1_data,A_2 = DG4$A_2_data, phi_1 = DG4$psi_1, phi_2 = DG4$psi_2,p1=DG4$p1,p2=DG4$p2)
#CC_grand4<-CC_function_grand(X_grand_GHS=DG4$X_grand,p1=DG4$p1,p2=DG4$p2)  ## CC grand is giving problem
Data_Set_4<-list(DG4,CC_cal_Data4)

### Option 5:
DG5<- Data_genaration2(no_patients=300,ex)
CC_cal_Data5<- CC_function(X_1=DG5$X_1,X_2 = DG5$X_2,A_1 =DG5$A_1_data,A_2 = DG5$A_2_data, phi_1 = DG5$psi_1, phi_2 = DG5$psi_2,p1=DG5$p1,p2=DG5$p2)
Data_Set_5<-list(DG5,CC_cal_Data5)

### Option 6:
DG6<- Data_genaration2(no_patients=50,ex)
CC_cal_Data6<- CC_function(X_1=DG6$X_1,X_2 = DG6$X_2,A_1 =DG6$A_1_data,A_2 = DG6$A_2_data, phi_1 = DG6$psi_1, phi_2 = DG6$psi_2,p1=DG6$p1,p2=DG6$p2)
Data_Set_6<-list(DG6,CC_cal_Data6)

### Option 7:
DG7<- Data_genaration1(no_patients=300,ex, scaling=.59)
CC_cal_Data7<- CC_function(X_1=DG7$X_1,X_2 = DG7$X_2,A_1 =DG7$A_1_data,A_2 = DG7$A_2_data, phi_1 = DG7$psi_1, phi_2 = DG7$psi_2,p1=DG7$p1,p2=DG7$p2)
Data_Set_7<-list(DG7,CC_cal_Data7)




