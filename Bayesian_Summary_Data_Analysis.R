

##### Library #####

library(coda)

##### Functions#####

norm_vec <- function(x) sqrt(sum(x^2))

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


HDI_function<-function(col_data,true_val, alpha)
{
  
  col_data<-Re(col_data)
  
  true_val<-Re(true_val)
  
  ci_hdi <- quantile(col_data,c(alpha/2, 1-(alpha/2)))   # Calculate Credible Interval
  
  Ind_CCA<-ifelse(ci_hdi[[1]]<true_val & true_val< ci_hdi[[2]],1,0)  ## count number of 0 
  
  return(list(HDI=ci_hdi,Count=Ind_CCA))
  
}  


##########################################################################

#input<-

#load("20NDFSM_BC_Output.Rdata")

#Output = GCM; shrinkage_threshold = 0.2; Data_Set_1=Data_Set_1; gene_chr_250 = gene_chr_250; plot_indicator=FALSE

Analysis_Summary<-function(Output,Data, shrinkage_threshold=0.2,plot_indicator, plot_components=10, alpha.CC, alpha.dir, shrinkage_check=TRUE)
{


## Set up parameters###

thin_seq<-1:length(Output$CCA_MCMC[,1,])

p1<-nrow(Output$Direction_CCA_Vec1_MCMC[,1,]);p2<-nrow(Output$Direction_CCA_Vec2_MCMC[,1,]); CCA_select<-ncol(Output$CCA_MCMC)

iter_max<-iter_saved<-length(Output$CCA_MCMC[,1,])  # As these are thinned 

## Storing Zero Counts###

Zero_Count_CCA<-array(dim=c(1,CCA_select))

Zero_count_bounds<-array(dim=c(CCA_select,2))

Zero_Count_Direction_vec_1<-array(dim=c(p1,CCA_select))

Zero_Count_Direction_vec_2<-array(dim=c(p2,CCA_select))

### storing HDI###

HDI_CCA <- array(dim = c(CCA_select, 2))

HDI_Direction_vector_1 <- array(dim = c(p1, 2, CCA_select))

HDI_Direction_vector_2 <- array(dim = c(p2, 2, CCA_select))

########### Parameter of Interest Storage #### 

CCA_thin <- array(dim = c(1, CCA_select, length(thin_seq)))

Direction_CCA_Vec1_thin <- array(dim = c(p1, CCA_select, length(thin_seq)))

Direction_CCA_Vec2_thin <- array(dim = c(p2, CCA_select, length(thin_seq)))

############################################################

### Assigning Vectors###

CCA_thin <- Output$CCA_MCMC

Direction_CCA_Vec1_thin <- Output$Direction_CCA_Vec1_MCMC

Direction_CCA_Vec2_thin <- Output$Direction_CCA_Vec2_MCMC


################ Post processing ###########

Direction_vector1_SC <- array(dim = c(p1, CCA_select, iter_max))

Direction_vector2_SC <- array(dim = c(p2, CCA_select, iter_max))

j <- 1  # Selecting first element

dir_vec <- 1  #selecting the first vector

j <- 1  # Selecting first element

dir_vec <- 1  #selecting the first vector

Direction_vec_1_mean <- c()

i <- 1

for(i in 1:p1)
{
  Direction_vec_1_mean[i] <-
    Re(mean(abs(Direction_CCA_Vec1_thin[i, 1, 1]), na.rm = TRUE))
  
  
}                         


Direction_vec_2_mean<-c()

i<-1

for(i in 1:p2)
{
  Direction_vec_2_mean[i]<-Re(mean(abs(Direction_CCA_Vec2_thin[i,1,]),na.rm=TRUE))
  
}   

W1<-which(Direction_vec_1_mean==max(Direction_vec_1_mean),arr.ind=TRUE)

W2<-which(Direction_vec_2_mean==max(Direction_vec_2_mean),arr.ind=TRUE)


iter <- 1

for(iter in 1:iter_saved)

{  
  if(max(Direction_vec_1_mean)>max(Direction_vec_2_mean))
    
  { vec_1 <- Re(Direction_CCA_Vec1_thin[W1, dir_vec, iter])
  
  }else{ vec_1 <- Re(Direction_CCA_Vec2_thin[W2, dir_vec, iter])
  
  
  }
  
  if(sign(vec_1)==-1)
  {   
    
    Direction_vector1_SC[,dir_vec,iter]<--1*(Re(Direction_CCA_Vec1_thin[,dir_vec,iter]))
    
    Direction_vector2_SC[,dir_vec,iter]<--1*(Re(Direction_CCA_Vec2_thin[,dir_vec,iter]))
    
  } else{
    
    Direction_vector1_SC[,dir_vec,iter]<-Re(Direction_CCA_Vec1_thin[,dir_vec,iter])
    
    Direction_vector2_SC[,dir_vec,iter]<-Re(Direction_CCA_Vec2_thin[,dir_vec,iter]) 
    
  }  
  
  
  


Direction_CCA_Vec1_thin[,1,iter]<-Direction_vector1_SC[,1,iter]

Direction_CCA_Vec2_thin[,1,iter]<-Direction_vector2_SC[,1,iter]

}

##########################################################################

Direction_CCA_Vec1_GHS_MCMC_mean <- array(dim = c(p1, CCA_select))

Direction_CCA_Vec1_GHS_MCMC_mean<-Re(apply(Direction_CCA_Vec1_thin,c(1,2),mean,na.rm=TRUE)) 

Direction_CCA_Vec2_GHS_MCMC_mean<-array(dim=c(p2,CCA_select))

Direction_CCA_Vec2_GHS_MCMC_mean<-Re(apply(Direction_CCA_Vec2_thin,c(1,2),mean,na.rm=TRUE))   

Direction_CCA_Vec1_GHS_summary <- array(dim = c(p1, CCA_select))

Direction_CCA_Vec2_GHS_summary <- array(dim = c(p2, CCA_select))


Direction_CCA_Vec1_MCMC <- Direction_CCA_Vec1_thin

Direction_CCA_Vec2_MCMC <- Direction_CCA_Vec2_thin

for(j in 1:ncol(Direction_CCA_Vec1_GHS_MCMC_mean))
{
  
  Direction_CCA_Vec1_GHS_summary[,j]<-Direction_CCA_Vec1_GHS_MCMC_mean[,j]/norm_vec(Direction_CCA_Vec1_GHS_MCMC_mean[,j])
  
  
}

for(j in 1:ncol(Direction_CCA_Vec2_GHS_MCMC_mean))
{
  
  Direction_CCA_Vec2_GHS_summary[,j]<-Direction_CCA_Vec2_GHS_MCMC_mean[,j]/norm_vec(Direction_CCA_Vec2_GHS_MCMC_mean[,j])
  
  
}


#### Counting  based on credible intervals#####

j<-1

for(j in 1:CCA_select)
{
 
  HDI_cal<-HDI_function(col_data=CCA_thin[,j,],true_val=0, alpha=alpha.CC)
  
  Zero_Count_CCA[, j] <- HDI_cal$Count

  HDI_CCA[j, 1] <- HDI_cal$HDI[[1]]
  
  HDI_CCA[j, 2] <- HDI_cal$HDI[[2]]
  
 }

### Direction Vector 1####  

for(i in 1:p1)
{
  
  for(j in 1:CCA_select)
  {
    
    
    HDI_Direction_Vector_1<-HDI_function(col_data=Direction_CCA_Vec1_thin[i,j,],true_val=0, alpha=alpha.dir)
   
    Zero_Count_Direction_vec_1[i, j] <- HDI_Direction_Vector_1$Count
    
    HDI_Direction_vector_1[i, 1, j] <- HDI_Direction_Vector_1$HDI[[1]]
    
    HDI_Direction_vector_1[i, 2, j] <- HDI_Direction_Vector_1$HDI[[2]]
    
    
  }  
  
}


#### Direction Vector 2####    

for(i in 1:p2)
{
  
  for(j in 1:CCA_select)
  {
    
    HDI_Direction_Vector_2<-HDI_function(col_data=Direction_CCA_Vec2_thin[i,j,],true_val=0, alpha=alpha.dir)
    
    Zero_Count_Direction_vec_2[i, j] <- HDI_Direction_Vector_2$Count
    
    HDI_Direction_vector_2[i, 1, j] <- HDI_Direction_Vector_2$HDI[[1]]
    
    HDI_Direction_vector_2[i, 2, j] <- HDI_Direction_Vector_2$HDI[[2]]
    
  }  
  
}


###############################################

CCA_MCMC_summary <- apply(CCA_thin, c(1, 2), mean, na.rm = TRUE) 


###############  Plots  ##################


if(plot_indicator==TRUE)
{
 par(mfrow=c(2,5))

 i<-1

for(i in 1:CCA_select)
{
   plot(x=1:iter_saved,y=CCA_thin[,i,1:iter_saved],xlab='',ylab='', main=paste0("CC"," ", i),type="l")

  abline(h=CCA_MCMC_summary[i],col="red")

 }


 par(mfrow=c(2,5))
 
 dim(Direction_vector1_SC)
  
 for(i in 1:plot_components)
 {
   plot(x=1:iter_saved,y=Direction_CCA_Vec1_thin[i,1,1:iter_saved],xlab='',ylab='',main=paste0("Dir.Vec1"," ", i),type="l")

   abline(h= Direction_CCA_Vec1_GHS_summary[i,1],col="red")

 }

 par(mfrow=c(2,5))
  
 for(i in 1:plot_components)
 {
   plot(x=1:iter_saved,y=Direction_CCA_Vec2_thin[i,1,1:iter_saved],xlab='',ylab='',main=paste0("Dir.Vec2", " ", i),type="l")

   abline(h=Direction_CCA_Vec2_GHS_summary[i,1],col="red")

 }

 ### plotting log likelihood###
 
 
 par(mfrow=c(1,1))
 
 plot(x=1:iter_saved, y=apply(Output$Log_likelihood_MCMC_Grand,2,sum), xlab="Iterations", ylab="Log Likelihood ", main=" Log Likelihood Plot", type="l")
 
 par(mfrow=c(1,1))
 
 plot(x=1:iter_saved, y=Output$log_det_MCMC[1,], xlab="Iterations", ylab="Log Determinant ", main="Log Determinant Plot", type="l")
 
 
}

### Final Summary ####

if(shrinkage_check==TRUE){ 
Shrinkage_calculator<-length(CCA_thin[which(CCA_thin[,1,]<shrinkage_threshold)])/length(thin_seq)

if(Shrinkage_calculator > 0.5) {
  message("Run DFSM Model; Overshrinkage is suspected")
} else{
  message("Continue with NDFSM; No evidence of overshrinkage")
}

}else{Shrinkage_calculator=NA}

effective_size_first_CCA <- effectiveSize(CCA_thin[, 1, ])

effective_size_log_likelihood<-effectiveSize(apply(Output$Log_likelihood_MCMC_Grand,2,sum))
  
CCA_MCMC_summary <- apply(CCA_thin, c(1, 2), mean, na.rm = TRUE)

#### CCA Summary 

CC.data.frame <- data.frame(cbind(t(CCA_MCMC_summary),HDI_CCA))

names(CC.data.frame) <- c("Estimated_value", "Lower_bound", "Upper_bound")

round(CC.data.frame,3)

#Direction Vectors Summary of View 1 for first CC

V1.dir.hat <- Direction_CCA_Vec1_GHS_summary[,1]

V1.dir.data.frame<-data.frame(V1.dir.hat=V1.dir.hat,HDI_Direction_vector_1[,,1])

colnames(V1.dir.data.frame)<-c("Estimated_value", "Lower_bound", "Upper_bound")

#Direction Vectors Summary of View 2 for 1st CC

V2.dir.hat <- Direction_CCA_Vec2_GHS_summary[,1]

V2.dir.data.frame <- data.frame( Direction_CCA_Vec2_GHS_summary[,1],
                                 HDI_Direction_vector_2[,,1])
colnames(V2.dir.data.frame)<-c("Estimated_value", "Lower_bound", "Upper_bound")
rownames(V2.dir.data.frame) <- colnames(Data_Set_1$X_2_rna)

# View 1 Significant features

first_direction_significant_features_View1 <- 1-Zero_Count_Direction_vec_1[,1]

# View 2 Significant features

first_direction_significant_features_View2<- 1-Zero_Count_Direction_vec_2[,1]

Direction_CCA_Vec1_thin <- Direction_CCA_Vec1_thin[,1,]
Direction_CCA_Vec2_thin <- Direction_CCA_Vec2_thin[,1,]


Summary_Analysis<-list(Shrinkage_calculator=Shrinkage_calculator,first_direction_significant_features_View1=first_direction_significant_features_View1,first_direction_significant_features_View2=first_direction_significant_features_View2, V1.dir.data.frame=V1.dir.data.frame, V2.dir.data.frame=V2.dir.data.frame, Direction_Vector_View_1=Direction_CCA_Vec1_thin, Direction_Vector_View_2=Direction_CCA_Vec2_thin, effective_size_first_CCA=effective_size_first_CCA, effective_size_log_likelihood=effective_size_log_likelihood,CC.data.frame=CC.data.frame)

return(Summary_Analysis)

}





















