

###
###   LOAD AN R LIBRARY (SIMULATES MULTIVARIATE NORMAL)
###

library(PMA)

###### parameters 

#load("1_Breast_cancer_20th_Dec_ch_1_50.Rdata")

Fused_lasso_analysis<-function(Data_Set_1=Data_Set_1, CCA_select=10, gene_chr_250=gene_chr_250)
{
X_1<-Data_Set_1$X_1_dna

X_2<-Data_Set_1$X_2_rna

no_patients<-nrow(X_1) #No of Experiments

no_experiment<-1

p1<-ncol(X_1)  # dimension of view 1

p2<-ncol(X_2)    # dimension of view 1

CCA_select<-CCA_select

gene_chr<-gene_chr_250

### Function

norm_vec <- function(x) sqrt(sum(x^2))

CCA_by_fused_lasso<-function(X_1,X_2,K)
{
  
  CCA_fused_lasso_cal_O<-CCA(X_1,X_2,typex = "ordered",typez = "ordered",K=K)
  
  CCA_fused_lasso_cal_O_CCA<-CCA_fused_lasso_cal_O$cors
  
  CCA_Direction_X_1_fused_lasso_O<-CCA_fused_lasso_cal_O$u  # These are Direction vectors
  
  CCA_Direction_X_2_fused_lasso_O<-CCA_fused_lasso_cal_O$v  # These  are direction vectors  
  
  CCA_fused_lasso_cal_std<-CCA(X_1,X_2,typex = "standard",typez = "standard",K=K)
  
  CCA_fused_lasso_cal_std_CCA<-CCA_fused_lasso_cal_std$cors
  
  CCA_Direction_X_1_fused_lasso_std<-CCA_fused_lasso_cal_std$u  # These are Direction vectors
  
  CCA_Direction_X_2_fused_lasso_std<-CCA_fused_lasso_cal_std$v  # These  are direction vectors  
  
  return(list(CCA_fused_lasso_cal_O_CCA=CCA_fused_lasso_cal_O_CCA,CCA_Direction_X_1_fused_lasso_O=CCA_Direction_X_1_fused_lasso_O,CCA_Direction_X_2_fused_lasso_O=CCA_Direction_X_2_fused_lasso_O,CCA_fused_lasso_cal_std_CCA=CCA_fused_lasso_cal_std_CCA,CCA_Direction_X_1_fused_lasso_std=CCA_Direction_X_1_fused_lasso_std,CCA_Direction_X_2_fused_lasso_std=CCA_Direction_X_2_fused_lasso_std))
  
}


#########################################################

CCA_Data<-array(dim=c(1,CCA_select,no_experiment))

#### Direction vector Storage###

Direction_CCA_Vec1_Data<-array(dim=c(p1,CCA_select,no_experiment))

Direction_CCA_Vec2_Data<-array(dim=c(p2,CCA_select,no_experiment))

#### Fused Lasso#### 

CCA_FL_O<-array(dim=c(1,CCA_select,no_experiment))      

CCA_FL_std<-array(dim=c(1,CCA_select,no_experiment))   


Direction_CCA_Vec1_FL_O<-array(dim=c(p1,CCA_select,no_experiment)) ## # writing according to the MCMC structure

Direction_CCA_Vec2_FL_O<-array(dim=c(p2,CCA_select,no_experiment)) ## # writing according to the MCMC structure


Direction_CCA_Vec1_FL_std<-array(dim=c(p1,CCA_select,no_experiment)) ## # writing according to the MCMC structure

Direction_CCA_Vec2_FL_std<-array(dim=c(p2,CCA_select,no_experiment)) ## # writing according to the MCMC structure


########################################


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

#########################################################

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

  ############# Function Run ###

    ex<-1
    
    FL<-CCA_by_fused_lasso(X_1=X_1,X_2=X_2,K=CCA_select)
    
    Direction_CCA_Vec1_FL_O[,,ex]<-FL$CCA_Direction_X_1_fused_lasso_O
    
    Direction_CCA_Vec2_FL_O[,,ex]<-FL$CCA_Direction_X_2_fused_lasso_O
    
    CCA_FL_O[,,ex]<-sort(FL$CCA_fused_lasso_cal_O_CCA, decreasing = TRUE)
  
    Direction_CCA_Vec1_FL_std[,,ex]<-FL$CCA_Direction_X_1_fused_lasso_std
    
    Direction_CCA_Vec2_FL_std[,,ex]<-FL$CCA_Direction_X_2_fused_lasso_std
    
    CCA_FL_std[,,ex]<-sort(FL$CCA_fused_lasso_cal_std_CCA, decreasing=TRUE)
    
    ################## Direction Vector Alignments###
    
    Direction_vector1_SC_O<-c()
    
    Direction_vector2_SC_O<-c()
    
    dir_vec<-1 
    
    W1<-which(FL$CCA_Direction_X_1_fused_lasso_O[,dir_vec]==max(Re(FL$CCA_Direction_X_1_fused_lasso_O[,dir_vec])),arr.ind=TRUE)
    
    W2<-which(FL$CCA_Direction_X_2_fused_lasso_O[,dir_vec]==max(Re(FL$CCA_Direction_X_2_fused_lasso_O[,dir_vec])),arr.ind=TRUE)
    
    
    if(max(Re(FL$CCA_Direction_X_1_fused_lasso_O[,dir_vec]))>max(Re(FL$CCA_Direction_X_2_fused_lasso_O[,dir_vec])))
      
    { vec_1<-Re(FL$CCA_Direction_X_1_fused_lasso_O[W1,dir_vec])
    
    }else{ vec_1<-Re(FL$CCA_Direction_X_2_fused_lasso_O[W2,dir_vec])
    
    
    }
    
    
    if(sign(vec_1[1])==-1)
      
    {   Direction_vector1_SC_O<--1*(Re(FL$CCA_Direction_X_1_fused_lasso_O[,dir_vec]))
    
    Direction_vector2_SC_O<--1*(Re(FL$CCA_Direction_X_2_fused_lasso_O[,dir_vec]))
    
    } else{
      
      Direction_vector1_SC_O<-(Re(FL$CCA_Direction_X_1_fused_lasso_O[,dir_vec]))
      
      Direction_vector2_SC_O<-(Re(FL$CCA_Direction_X_2_fused_lasso_O[,dir_vec]))
      
    }  
    
    
    
    ##### Direction vector summary####
    
    Direction_CCA_Vec1_fused_lasso_O_mean<-array(dim=c(p1,p1))
    
    Direction_CCA_Vec1_fused_lasso_O_mean<-Re(FL$CCA_Direction_X_1_fused_lasso_O)
    
    Direction_CCA_Vec1_fused_lasso_O_mean[,1]<-Direction_vector1_SC_O
    
    
    #############  Direction Vector 2 Summary ########
    
    Direction_CCA_Vec2_fused_lasso_O_mean<-array(dim=c(p2,p2))
    
    Direction_CCA_Vec2_fused_lasso_O_mean<-Re(FL$CCA_Direction_X_1_fused_lasso_O)
    
    Direction_CCA_Vec2_fused_lasso_O_mean[,1]<-Direction_vector1_SC_O
    
    
    
  ##################################################################################
    
    Direction_vector1_SC_std<-c()
    
    Direction_vector2_SC_std<-c()
    
    dir_vec<-1 
    
    W1<-which(FL$CCA_Direction_X_1_fused_lasso_std[,dir_vec]==max(Re(FL$CCA_Direction_X_1_fused_lasso_std[,dir_vec])),arr.ind=TRUE)
    
    W2<-which(FL$CCA_Direction_X_2_fused_lasso_std[,dir_vec]==max(Re(FL$CCA_Direction_X_2_fused_lasso_std[,dir_vec])),arr.ind=TRUE)
    
    if(max(Re(FL$CCA_Direction_X_1_fused_lasso_std[,dir_vec]))>max(Re(FL$CCA_Direction_X_2_fused_lasso_std[,dir_vec])))
      
    { vec_1<-Re(FL$CCA_Direction_X_1_fused_lasso_std[W1,dir_vec])
    
    }else{ vec_1<-Re(FL$CCA_Direction_X_2_fused_lasso_std[W2,dir_vec])
    
    
    }
    
    
    if(sign(vec_1[1])==-1)
      
    {   Direction_vector1_SC_std<--1*(Re(FL$CCA_Direction_X_1_fused_lasso_std[,dir_vec]))
    
    Direction_vector2_SC_std<--1*(Re(FL$CCA_Direction_X_2_fused_lasso_std[,dir_vec]))
    
    } else{
      
      Direction_vector1_SC_std<-(Re(FL$CCA_Direction_X_1_fused_lasso_std[,dir_vec]))
      
      Direction_vector2_SC_std<-(Re(FL$CCA_Direction_X_2_fused_lasso_std[,dir_vec]))
      
    }  
    
    ##### Direction vector summary####
    
    Direction_CCA_Vec1_fused_lasso_std_mean<-array(dim=c(p1,p1))
    
    Direction_CCA_Vec1_fused_lasso_std_mean<-Re(FL$CCA_Direction_X_1_fused_lasso_std)
    
    Direction_CCA_Vec1_fused_lasso_std_mean[,1]<-Direction_vector1_SC_std
    
    
    #############  Direction Vector 2 Summary ########
    
    Direction_CCA_Vec2_fused_lasso_std_mean<-array(dim=c(p2,p2))
    
    Direction_CCA_Vec2_fused_lasso_std_mean<-Re(FL$CCA_Direction_X_2_fused_lasso_std)
    
    Direction_CCA_Vec2_fused_lasso_std_mean[,1]<-Direction_vector2_SC_std
    
    Direction_CCA_Vec2_fused_lasso_std_mean<-cbind(Direction_CCA_Vec2_fused_lasso_std_mean,gene_chr)
    
    #head( Direction_CCA_Vec2_fused_lasso_std_mean)
    
    Direction_CCA_Vec2_fused_lasso_std_mean_chr<-subset(Direction_CCA_Vec2_fused_lasso_std_mean, Direction_CCA_Vec2_fused_lasso_std_mean$gene_chr==1)
    
    Direction_CCA_Vec2_fused_lasso_std_mean_chr_significant<-which(Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1]!=0)
    
    #length(Direction_CCA_Vec2_fused_lasso_std_mean_chr_significant)
    
    Direction_CCA_Vec2_fused_lasso_std_mean_chr_1_normed<-Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1]/norm_vec(Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1])
    
     #norm_vec(Direction_CCA_Vec2_fused_lasso_std_mean_chr_normed)
    # norm_vec(Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1])
    # 
    # length(Direction_CCA_Vec2_fused_lasso_std_mean_chr_significant)
    
    # 
    ##### Direction vector summary Ordered####
    
    Direction_CCA_Vec1_fused_lasso_O_mean<-array(dim=c(p1,p1))
    
    Direction_CCA_Vec1_fused_lasso_O_mean<-Re(FL$CCA_Direction_X_1_fused_lasso_std)
    
    Direction_CCA_Vec1_fused_lasso_O_mean[,1]<-Direction_vector1_SC_O
    
    #############  Direction Vector 2 Summary ########
    
    Direction_CCA_Vec2_fused_lasso_O_mean<-array(dim=c(p2,p2))
    
    Direction_CCA_Vec2_fused_lasso_O_mean<-Re(FL$CCA_Direction_X_2_fused_lasso_O)
    
    Direction_CCA_Vec2_fused_lasso_O_mean[,1]<-Direction_vector2_SC_O
    
    Direction_CCA_Vec2_fused_lasso_O_mean<-cbind(Direction_CCA_Vec2_fused_lasso_O_mean,gene_chr)
    
    # head( Direction_CCA_Vec2_fused_lasso_O_mean)
    
    Direction_CCA_Vec2_fused_lasso_O_chr<-subset(Direction_CCA_Vec2_fused_lasso_O_mean,gene_chr==1)
    
    ##sum(Direction_CCA_Vec2_fused_lasso_O_chr[,1]^{2})
    
    Direction_CCA_Vec2_fused_lasso_O_mean_chr_1_normed<-Direction_CCA_Vec2_fused_lasso_O_chr[,1]/norm_vec(Direction_CCA_Vec2_fused_lasso_O_chr[,1])
    
    
    #sum(Direction_CCA_Vec2_fused_lasso_O_mean_chr_1_normed^{2})
    
    #length(Direction_CCA_Vec2_fused_lasso_O_mean_chr_significant)
    
    #############  Summary  ######
    
    Significant_Copy_Number_Std<-length(Direction_vector1_SC_std[Direction_vector1_SC_std!=0])
    
    Significant_Gene_Count_Std<-length(Direction_vector2_SC_std[Direction_vector2_SC_std!=0])
    
    Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_Std<-round(sum(Direction_CCA_Vec2_fused_lasso_std_mean_chr[,1]^{2}),4)*100
    
    Significant_Copy_Number_O<-length(Direction_vector1_SC_O[Direction_vector1_SC_O!=0])
    
    Significant_Gene_Count_O<-length(Direction_vector2_SC_O[Direction_vector2_SC_O!=0])
    
    Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_O<-round(sum(Direction_CCA_Vec2_fused_lasso_O_chr[,1]^{2}),4)*100
    
    Summary_Analysis<-list(CCA_FL_std=CCA_FL_std,CCA_FL_O=CCA_FL_O, Significant_Copy_Number_Std=Significant_Copy_Number_Std,Significant_Gene_Count_Std= Significant_Gene_Count_Std, Significant_Copy_Number_O= Significant_Copy_Number_O, Significant_Gene_Count_O=Significant_Gene_Count_O,Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_Std=   Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_Std,  Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_O= Percentage_Weightage_of_Significant_Genes_at_Chromosome_1_O) 
    
    return(Summary_Analysis)
    
    }
    
    
    
    
    
    
    
    
    
    
















