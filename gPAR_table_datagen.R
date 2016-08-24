##########################################################################################################
# 1. qU_to_qA -- generate data table using gPAR
########################################################################################################## 
# 公式:
#		qA = r*qU/(1+(r-1)*qU)
#		r = alpha/((1-alpha)*qU) + 1
# 输入:
#		qU: MAF in the unaffected population (Wright's formular)
#		alpha: PAR
# 输出:
#		qA: MAF in the affected population
##########################################################################################################

qU_to_qA <- function(qU, alpha){
  qA <- qU*(1-alpha) + alpha
  return(qA)
}

##########################################################################################################
# 1. gPAR_table_datagen -- generate data table using gPAR
########################################################################################################## 
# 方法:
#	给定 qU, alpha, 得到 data table
#	alpha=0 neutral, alpha>0 risk, alpha<0 protective
#
# 输入:
#	qU: length M
#	alpha: length M
#	N_case:
#	N_control:
# 输出:
#	data_table
########################################################################################################## 


gPAR_table_datagen <- function(M, qU, alpha, N_case, N_control){
  
  genotypes <- matrix(0, nrow=2, ncol=M, dimnames=list( c("case", "control"), paste("G",1:M,sep="")) )
  for(j in 1:M){
    genotypes["case", j] <- rbinom(n = 1, size = N_case,  prob = qA[j])
    genotypes["control", j] <- rbinom(n = 1, size = N_case,  prob = qU[j])
  }
  
  return(genotypes)
}


# 
# gPAR <- 0.04
# M <- 50 # total number of variants
# m <- 10 # number of risk allele
# qU <- sample(q_pool, size = M)
# alpha <- c(rep(gPAR/m, m) , rep(0, M-m) )
# qA <- qU_to_qA(qU, alpha)
# gPAR_table_datagen(M = M, qU = qU, alpha = alpha, N_case = N_case, N_control = N_control)
