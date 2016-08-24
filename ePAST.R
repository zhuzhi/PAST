##########################################################################################################
# 1. ePAST -- ePAST
########################################################################################################## 
#
# 方法：
# 		针对每个位点计算 T_j， 通过阈值 C 来选取 risk set D={j| T_j >= C_1 }.
#		等价的，通过解方程可得到方程:
#			X_case_j/X_control_j <= 1/(  
#										N_case/N_control + C_1^2/(2*X_control_j) + 
#										sqrt(  
#												(N_case/N_control + N_case^2/N_control^2)*C_1^2/X_control_j + 
#												C_1^4 /(4*X_control_j^2 )
#											)      
#									  )
#
# 公式:
#		T_j = (X_case_j - N_case/N_control*X_control_j)/sqrt(X_case_j + N_case^2/N_control^2*X_control_j) (*1)
#		C_2 = 1/(  
#				 N_case/N_control + C_1^2/(4*N_control*p_hat) + 
#				 sqrt(  
#						(N_case/N_control + N_case^2/N_control^2)*C_1^2/(2*N_control*p_hat) + 
#						 C_1^4 /(16*N_control^2*p_hat^2)   
#					 )      
#				)
#		T_D = (X_case - N_Case/N_control*X_control)/sqrt(X_case + N_Case^2/N_control^2*X_control)
#
# 其中： D is the risk allele set determined by C.
#		p_hat 是对照组中 allele freq 的一个总体估计， 2*N_control*p_hat 是 X_control_j 的一个估计
# 输入: 
# 		data_table: case/control table, row names 为 "case", "control"
#		alpha_level: significance level, 等价于 给出了 C_1
#		N_case: case组的总人数
#		N_control: control组的总人数
# 输出:
#		ePAST结果
# 流程:
# 		1. sorting: 决定各个 D_k
#		2. 由公式计算各个 T_k，并且得到最终 T_k
# 注意事项:
#		1. p_hat 理论上可能为 0，但这里不予考虑，几乎不可能发生。p_hat是先验的，也可以是计算的？
#		2. p_hat 去一个固定值，则 X_control_j/(2*N_control)<p_hat 的位点，被选取的概率增大，反之则减小。
#			即一个固定的 p 有增加 稀有基因 被挑选的概率的好处。
#		3. X_control_j/X_case_j 有可能出现 X_case_j 为0的情形，此时，直接视为无穷大，排除之。
##########################################################################################################



ePAST <- function(data_table, N_case, N_control, alpha_level=0.05){

  C_1 <- qnorm(p = 1-alpha_level)
  M <- ncol(data_table)
#   C_2 <- 1/(  
#     N_case/N_control + C_1^2/(4*N_control*p_hat) + 
#       sqrt(  
#         (N_case/N_control + N_case^2/N_control^2)*C_1^2/(2*N_control*p_hat) + 
#           C_1^4 /(16*N_control^2*p_hat^2)   
#       )      
#   )
  T_result <- numeric(M)
  D_set <- numeric(0)
  for(j in 1:M){
    
    X_case <- sum(data_table["case", j])
    X_control <- sum(data_table["control", j])
    if(X_case==0 & X_control==0){
      next
    }
    T_result[j] <- (X_case - N_case/N_control*X_control)/sqrt(X_case + N_case^2/N_control^2*X_control)
    if(T_result[j] >= C_1){
      D_set <- c(D_set, j)
    }
  }
  
  if(length(D_set)==0){
    
    T_result <- 0
    return(T_result)
  }else{
    X_case <- sum(data_table["case", D_set])
    X_control <- sum(data_table["control", D_set])
    T_result <- (X_case - N_case/N_control*X_control)/sqrt(X_case + N_case^2/N_control^2*X_control)
    
    return(T_result)
    
  }

}







# 
# 预先给定的p并不能控制 type-I error。
# 
# 
# ePAST <- function(data_table, N_case, N_control, alpha_level=0.05, p_hat=0.01){
#   
#   C_1 <- qnorm(p = 1-alpha_level)
#   M <- ncol(data_table)
#   C_2 <- 1/(  
#     N_case/N_control + C_1^2/(4*N_control*p_hat) + 
#       sqrt(  
#         (N_case/N_control + N_case^2/N_control^2)*C_1^2/(2*N_control*p_hat) + 
#           C_1^4 /(16*N_control^2*p_hat^2)   
#       )      
#   )
#   
#   D_set <- numeric(0)
#   
#   for(j in 1:M){
#     
#     X_case <- data_table["case", j]
#     X_control <- data_table["control", j]
#     
#     if(X_case==0){
#       next
#     }
#     
#     if(X_control/X_case<=C_2){
#       D_set <- c(D_set, j)
#     }
#     
#   }
#   
#   if(length(D_set)==0){
#     
#     T_result <- 0
#     return(T_result)
#   }
#   
#   
#   X_case <- sum(data_table["case", D_set])
#   X_control <- sum(data_table["control", D_set])
#   T_result <- (X_case - N_case/N_control*X_control)/sqrt(X_case + N_case^2/N_control^2*X_control)
#   
#   return(T_result)
#   
# }






##########################################################################################################
# setwd("/home/zhuzhi/projects/PAST/")
# source("table_datagen.R")
# source("Wright's formula.R")
# q_pool <- replicate(10000, r_wright())
# 
# M <- 10
# N_case <- 1000
# N_control <- 1000
# 
# q <- sample(q_pool, M)
# alpha <- rep(0.01, times=M)
# beta <- rep(0, M)
# noise <- runif(n=M, min=-1e-2, max=1e-2)
# beta <- beta+noise
# beta <- rep(0, M)
# 
# N <- 1000
# T_results <- numeric(N)
# for(k in 1:N){
#   data_table <- table_bayes_datagen(M=M, N_case=N_case, N_control=N_control, q=q, alpha=alpha, beta=beta)
#   T_results[k] <- ePAST(data_table=data_table, N_case=N_case, N_control=N_control)
# }
# 
# sum(T_results>qnorm(0.95), na.rm = TRUE)/N
# mean(T_results>qnorm(0.95))



