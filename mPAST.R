##########################################################################################################
# 1. mPAST -- mPAST
########################################################################################################## 
#
# 方法：
# 		T_max is taken to be the maximal value of T_k(T_D_k)
# 公式:
#		T_D = (X_case - N_case/N_control*X_control)/sqrt(X_case + N_case^2/N_control^2*X_control)
# 其中： D_k is the set containing the variants with the k smallest risk ratios.

# 输入: 
# 		data_table: case/control table, row names 为 "case", "control"
#		N_case: case组的总人数
#		N_control: control组的总人数
# 输出:
#		mPAST结果
# 流程:
# 		1. sorting: 决定各个 D_k
#		2. 由公式计算各个 T_k，并且得到最终 T_k
##########################################################################################################
mPAST <- function(data_table, N_case, N_control){
  M <- ncol(data_table)
  
  risk_ratio <- numeric(M)
  for(j in 1:M){
    if(data_table["control", j] == 0){
      risk_ratio[j] <- Inf
    }else{
      risk_ratio[j] <- (data_table["case", j]/N_case)/(data_table["control", j]/N_control) 
    }
  }
  
  data_table_reorder <- data_table[, order(risk_ratio, decreasing = FALSE)]
  
  T_result <- numeric(M)
  for(j in 1:M){
    
    X_case <- sum(data_table_reorder["case", 1:j])
    X_control <- sum(data_table_reorder["control", 1:j])
    T_result[j] <- (X_case - N_case/N_control*X_control)/sqrt(X_case + N_case^2/N_control^2*X_control)
    
  }
  
  T_max <- max(T_result)
  
  return(T_max)
}




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
# 
# 
# N <- 1000
# T_results <- numeric(N)
# for(k in 1:N){
#   data_table <- table_bayes_datagen(M=M, N_case=N_case, N_control=N_control, q=q, alpha=alpha, beta=beta)
#   T_results[k] <- mPAST(data_table=data_table, N_case=N_case, N_control=N_control)
# }
# 
# mean(T_results>qnorm(0.95))
# mean(T_results>qnorm(0.99))
# 
# 
# 
# 
# y <- 10
# plot( function(x) (x-y)/sqrt(x+y) , xlim = c(-y,2*y))
# 
# 
# x <- 10
# plot( function(y) (x-y)/sqrt(x+y) , xlim = c(-x,2*x))
# 
# 
# 
