M <- 10
N_case <- 1000
N_control <- 1000

q <- sample(q_pool, M)
alpha <- rep(0.01, times=M)
beta <- rep(0, M)
noise <- runif(n=M, min=-1e-2, max=1e-2)
beta <- beta+noise

data_table <- table_bayes_datagen(M=M, N_case=N_case, N_control=N_control, q=q, alpha=alpha, beta=beta)


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
  T_result
}else{
  X_case <- sum(data_table["case", D_set])
  X_control <- sum(data_table["control", D_set])
  T_result <- (X_case - N_case/N_control*X_control)/sqrt(X_case + N_case^2/N_control^2*X_control)
  
  T_result
  
}

