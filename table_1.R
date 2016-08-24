# sample size
N0 = 1000
N1 = 1000

# number of permutation
K = 1000

# repeat number of null hypothesis
K_H0 = 5000

# repeat number of alternative hypothesis
K_H1 = 1000

# allele number M
M <- seq(10, 100, by=10)

# significance level
alpha <- c(0,01, 0,05)

# threshold of ePAST
threshold <- function(C_1, N0, N1, p_hat){
  result = 1/(   N1/N0 + C_1^2/(4*N0*p_hat)  + sqrt( (N1/N0 + N1^2/N0^2)*C_1^2/(2*N0*p_hat)  + C_1^4/(16*N0^2*p_hat^2)  ) )
  return(result)
}

T_test <- function(X, Y, N1, N0){
  
  result <- (X - N1/N0*Y)/sqrt(X + N1^2/N0^2*Y)
  
  return(result)
  
}




threshold(1.95, N0 = 1000, N1 = 1000, 0.01)
threshold(1.95, N0 = 1000, N1 = 1000, 0.001)
threshold(1.95, N0 = 1000, N1 = 1000, 0.02)

C_1 = 1.65
p_hat=0.02
1/(1 + C_1^2/(4*1000*p_hat)  + sqrt( 2*C_1^2/(2*1000*p_hat)  + C_1^4/(16*1000^2*p_hat^2)  ))

### tiny try
# null hypothesis
# M = 10



# ePAST (explicit)
M <- 10
alpha <- 0.01
C_1 <- qnorm(1-alpha)
N0 = 1000
N1 = 1000

result_list <- NULL

for(k in 1:1000){
  cat(k, fill = TRUE)
  q_U <- replicate(M, r_wright(...))
  q_A <- q_U
  
  
  table_case_control <- matrix(-1, nrow = 2, ncol = M, dimnames = list(c("case", "control"), 1:10))
  
  
  for(j in 1:M){
    
    q <- q_U[j]
    p <- q
    # control
    table_case_control["control", j] <- rbinom(1, size = 2*N0, prob = q)
    # case
    table_case_control["case", j] <- rbinom(1, size = 2*N1, prob = p)
    
  }
  
  p_hat <- sum(table_case_control["control", ])/(2*N0*M)
  
  C_threshold <- threshold(C_1, N0, N1, p_hat)
  
  D_score <- table_case_control["control",] / (table_case_control["case", ] + 0.001)
  D_flag <- D_score <= C_threshold
  
  Y <- table_case_control["control", D_flag]
  X <- table_case_control["case", D_flag]
  
  
  result <- T_test(X = X,Y = Y, N1=N1, N0=N0)
  
  result_list <- append(result_list, result)
}



sum(result_list > C_1)

table_case_control["control", D_flag]/table_case_control["case", D_flag]
