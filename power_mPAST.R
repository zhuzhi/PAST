source(file = "gPAR_table_datagen.R")
source(file = "mPAST.R")
source(file = "Wright's formula.R")

file_name = "rare_variants_only"

q_pool <- replicate(10000, r_wright())


gPAR_list <- seq(0, 0.1, by=0.01)
power_list <- NULL

K = 1000 
J = 1000
M <- 100 # total number of variants
m <- 10 # number of risk allele
N_case = 1000
N_control = 1000



for(gPAR in gPAR_list){
  
  power_results <- numeric(K)
  q_pool_e2 <- q_pool[q_pool<1e-2]
  
  for(i in 1:K){
    
    qU <- sample(q_pool_e2, size = M)
    alpha <- c(rep(gPAR/m, m) , rep(0, M-m) )
    qA <- qU_to_qA(qU, alpha)
    
    T_results <- numeric(J)
    for(j in 1:J){
      cat(i, j , sep = ",", fill = TRUE)
      data_table <- gPAR_table_datagen(M = M, qU = qU, alpha = alpha, N_case = N_case, N_control = N_control)
      T_results[j] <- mPAST(data_table = data_table, N_case =  N_case, N_control = N_control)
    }
    power_results[i] <- mean(T_results > 1.65)
    
  }

  power_list <- c(power_list, mean(power_results) )
}

write.table(x = data.frame(gPAR=gPAR_list, power = power_list, stringsAsFactors = FALSE),
            file = file_name,
            sep=",",
            row.names = FALSE
            )









