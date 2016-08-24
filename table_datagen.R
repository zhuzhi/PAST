##########################################################################################################
# 1. 产生某个位点的case/control数据 -- locus_datagen 
########################################################################################################## 
#
# 公式依据：
# 		P(Y=1) = exp(alpha+beta*G) / (1 + exp(alpha+beta*G))  	(*1)
# 其中： G为基因型，alpha为预设参数，beta代表基因型对于基因的影响
# 		beta越大，越导致疾病，beta为负的时候，可以看成是保护基因。
#
# 输入: 
# 		N：样本总数
#		alpha: 调节参数
#		beta: 基因参数
# 输出:
#		case/control样本数
# 流程:
# 	个人：	q给定 --> 由q计算G --> 由(*1)确定致病概率 --> 随机归属case/control
#		固定数量的case/control: 不断重复个人流程，直到达到固定数量的case/control。
##########################################################################################################
inner_loc_datagen <-function(q, alpha, beta){
  G <- rbinom(n=1, size=2, prob=q)
  case_probs <- exp(alpha+beta*G)/(1+exp(alpha+beta*G))
  Y=as.integer(runif(1)<=case_probs)
  return(c(Y=Y,G=G))
}
# tiny example
# q=0.1
# alpha=0.01
# beta=1
# N=10000
# tiny_sample <- replicate(10000, inner_loc_datagen(q, alpha, beta))
# by(tiny_sample["G",], INDICES=tiny_sample["Y",], FUN=table)

locus_datagen <- function(N_case, N_control, q, alpha, beta, detail=FALSE){
  
  N <- N_case+N_control
    
  genotypes_case <- NULL
  genotypes_control <- NULL
  
  n_case <- 0
  n_control <- 0
  while(n_case < N_case| n_control < N_control){
    
    genotypes <- rbinom(n=N, size=2, prob=q)
    case_prob <- exp(alpha+beta*genotypes) / (1 + exp(alpha+beta*genotypes)) 
    Y <- as.integer(runif(N)<=case_prob)
    genotypes_case <- c(genotypes_case, genotypes[Y==1])
    genotypes_control <- c(genotypes_control, genotypes[Y==0])
    n_case <- length(genotypes_case)
    n_control <- length(genotypes_control)
  }
  
  genotypes_case <- genotypes_case[1:N_case]
  genotypes_control <- genotypes_control[1:N_control]
  
  if(detail){
    return(list(genotypes_case=genotypes_case, genotypes_control=genotypes_control))
  }else{
    return(c(G_case=sum(genotypes_case), G_control=sum(genotypes_control)))
  }
  
}

# q=0.1
# alpha=0.01
# beta=1
# N=1000
# locus_datagen(N_case=N, N_control=N, q=q, alpha=alpha, beta=beta)
# tiny_sample <- locus_datagen(N_case=N, N_control=N, q=q, alpha=alpha, beta=beta, TRUE)
# table(tiny_sample$genotypes_case)
# table(tiny_sample$genotypes_control)



##########################################################################################################
# 2. bayes方法产生某个位点的case/control数据 -- locus_bayes_datagen 
########################################################################################################## 
#
# 公式依据：
# 		P(Y=1) = exp(alpha+beta*G) / (1 + exp(alpha+beta*G))  	(*1)
# 其中： G为基因型，alpha为预设参数，beta代表基因型对于基因的影响
# 		beta越大，越导致疾病，beta为负的时候，可以看成是保护基因。
#
# 输入: 
# 		N：样本总数
#		alpha: 调节参数
#		beta: 基因参数
# 输出:
#		case/control样本数
# 流程:
#		P(G=g|Y=1) = P(Y=1|G=g)*P(G=g) / {P(Y=1|G=0)*P(G=0) + P(Y=1|G=1)*P(G=1) + P(Y=1|G=2)*P(G=2)}
#		即case组是由 P(G=g|Y=1), g=c(1,2,3)给定的人multinomial分布。
# 		同理可得control组的数据。
# 		由此可得数据
##########################################################################################################

locus_bayes_datagen <- function(N_case, N_control, q, alpha, beta, detail=FALSE){
  
  G <- 0:2
  P_G <- c((1-q)^2, 2*q*(1-q), q^2)
  P_case_if_G <- exp(alpha+beta*G) / (1 + exp(alpha+beta*G))  
  P_G_if_case <- P_case_if_G*P_G/sum(P_case_if_G*P_G)
  
  P_control_if_G <- 1 - P_case_if_G
  P_G_if_control <- P_control_if_G*P_G/sum(P_control_if_G*P_G)
  
  if(detail){
    
    genotypes_case <- rmultinom(n=N_case, size=1, prob=P_G_if_case)
    genotypes_case <- genotypes_case[2, ] + 2*genotypes_case[3, ]
    
    genotypes_control <- rmultinom(n=N_control, size=1, prob=P_G_if_control)
    genotypes_control <- genotypes_control[2, ] + 2*genotypes_control[3,]
    
    return(list(genotypes_case=genotypes_case, genotypes_control=genotypes_control))
  }else{
    
    G_case <- drop( rmultinom(n=1, size=N_case, prob=P_G_if_case) )
    G_case <- sum(G_case*G)
    G_control <- drop(rmultinom(n=1, size=N_control, prob=P_G_if_control)  )
    G_control <- sum(G_control*G)

    return(c(G_case=G_case, G_control=G_control) )
  }
}

# q=0.1
# alpha=0.01
# beta=1
# N=1000
# locus_bayes_datagen(N_case=N, N_control=N, q=q, alpha=alpha, beta=beta)
# tiny_sample <- locus_bayes_datagen(N_case=N, N_control=N, q=q, alpha=alpha, beta=beta, TRUE)
# table(tiny_sample$genotypes_case)
# table(tiny_sample$genotypes_control)

# compare 
# 重复10000次，原方法需要 9.109561秒， Bayes方法需要 0.2520821。
# 两次结果均值误差在小数点后两位
# q=0.1
# alpha=0.01
# beta=1
# N=1000
# start_time <- Sys.time()
# naive_result <- replicate(10000, locus_datagen(N_case=N, N_control=N, q=q, alpha=alpha, beta=beta))
# end_time <- Sys.time()
# naive_time <- end_time-start_time
# start_time <- Sys.time()
# bayes_result <- replicate(10000, locus_bayes_datagen(N_case=N, N_control=N, q=q, alpha=alpha, beta=beta))
# end_time <- Sys.time()
# bayes_time <- end_time-start_time
# cat("Naive time: ", naive_time, ", Bayes time: ", bayes_time, "." , sep = "")
# cat("Naive mean: ", apply(naive_result,1,mean), ", Bayes result: ", apply(bayes_result,1,mean), " .")


##########################################################################################################
# 3. bayes方法产生case/control数据table -- table_bayes_datagen
########################################################################################################## 
#
# 方法depend：
# 		locus_bayes_datagen  (bayes方法产生某个位点的case/control数据，由于bayes方法显著好于naive的方法)
#
# 输入: 
# 		N_case： case组人数
#		N_control: 控制组人数
#		M: 位点个数
#		q: MAF, 向量，长度为M
#		alpha: 调节参数，向量，长度为M
#		beta: 基因参数，向量，长度为M
# 输出:
#		case/control table (Matrix, 2xM)
# 流程:
#		M 个位点，每个应用下locus_bayes_datagen
##########################################################################################################

table_bayes_datagen <- function(M, N_case, N_control, q, alpha, beta){
  
  genotypes <- matrix(0, nrow=2, ncol=M, dimnames=list( c("case", "control"), paste("G",1:M,sep="")) )
  for(j in 1:M){
    
    locus_result <- locus_bayes_datagen(N_case=N_case, N_control=N_control, 
                                        q=q[j], alpha=alpha[j], beta=beta[j])
    genotypes["case", j] <- locus_result["G_case"]
    genotypes["control", j] <- locus_result["G_control"]
  }
  
  return(genotypes)
  
}


# tiny example
# M <- 10
# N_case <- 1000
# N_control <- 1000
# q <- replicate(M, r_wright())
# alpha <- rep(0.01, times=M)
# beta <- runif(n=M, min=-1, max=1)
# beta <- sort(beta)
# 
# table_bayes_datagen(M=M, N_case=N_case, N_control=N_control, q=q, alpha=alpha, beta=beta)

