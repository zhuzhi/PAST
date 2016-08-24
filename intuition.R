 
++z # 1. q_A and q_U relationship

q_A <- function(r, q_U= 0.5){
  r*q_U/(1+ (r-1)*q_U)
}

q_A_div_U <- function(r, q_U= 0.5){
  r/(1+ (r-1)*q_U)
}

plot(q_A, xlim = c(0, 2))
plot(q_A_div_U, xlim = c(0, 2))
