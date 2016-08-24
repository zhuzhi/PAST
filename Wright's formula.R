# setwd("/home/zhuzhi/projects/PAST/")
# Wright's formula:
# f(p| 0<p<1) = c*p^(beta_s-1)(1-p)^(beta_N-1)exp(s(1-p))
# PAST choose:
# 1. beta_s=0.001, beta_N=beta_s/3, s=12
# 2. disclard those p smaller than 1e-3. e.g. f(p|1e-3<p<1)

# f( p | 1e-3<p<1 )
#                   = c_new * p^(0.001-1) * (1-p)^(0.001/3-1) * exp(12(1-p))

# generate function by Hongyan Fang :)
# 感谢师姐
generate <- function(a){
  p = 0
  while(p == 0){
    ky = rbeta(1, 0.001, 0.001/3)
    u = runif(1)
    if(ky>=0.001 & ky <= 0.5 & u<=exp(-12*ky)){
      p = ky
    }
  }
  return(p)
}


# rejection sampling
# truncated g(p) = c_beta * p^(shape1-1) * (1-p)^(shape2-1)
# truncated f(p) = c * exp(-s*p) * p^(shape1-1) * (1-p)^(shape2-1)
#                = c_new * exp(-s*p) * g(p)
# f(p)/g(p) = c_new * exp(-s*p) <= c_new * exp(-s*left_truncation)
# let M equal c_new * exp(-s*left_truncation), so the threshold should be 
#         f(p)/(M*g(p)) = exp(-s * (p-left_truncation) )

r_wright <- function(..., left_trunction=0.001, right_trunction=0.05, 
                     beta_s=0.001, beta_N=0.001/3, s=12){
  
  flag = TRUE
  while(flag){
    x <- rbeta(1, shape1=beta_s, shape2=beta_N)
  
    if(x <= right_trunction & x >= left_trunction){
      u_x <- exp(-s*(x - left_trunction))
      flag = runif(1) > u_x
      }
  }
  return(x)
}

 