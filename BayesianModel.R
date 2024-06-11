# Require Input: Y, S
n = nrow(Y)
J = ncol(Y)
X = rowSums(Y)

#--cutesv, pbsv, sniffles, debreak, svim
a.censor = c(-Inf, -Inf, -Inf, -Inf, 0)
b.censor = c(Inf, Inf, 60, 60, Inf)

dirichlet_prior = rep(1,J)
sigma_prior = c(NA, NA, 1, 1, 1)
mu_prior = c(NA, NA, 60, 60, 0)

score_set = which(colSums(is.na(S)==FALSE)!=0)

library(MCMCpack)

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
}

# Sampler:
Iter = 1000
Z_iter = matrix(NA, n, Iter); Z_iter[,1] = apply(Y, 1, function(Yi) as.numeric(sum(Yi) >= J/2))
P_iter = matrix(NA, n, Iter); P_iter[,1] = rep(0.5, n)
S_iter = array(NA, dim = c(n, J, Iter)); S_iter[,,1:Iter] = S
lambda_iter = rep(NA, Iter); lambda_iter[1] = 0.5

Omega1_iter = matrix(NA, J, Iter); Omega1_iter[,1] = seq(0.1, 0.9, length = J) / sum(seq(0.1, 0.9, length = J))
Omega0_iter = matrix(NA, J, Iter); Omega0_iter[,1] = seq(0.9, 0.1, length = J) / sum(seq(0.9, 0.1, length = J))

gamma1_iter = matrix(NA, J, Iter); gamma1_iter[,1] = Omega1_iter[,1] * sum(Z_iter[,1]==1)
gamma0_iter = matrix(NA, J, Iter); gamma0_iter[,1] = Omega0_iter[,1] * sum(Z_iter[,1]==0)

mu1_iter = matrix(NA, J, Iter); mu1_iter[,1] = rep(1, J)
mu0_iter = matrix(NA, J, Iter); mu0_iter[,1] = rep(-1, J)
sigma1_iter = matrix(NA, J, Iter); sigma1_iter[,1] = rep(1, J)
sigma0_iter = matrix(NA, J, Iter); sigma0_iter[,1] = rep(1, J)

likelihood.censor.normal <- function(y,a,b,mu,sigma){
  if(y<=a)
    return(pnorm((a-mu)/sigma))
  else if(y>=b)
    return(1-pnorm((b-mu)/sigma))
  else
    return(dnorm(y, mu, sigma))
}

for(t in 2:Iter){
  #------Update mu1, sigma1, mu0, sigma0
  for(j in score_set){
    tmp1 = Z_iter[,t-1]==1 & Y[,j]==1
    mu1_iter[j,t] = rnorm(1, (mean(S_iter[tmp1,j,t-1])*sum(tmp1)/sigma1_iter[j,t-1]^2 + mu_prior[j]/sigma_prior[j])/(sum(tmp1)/sigma1_iter[j,t-1]^2+1/sigma_prior[j]), 1/sqrt(sum(tmp1)/sigma1_iter[j,t-1]^2+1/sigma_prior[j]))
    
    tmp0 = Z_iter[,t-1]==0 & Y[,j]==1
    mu0_iter[j,t] = rnorm(1, (mean(S_iter[tmp0,j,t-1])*sum(tmp0)/sigma0_iter[j,t-1]^2 + mu_prior[j]/sigma_prior[j])/(sum(tmp0)/sigma0_iter[j,t-1]^2+1/sigma_prior[j]), 1/sqrt(sum(tmp0)/sigma0_iter[j,t-1]^2+1/sigma_prior[j]))
    
    if(mu0_iter[j,t] > mu1_iter[j,t]){
      tmp = mu0_iter[j,t]
      mu0_iter[j,t] = mu1_iter[j,t]
      mu1_iter[j,t] = tmp
    }
    
    sigma1_iter[j,t] = sqrt(1/rgamma(1, sum(tmp1)/2 + 2, sum((S_iter[tmp1,j,t-1]-mu1_iter[j,t])^2)/2 + 1))
    sigma0_iter[j,t] = sqrt(1/rgamma(1, sum(tmp0)/2 + 2, sum((S_iter[tmp0,j,t-1]-mu0_iter[j,t])^2)/2 + 1))
  }
  
  #------Update uncensored score S_iter
  for(j in score_set){
    tmp1 = Z_iter[,t-1]==1 & Y[,j]==1 & S[,j] == b.censor[j]
    S_iter[tmp1,j,t] = rtrunc(sum(tmp1), "norm", b.censor[j], Inf, mean = mu1_iter[j,t], sd = sigma1_iter[j,t])
    
    tmp0 = Z_iter[,t-1]==0 & Y[,j]==1 & S[,j] == b.censor[j]
    S_iter[tmp0,j,t] = rtrunc(sum(tmp0), "norm", b.censor[j], Inf, mean = mu0_iter[j,t], sd = sigma0_iter[j,t])
    
    tmp1 = Z_iter[,t-1]==1 & Y[,j]==1 & S[,j] == a.censor[j]
    S_iter[tmp1,j,t] = rtrunc(sum(tmp1), "norm", -Inf, a.censor[j], mean = mu1_iter[j,t], sd = sigma1_iter[j,t])
    
    tmp0 = Z_iter[,t-1]==0 & Y[,j]==1 & S[,j] == a.censor[j]
    S_iter[tmp0,j,t] = rtrunc(sum(tmp0), "norm", -Inf, a.censor[j], mean = mu0_iter[j,t], sd = sigma0_iter[j,t])
  }
  
  #------Update Omega
  tmp1 = sapply(1:J, function(k) return(sum(X == k & Z_iter[, t-1] == 1))) + dirichlet_prior
  tmp0 = sapply(1:J, function(k) return(sum(X == k & Z_iter[, t-1] == 0))) + dirichlet_prior
  
  gamma1_iter[1,t] = rtrunc(1, "gamma", -Inf, gamma1_iter[2,t-1], shape=tmp1[1], scale=1)
  if(gamma1_iter[1,t] == 0){
    gamma1_iter[1,t] = 0.01
  }else if(gamma1_iter[1,t] == Inf){
    gamma1_iter[1,t] = gamma1_iter[2,t-1]
  }
  
  for(k in 2:(J-1)){
    gamma1_iter[k,t] = rtrunc(1, "gamma", gamma1_iter[k-1,t], gamma1_iter[k+1,t-1], shape=tmp1[k], scale=1)
    if(gamma1_iter[k,t] == 0){
      gamma1_iter[k,t] = gamma1_iter[k-1,t]
    }else if(gamma1_iter[k,t] == Inf){
      gamma1_iter[k,t] = gamma1_iter[k+1,t-1]
    }
  }
  gamma1_iter[J,t] = rtrunc(1, "gamma", gamma1_iter[J-1,t], Inf, shape=tmp1[J], scale=1)
  if(gamma1_iter[J,t] == 0){
    gamma1_iter[J,t] = gamma1_iter[J-1,t]
  }
  
  Omega1_iter[,t] = gamma1_iter[,t]/sum(gamma1_iter[,t])
  
  
  gamma0_iter[1,t] = rtrunc(1, "gamma", gamma0_iter[2,t-1], Inf, shape=tmp0[1], scale=1)
  if(gamma0_iter[1,t] == 0){
    gamma0_iter[1,t] = gamma0_iter[2,t-1]
  }
  for(k in 2:(J-1)){
    gamma0_iter[k,t] = rtrunc(1, "gamma", gamma0_iter[k+1,t-1], gamma0_iter[k-1,t], shape=tmp0[k], scale=1)
    if(gamma0_iter[k,t] == 0){
      gamma0_iter[k,t] = gamma0_iter[k+1,t-1]
    }else if(gamma0_iter[k,t] == Inf){
      gamma0_iter[k,t] = gamma0_iter[k-1,t]
    }
  }
  gamma0_iter[J,t] = rtrunc(1, "gamma", -Inf, gamma0_iter[J-1,t], shape=tmp0[J], scale=1)
  if(gamma0_iter[J,t] == 0){
    gamma0_iter[J,t] = 0.01
  }else if(gamma0_iter[J,t] == Inf){
    gamma0_iter[J-1,t]
  }
  
  Omega0_iter[,t] = gamma0_iter[,t]/sum(gamma0_iter[,t])
  
  #------Update Z_i
  for(i in 1:n){
    tmp1 = prod(sapply(score_set, 
                       function(j) ifelse(Y[i,j]==1, likelihood.censor.normal(S[i,j], a.censor[j], b.censor[j], mu1_iter[j,t], sigma1_iter[j,t]), 1))) *
      Omega1_iter[X[i],t]
    
    tmp0 = prod(sapply(score_set, 
                       function(j) ifelse(Y[i,j]==1, likelihood.censor.normal(S[i,j], a.censor[j], b.censor[j], mu0_iter[j,t], sigma0_iter[j,t]), 1))) *
      Omega0_iter[X[i],t]
    
    P_iter[i, t] = lambda_iter[t-1]*prod(tmp1)/(lambda_iter[t-1]*prod(tmp1)+(1-lambda_iter[t-1])*prod(tmp0))
    Z_iter[i, t] = rbinom(1, 1, P_iter[i, t])
  }
  
  #------Update lambda
  lambda_iter[t] = rbeta(1, sum(Z_iter[,t])+1, sum(1-Z_iter[,t])+1)
  
  #progress hint
  if(0 == t%%200)
    cat(t,",") 
  
}

idx = 801:1000
P_esi = rowMeans(P_iter[,idx])

FDR <- function(tau, P_esi){
  D = as.numeric(P_esi>tau)
  return(1-mean(P_esi[D==1]))
}

x=seq(0,1,by=0.001)
y=sapply(x, function(i) FDR(i, P_esi))
plot(x,y, type = "l", main = "False Discovery Rate Control", xlab = "PEP Threshold", ylab = "FDR")
abline(h=0.05, col=2)
abline(h=0.01, col=3)
legend("topright", legend = c("FDR of 0.05","FDR of 0.01"), lty=1, col = c(2,3))


#Model-0.950
D_0.95 = P_esi>x[which.min(abs(y-0.05))]
#Model-0.990
D_0.99 = P_esi>x[which.min(abs(y-0.01))]
#Model-0.999
D_0.999 = P_esi>x[which.min(abs(y-0.001))]
