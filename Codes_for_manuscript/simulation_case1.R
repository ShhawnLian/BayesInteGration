library(MCMCpack)

set.seed(123)
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

FDR <- function(tau, P_esi){
  D = as.numeric(P_esi>tau)
  return(1-mean(P_esi[D==1]))
}

sim1 <- function(){
  J = 5
  n = 2000
  lambda = 0.70
  
  # increasin ratio
  Omega1 = c(0.05, 0.10, 0.15, 0.25, 0.45)
  Omega0 = c(0.40, 0.25, 0.15, 0.12, 0.08)
  
  Z = rbinom(n, 1, lambda)
  X = rep(NA, n)
  for(i in 1:n){
    if(Z[i]==1)
      X[i] = sample(1:J, 1, prob = Omega1)
    else
      X[i] = sample(1:J, 1, prob = Omega0)
  }
  
  Y = matrix(0, n, J)
  tmp = sapply(X, function(x) sample(1:J, x))
  for(i in 1:n){
    Y[i,tmp[[i]]]=1
  }
  
  mu1 = c(0 ,0, 2, 2, 2)
  mu0 = c(0, 0, -2, -2, -2)
  sigma1 = c(1, 1, 1 ,1 ,1)
  sigma0 = c(1, 1, 1, 1, 1)
  
  a.censor = c(-Inf, -Inf, -Inf, 0, -Inf)
  b.censor = c(Inf, Inf, Inf, Inf, 0)
  
  score_set = 3:5
  
  #generate the scores:
  S = matrix(NA, n, J)
  for(i in 1:n){
    for(j in score_set){
      if(Y[i,j]==1 & Z[i]==1){
        S[i,j] = rnorm(1, mu1[j], sigma1[j])
        if(S[i,j] > b.censor[j])
          S[i,j] = b.censor[j]
        else if(S[i,j] < a.censor[j])
          S[i,j] = a.censor[j]
      }else if(Y[i,j]==1 & Z[i]==0){
        S[i,j] = rnorm(1, mu0[j], sigma0[j])
        if(S[i,j] > b.censor[j])
          S[i,j] = b.censor[j]
        else if(S[i,j] < a.censor[j])
          S[i,j] = a.censor[j]
      }else{
        S[i,j] = NA
      }
    }
  }
  
  return(list(Y=Y, S=S, Z=Z, X=X, lambda=lambda, Omega0=Omega0, Omega1=Omega1, mu1=mu1, mu0=mu0, sigma1=sigma1, sigma0=sigma0))
}


fdrj = matrix(NA, 10, 10)
precision = matrix(NA, 10, 10)
recall = matrix(NA, 10, 10)
f1 = matrix(NA, 10, 10)

for(repl in 1:10){
  res = sim1()
  Y=res$Y
  S=res$S
  
  #Input: X, Y, S
  n = nrow(Y)
  J = ncol(Y)
  X = rowSums(Y)
  
  a.censor = c(-Inf, -Inf, -Inf, 0, -Inf)
  b.censor = c(Inf, Inf, Inf, Inf, 0)
  
  #score_set = 3:5
  score_set = which(colSums(is.na(S)==FALSE)!=0)
  
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
  
  dirichlet_prior = rep(1,J)
  
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
      mu1_iter[j,t] = rnorm(1, mean(S_iter[tmp1,j,t-1]), sigma1_iter[j,t-1]/sqrt(sum(tmp1)))
      
      tmp0 = Z_iter[,t-1]==0 & Y[,j]==1
      mu0_iter[j,t] = rnorm(1, mean(S_iter[tmp0,j,t-1]), sigma0_iter[j,t-1]/sqrt(sum(tmp0)))
      
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
    }#else if(gamma1_iter[J,t] == Inf){
    #  gamma1_iter[J,t] = n
    #}
    
    Omega1_iter[,t] = gamma1_iter[,t]/sum(gamma1_iter[,t])
    
    
    gamma0_iter[1,t] = rtrunc(1, "gamma", gamma0_iter[2,t-1], Inf, shape=tmp0[1], scale=1)
    if(gamma0_iter[1,t] == 0){
      gamma0_iter[1,t] = gamma0_iter[2,t-1]
    }#else if(gamma0_iter[1,t] == Inf){
    #  
    #}
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
  
  idx = 501:1000
  P_esi = rowMeans(P_iter[,idx])
  
  hist(P_esi, breaks = 50, main = "Posterior Mean of P(z=1|.)", xlab = "P", freq = F)
  
  x=seq(0,1,by=0.001)
  y=sapply(x, function(i) FDR(i, P_esi))
  plot(x,y, type = "l", main = "Posterior predictive FDR", xlab = "Threshold (tau)", ylab = "FDR")
  abline(h=0.05)
  abline(h=0.01)
  abline(h=0.001)
  
  D_0.95 = P_esi>x[which.min(abs(y-0.05))]
  res_D_0.95 = as.numeric(D_0.95)
  
  D_0.99 = P_esi>x[which.min(abs(y-0.01))]
  res_D_0.99 = as.numeric(D_0.99)
  
  D_0.999 = P_esi>x[which.min(abs(y-0.001))]
  res_D_0.999 = as.numeric(D_0.999)
  
  D_vote1 = apply(Y, 1, function(x) as.numeric(sum(x)>J/2))
  D_vote2 = apply(Y, 1, function(x) as.numeric(sum(x)==J))
  
  true = res$Z
  
  FDR1 = rep(NA, J)
  for(j in 1:J){
    FDR1[j] = 1-mean(P_esi[Y[,j]==1])
  }
  
  fdrj[1:5,repl] = FDR1
  fdrj[6:8,repl] = c(0.05, 0.01, 0.001)
  fdrj[9, repl] = 1-mean(P_esi[D_vote1==1])
  fdrj[10, repl] = 1-mean(P_esi[D_vote2==1])
  
  for(j in 1:5){
    precision[j,repl] = sum(true & Y[,j])/sum(Y[,j])
    recall[j,repl] = sum(true & Y[,j])/sum(true)
  }
  precision[6,repl] = sum(true & res_D_0.95)/sum(res_D_0.95)
  recall[6,repl] = sum(true & res_D_0.95)/sum(true)
  precision[7,repl] = sum(true & res_D_0.99)/sum(res_D_0.99)
  recall[7,repl] = sum(true & res_D_0.99)/sum(true)
  precision[8,repl] = sum(true & res_D_0.999)/sum(res_D_0.999)
  recall[8,repl] = sum(true & res_D_0.999)/sum(true)
  
  precision[9,repl] = sum(true & D_vote1)/sum(D_vote1)
  recall[9,repl] = sum(true & D_vote1)/sum(true)
  precision[10,repl] = sum(true & D_vote2)/sum(D_vote2)
  recall[10,repl] = sum(true & D_vote2)/sum(true)
  
  for(j in 1:10){
    f1[j,repl] = 2*precision[j,repl]*recall[j,repl] / (precision[j,repl]+recall[j,repl])
  }
  
}


# Output order: Tool 1-Tool 5, Model-0.950, Model-0.990, Model-0.999, Vote 1, Vote 2
cat("1. Estimated FDR (mean and sd):\n")
cat(apply(1-fdrj, 1, mean), "\n", apply(1-fdrj, 1, sd), '\n')
cat("2. Precision (mean and sd):\n")
cat(apply(precision, 1, mean), "\n", apply(precision, 1, sd), '\n')
cat("3. Recall (mean and sd):\n")
cat(apply(recall, 1, mean), "\n", apply(recall, 1, sd), '\n')
cat("4. F1 score (mean and sd):\n")
cat(apply(f1, 1, mean), "\n", apply(f1, 1, sd))


row_names <- c("Tool1","Tool2","Tool3","Tool4","Tool5","Model-0.950","Model-0.990","Model-0.999","Vote 1","Vote 2")
col_names <- c("1-FDR", "Precision", "Recall", "F1 score")

simulation1_results <- 
cbind(
paste(format(apply(1-fdrj, 1, mean),digits=3), "(", format(apply(1-fdrj, 1, sd),nsmall = 3,digits=1), ")", sep = ""),
paste(format(apply(precision, 1, mean),digits=3), "(", format(apply(precision, 1, sd),nsmall = 3,digits=1), ")", sep = ""),
paste(format(apply(recall, 1, mean),digits=3), "(", format(apply(recall, 1, sd),nsmall = 3,digits=1), ")", sep = ""),
paste(format(apply(f1, 1, mean),digits=3), "(", format(apply(f1, 1, sd),nsmall = 3,digits=1), ")", sep = ""))


data <- data.frame(simulation1_results)
colnames(data) <- col_names 
rownames(data) <- row_names

write.csv(data, file = "table1.csv")

