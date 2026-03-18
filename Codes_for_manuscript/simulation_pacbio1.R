## Inference for the simulated pacbio-like sequencing data

set.seed(123)
# read S
Y = read.table("./result_simulated_1/Indexd.csv", header=T, sep=',')
method = Y[,1]
Y = t(Y[,-1])

# read S
S = read.table("./result_simulated_1/Scored.csv", header=T, sep=',')
S = t(S[,-1])

# re-order
reset = sapply(c("cutesv", "pbsv", "sniffles", "debreak", "svim"), function(c) return(which(method==c)))
Y = Y[,reset]
S = S[,reset]
method = method[reset]

## --------------------------MCMC Sampler--------------------------
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

#Input: X, Y, S
n = nrow(Y)
J = ncol(Y)
X = rowSums(Y)

#--cutesv, pbsv, sniffles, debreak, svim
a.censor = c(-Inf, -Inf, -Inf, -Inf, 0)
b.censor = c(Inf, Inf, 60, 60, Inf)

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

## --------------------------Posterior Inference--------------------------
idx = 801:1000
P_esi = rowMeans(P_iter[,idx])

omega1 = rowMeans(Omega1_iter[,idx])
omega0 = rowMeans(Omega0_iter[,idx])

info = read.table("./result_simulated_1/IndexSVInfo.csv", header=T, sep=',')

FDR1 = rep(NA, J)
for(j in 1:J){
  FDR1[j] = 1-mean(P_esi[Y[,j]==1])
}
# 1-FDR1

FDR <- function(tau, P_esi){
  D = as.numeric(P_esi>tau)
  return(1-mean(P_esi[D==1]))
}

par(mfrow=c(1,2))
hist(1-P_esi, breaks = 50, main = "Histogram of Posterior Error Probabilities", xlab = "PEP", freq = F)

x=seq(0,1,by=0.001)
y=sapply(x, function(i) FDR(i, P_esi))
plot(x,y, type = "l", main = "False Discovery Rate Control", xlab = "PEP Threshold", ylab = "FDR")
abline(h=0.05, col=2)
abline(h=0.01, col=3)
#abline(h=0.001)
legend("topright", legend = c("FDR of 0.05","FDR of 0.01"), lty=1, col = c(2,3))


D_0.95 = P_esi>x[which.min(abs(y-0.05))]
res_D_0.95 = info[D_0.95,]
# write.csv(res_D_0.95, './arrayfile_simulated/res_D_0.95.csv', row.names = F)

D_0.99 = P_esi>x[which.min(abs(y-0.01))]
res_D_0.99 = info[D_0.99,]
# write.csv(res_D_0.99, './arrayfile_simulated/res_D_0.99.csv', row.names = F)

D_0.999 = P_esi>x[which.min(abs(y-0.001))]
res_D_0.999 = info[D_0.999,]
# write.csv(res_D_0.999, './arrayfile_simulated/res_D_0.999.csv', row.names = F)

D_vote1 = apply(Y, 1, function(x) sum(x)>J/2)
D_vote2 = apply(Y, 1, function(x) sum(x)==J)

res_D_vote1 = info[D_vote1,]
res_D_vote2 = info[D_vote2,]

# mean(P_esi[D_vote1])
# mean(P_esi[D_vote2])

# write.csv(res_D_vote1, './arrayfile_simulated/res_D_vote1.csv', row.names = F)
# write.csv(res_D_vote2, './arrayfile_simulated/res_D_vote2.csv', row.names = F)


## --------------------------compare with truth--------------------------
compare <- function(df1, df2, Span_loc = 500, Ratio_len = 2.0){
  # df2 is the ground truth data frame
  
  df1$lb = df1$start - Span_loc + 1
  df1$rb = df1$start + Span_loc
  df2$lb = df2$start - Span_loc + 1
  df2$rb = df2$start + Span_loc
  
  type_list = unique(c(df1$type, df2$type))
  
  count_joint1 = 0
  count_joint2 = 0
  
  for(k in 1:length(type_list)){
    df1_sub = df1[df1$type == type_list[k],]
    df2_sub = df2[df2$type == type_list[k],]
    # do not count for unique subtype
    if(nrow(df1_sub) == 0 | nrow(df2_sub) == 0){
      next
    }else{
      df1_sub = df1_sub[order(df1_sub$chro, df1_sub$start, df1_sub$length), ]
      df2_sub = df2_sub[order(df2_sub$chro, df2_sub$start, df2_sub$length), ]
      
      for(i in 1:nrow(df1_sub)){
        # Extract all the records that satisfy the position condition
        tmp_dfs = df2_sub[df2_sub$start >= df1_sub$lb[i] & df2_sub$start <= df1_sub$rb[i], ]
        if(nrow(tmp_dfs) > 0){
          # Find the record with the closest length
          tmp_idx = which.min(abs(tmp_dfs$length - df1_sub$length[i]))
          tmp_ratio = tmp_dfs$length[tmp_idx] / df1_sub$length[i]
          if(tmp_ratio <= Ratio_len & tmp_ratio >= 1/Ratio_len){
            count_joint1 = count_joint1 + 1
          }
        }
      }
      
      for(i in 1:nrow(df2_sub)){
        tmp_dfs = df1_sub[df1_sub$start >= df2_sub$lb[i] & df1_sub$start <= df2_sub$rb[i], ]
        if(nrow(tmp_dfs) > 0){
          tmp_idx = which.min(abs(tmp_dfs$length - df2_sub$length[i]))
          tmp_ratio = tmp_dfs$length[tmp_idx] / df2_sub$length[i]
          if(tmp_ratio <= Ratio_len & tmp_ratio >= 1/Ratio_len){
            count_joint2 = count_joint2 + 1
          }
        }
      }
    }
  }
  
  precision = count_joint1 / nrow(df1)
  recall = count_joint2 / nrow(df2)
  
  f1 = 2*precision*recall/(precision+recall)
  cat("number", "precision", "recall", "f1", "\n")
  cat(nrow(df1), precision, recall, f1, "\n")
  res <- c(nrow(df1), precision, recall, f1)
  return(res)
}


path = "./arrayfile_simulated/"
truth = read.table(paste(path, "truth.csv", sep = ''), header=T, sep=',')
# my0.95 = read.table(paste(path, "res_D_0.95.csv", sep = ''), header=T, sep=',')
# my0.99 = read.table(paste(path, "res_D_0.99.csv", sep = ''), header=T, sep=',')
# my0.999 = read.table(paste(path, "res_D_0.999.csv", sep = ''), header=T, sep=',')

combisv = read.table(paste(path, "combisv.csv", sep = ''), header=T, sep=',')
combisv = combisv[complete.cases(combisv),]

# D_vote1 = read.table(paste(path, "res_D_vote1.csv", sep = ''), header=T, sep=',')
# D_vote2 = read.table(paste(path, "res_D_vote2.csv", sep = ''), header=T, sep=',')

merged = read.table(paste(path, "merged_1.csv", sep = ''), header=T, sep=',')
method = c("cutesv", "pbsv", "sniffles", "debreak", "svim")
simulation_results <- NULL
for(m in method){
  simulation_results <- rbind(simulation_results,format(compare(merged[merged$Method==m,], truth),digits=3))

}

# compare(my0.95, truth)
# compare(my0.99, truth)
# compare(my0.999, truth)
simulation_results <- rbind(simulation_results,
                            format(compare(res_D_0.95, truth),digits=3),
                            format(compare(res_D_0.99, truth),digits=3),
                            format(compare(res_D_0.999, truth),digits=3),
                            format(compare(combisv, truth),digits=3),
                            format(compare(res_D_vote1, truth),digits=3),
                            format(compare(res_D_vote2, truth),digits=3))

idx <- function(df1, df2, Span_loc = 500, Ratio_len = 1.5){
  df1$original_rowid <- 1:nrow(df1)
  df2$original_rowid <- 1:nrow(df2)
  
  df1$lb = df1$start - Span_loc + 1
  df1$rb = df1$start + Span_loc
  df2$lb = df2$start - Span_loc + 1
  df2$rb = df2$start + Span_loc
  
  id1 <- integer(nrow(df1))
  id2 <- integer(nrow(df2))
  
  type_list = unique(c(df1$type, df2$type))
  
  for(k in 1:length(type_list)){
    df1_sub = df1[df1$type == type_list[k],]
    df2_sub = df2[df2$type == type_list[k],]
    # do not count for unique subtype
    if(nrow(df1_sub) == 0 | nrow(df2_sub) == 0){
      next
    }else{
      df1_sub = df1_sub[order(df1_sub$chro, df1_sub$start, as.numeric(df1_sub$length)), ]
      df2_sub = df2_sub[order(df2_sub$chro, df2_sub$start, as.numeric(df2_sub$length)), ]
      
      for(i in 1:nrow(df1_sub)){
        # Extract all the records that satisfy the position condition
        tmp_dfs = df2_sub[df2_sub$start >= df1_sub$lb[i] & df2_sub$start <= df1_sub$rb[i], ]
        if(nrow(tmp_dfs) > 0){
          # Find the record with the closest length
          tmp_idx = which.min(abs(as.numeric(tmp_dfs$length) - as.numeric(df1_sub$length)[i]))
          tmp_ratio = as.numeric(tmp_dfs$length)[tmp_idx] / as.numeric(df1_sub$length)[i]
          if(tmp_ratio <= Ratio_len & tmp_ratio >= 1/Ratio_len){
            id1[df1_sub$original_rowid[i]] <- 1
          }else{
            id1[df1_sub$original_rowid[i]] <- 0
          }
        }else{
          id1[df1_sub$original_rowid[i]] <- 0
        }
      }
      
      for(i in 1:nrow(df2_sub)){
        tmp_dfs = df1_sub[df1_sub$start >= df2_sub$lb[i] & df1_sub$start <= df2_sub$rb[i], ]
        if(nrow(tmp_dfs) > 0){
          tmp_idx = which.min(abs(as.numeric(tmp_dfs$length) - as.numeric(df2_sub$length)[i]))
          tmp_ratio = as.numeric(tmp_dfs$length)[tmp_idx] / as.numeric(df2_sub$length)[i]
          if(tmp_ratio <= Ratio_len & tmp_ratio >= 1/Ratio_len){
            id2[df2_sub$original_rowid[i]] <- 1 
          }else{
            id2[df2_sub$original_rowid[i]] <- 0
          }
        }else{
          id2[df2_sub$original_rowid[i]] <- 0
        }
      }
    }
  }
  return(list(id1, id2))
}

# the value
rr_combisv = idx(info, combisv, Span_loc = 500, Ratio_len = 2.0)


fdr_pacbio1 <- c(1-FDR1, 0.95,0.99,0.999,mean(P_esi[rr_combisv[[1]]==1]),mean(P_esi[D_vote1]), mean(P_esi[D_vote2]))
#format(compare(res_D_0.95, truth),digits=3)
#format(compare(res_D_0.99, truth),digits=3)
#format(compare(res_D_0.999, truth),digits=3)

#format(compare(combisv, truth),digits=3)
# compare(D_vote1, truth)
# compare(D_vote2, truth)
#format(compare(res_D_vote1, truth),digits=3)
#format(compare(res_D_vote2, truth),digits=3)


row_names <- c("cuteSV","pbsv","Sniffles","DeBreak","SVIM","Model-0.950","Model-0.990","Model-0.999","combiSV","Vote 1","Vote 2")
col_names <- c("number of SVs","1-FDR", "Precision", "Recall", "F1 score")  

simulation_results <- cbind(simulation_results[,1],fdr_pacbio1,simulation_results[,2:4])

data <- data.frame(simulation_results)
colnames(data) <- col_names  
rownames(data) <- row_names

write.csv(data, file = "table4.csv")
