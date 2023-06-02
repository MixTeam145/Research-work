#source("mcssa_utils_zenodo.R")
library(foreach)
library(doParallel)
library(doRNG)
source("toeplitz_mssa.R")
source("mcmssa_utils.R")

N <- 100
L <- 30
G <- 1000
M <- 1000
omega <- 0.1
signal <- sin(2*pi*(1:N) * omega)
plot(signal, type = "l")

cores <- detectCores()
cluster <- makeCluster(cores - 1)
registerDoSNOW(cluster)

rejected <- foreach(i=1:M, .export=c('ssa', 'nu', 'parestimate', 'Norm','rowQuantiles'), .combine='+', .options.snow = opts) %dopar% {
  f <- generate(model, 0)
  res.ev <- MonteCarloSSA(f = f, L = L, basis = "ev", model = model, G = G, level.conf = levels.conf)
  res.ev.t <- MonteCarloSSA(f = f, L = L, basis = "ev.t", model = model, G = G, level.conf = levels.conf)
  res.t <- MonteCarloSSA(f = f, L = L, basis = "t", model = model, G = G, level.conf = levels.conf)
  rbind(res.ev$reject, res.ev.t$reject, res.t$reject)
}
ecdf <- list(ev = rejected[1,]/M, ev.t = rejected[2,]/M, t = rejected[3,]/M)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)
pvals <- foreach(i=1:100, .export = c('ssa', 'nu', 'parestimate', 'Norm', 'rowQuantiles', 'adiag', 'Toeplitz', 'toeplitz.mssa'), .combine='c') %dorng% {
  f <- generate(model, signal, D)
  (res <- MonteCarloSSA(f = f, L = 20, model = model, basis = "ev", kind = "fa", D = D, G = G, level.conf = NULL))
  res$p.value
}
plot(punif)
lines(ecdf(pvals), col = "green")
stopCluster(cl)




varphi <- 0.7
delta <- 1
omega <- 0.075
N <- 100
Ls <- c(10, 20, 50, 80, 90)
L_idx <- 1:length(Ls)
D <- 2 
G <- 1000
M <- 1000
model <- list(list(varphi = varphi,
                   delta = delta,
                   N = N),
              list(varphi = varphi,
                   delta = delta,
                   N = N))
signal <- replicate(D, 0, simplify=F)

set.seed(5)
p.values_noise.me1block.ev <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "ev", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1block.ev[[idx]] <- pvals
}


p.values_noise.me1block.fa <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "fa", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1block.fa[[idx]] <- pvals
}


pval.noise.me1b.fa <- 
  foreach (i=1:10000, .combine = 'c') %do% {
  f <- generate(model, signal, D)
  res <- MonteCarloSSA(f = f, L = 10, model = model, basis = "ev", kind = "fa", D = D, G = G, level.conf = NULL)
  res$p.value
}
#p.values_noise.me1block.fa[[idx]] <- pvals

signal <- list(signal.one.channel(model[[1]]$N, omega), 
               signal.one.channel(model[[1]]$N, omega))
pval.signal.me1b.fa <- 
  foreach (i=1:10000, .combine = 'c') %do% {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = 10, model = model, basis = "ev", kind = "fa", D = D, G = G, level.conf = NULL)
    res$p.value
  }



signal <- list(signal.one.channel(model[[1]]$N, omega), 
               signal.one.channel(model[[1]]$N, omega))
p.values_signal.me1block.ev <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "ev", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_signal.me1block.ev[[idx]] <- pvals
}

p.values_signal.me1block.fa <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "fa", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_signal.me1block.fa[[idx]] <- pvals
}



signal <- replicate(D, 0, simplify=F)

p.values_noise.me1sum.ev <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "ev", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1sum.ev[[idx]] <- pvals
}


p.values_noise.me1sum.fa <- list()
for (idx in L_idx) {
  pvals <- c()
  for (i in 1:M) {
    f <- generate(model, signal, D)
    res <- MonteCarloSSA(f = f, L = Ls[idx], model = model, basis = "ev", kind = "fa", D = D, G = G, level.conf = NULL)
    pvals <- append(pvals, res$p.value)
  }
  p.values_noise.me1sum.fa[[idx]] <- pvals
}




alphas <- 0:1000/1000
alphas_idx <- 1:length(alphas)
roc.me1sum.fa <- list()
alpha_1.me1sum.fa <- list()
beta.me1sum.fa <- list()


alphaI <- sapply(alphas, function(a) sum(pval.noise.me1b.fa < a)/10000)
beta <- sapply(alphas, function(a) sum(pval.signal.me1b.fa < a)/10000)
for (l in L_idx) {
  alphaI[[l]] <- sapply(alphas, function(a) sum(pval.noise.me1b.fa[[l]] < a)/M)
  beta[[l]] <- sapply(alphas, function(a) sum(pval.signal.me1b.fa[[l]] < a)/M)
}

plot(alphaI[[1]], beta[[1]], type = "l")
lines(alphaI[[2]], beta[[2]], col = "red")
lines(alphaI[[3]], beta[[3]], col = "green")
lines(alphaI[[4]], beta[[4]], col = "orange")
lines(alphaI[[5]], beta[[5]], col = "purple")

for (l in L_idx)
{
  roc.me1sum.fa[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1sum.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1sum.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1sum.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.fa[[l]] < alpha)/M
    alpha_1.me1sum.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.fa[[l]] < alpha)/M
    roc.me1sum.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.fa[[l]] < alpha)/M
    beta.me1sum.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.fa[[l]] < alpha)/M
    roc.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1sum.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

clrs <- c('black', 'red', 'green', 'orange', 'purple')
lwds <- c(2, 1, 1, 1, 1)

#jpeg(filename = "C://Temp//text_img//roc_me1_fa_full_sgl.jpg", height = 300, quality = 100)
plot(roc.me1sum.fa[[1]]$fpr, roc.me1sum.fa[[1]]$tpr, type = 'l', lwd = lwds[1], xlim=c(0, 1), ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
for (l in L_idx[-1])
  lines(roc.me1sum.fa[[l]]$fpr, roc.me1sum.fa[[l]]$tpr, col = clrs[l], type = 'l', lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
#dev.off()

plot(alpha_1.me1sum.fa[[1]]$alpha, alpha_1.me1sum.fa[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1sum.fa[[l]]$alpha, alpha_1.me1sum.fa[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)


roc.me1sum.ev <- list()
alpha_1.me1sum.ev <- list()
beta.me1sum.ev <- list()
for (l in L_idx)
{
  roc.me1sum.ev[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1sum.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1sum.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1sum.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.ev[[l]] < alpha)/M
    alpha_1.me1sum.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1sum.ev[[l]] < alpha)/M
    roc.me1sum.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.ev[[l]] < alpha)/M
    beta.me1sum.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1sum.ev[[l]] < alpha)/M
    roc.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1sum.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

clrs <- c('black', 'red', 'green', 'orange', 'purple')
lwds <- c(2, 1, 1, 1, 1)

#jpeg(filename = "C://Temp//text_img//roc_me1_fa_full_sgl.jpg", height = 300, quality = 100)
plot(roc.me1sum.ev[[1]]$fpr, roc.me1sum.ev[[1]]$tpr, type = 'l', lwd = lwds[1], xlim=c(0, 1), ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
for (l in L_idx[-1])
  lines(roc.me1sum.ev[[l]]$fpr, roc.me1sum.ev[[l]]$tpr, col = clrs[l], type = 'l', lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
#dev.off()

plot(alpha_1.me1sum.ev[[1]]$alpha, alpha_1.me1sum.ev[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1sum.ev[[l]]$alpha, alpha_1.me1sum.ev[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)



roc.me1block.ev <- list()
alpha_1.me1block.ev <- list()
beta.me1block.ev <- list()
for (l in L_idx)
{
  roc.me1block.ev[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1block.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1block.ev[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1block.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.ev[[l]] < alpha)/M
    alpha_1.me1block.ev[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.ev[[l]] < alpha)/M
    roc.me1block.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.ev[[l]] < alpha)/M
    beta.me1block.ev[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.ev[[l]] < alpha)/M
    roc.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1block.ev[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

clrs <- c('black', 'red', 'green', 'orange', 'purple')
lwds <- c(2, 1, 1, 1, 1)

#jpeg(filename = "C://Temp//text_img//roc_me1_fa_full_sgl.jpg", height = 300, quality = 100)
plot(roc.me1block.ev[[1]]$fpr, roc.me1block.ev[[1]]$tpr, type = 'l', lwd = lwds[1], xlim=c(0, 1), ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
for (l in L_idx[-1])
  lines(roc.me1block.ev[[l]]$fpr, roc.me1block.ev[[l]]$tpr, col = clrs[l], type = 'l', lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
#dev.off()

plot(alpha_1.me1block.ev[[1]]$alpha, alpha_1.me1block.ev[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1block.ev[[l]]$alpha, alpha_1.me1block.ev[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)



roc.me1block.fa <- list()
alpha_1.me1block.fa <- list()
beta.me1block.fa <- list()
for (l in L_idx)
{
  roc.me1block.fa[[as.character(Ls[l])]] <- data.frame(fpr=numeric(length(alphas)), tpr=numeric(length(alphas)), alpha=numeric(length(alphas)))
  alpha_1.me1block.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), fpr=numeric(length(alphas)))
  beta.me1block.fa[[as.character(Ls[l])]] <- data.frame(alpha=numeric(length(alphas)), tpr=numeric(length(alphas)))
  
  for (i in alphas_idx) 
  {
    alpha <- alphas[i]
    
    roc.me1block.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.fa[[l]] < alpha)/M
    alpha_1.me1block.fa[[as.character(Ls[l])]]$fpr[i] <- sum(p.values_noise.me1block.fa[[l]] < alpha)/M
    roc.me1block.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.fa[[l]] < alpha)/M
    beta.me1block.fa[[as.character(Ls[l])]]$tpr[i] <- sum(p.values_signal.me1block.fa[[l]] < alpha)/M
    roc.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    alpha_1.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
    beta.me1block.fa[[as.character(Ls[l])]]$alpha[i] <- alpha
  }
}

clrs <- c('black', 'red', 'green', 'orange', 'purple')
lwds <- c(2, 1, 1, 1, 1)

#jpeg(filename = "C://Temp//text_img//roc_me1_fa_full_sgl.jpg", height = 300, quality = 100)
plot(roc.me1block.fa[[1]]$fpr, roc.me1block.fa[[1]]$tpr, type = 'l', lwd = lwds[1], xlim=c(0, 1), ylim=c(0, 1), main = 'ROC curve', xlab = 'type I error', ylab = 'power')
for (l in L_idx[-1])
  lines(roc.me1block.fa[[l]]$fpr, roc.me1block.fa[[l]]$tpr, col = clrs[l], type = 'l', lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)
#dev.off()

plot(alpha_1.me1block.fa[[1]]$alpha, alpha_1.me1block.fa[[1]]$fpr, type = 'l', lwd = lwds[1], main = 'Type I error', xlab = 'significance level', ylab = 'type I error')
for (l in L_idx[-1])
  lines(alpha_1.me1block.fa[[l]]$alpha, alpha_1.me1block.fa[[l]]$fpr, col = clrs[l], lwd = lwds[l])
lines(c(0, 1), c(0, 1), type='l', lty = 3, col = 'blue')
legend(x = "bottomright", as.character(Ls), col = clrs, lty = 1, lwd = lwds)



























plot(punif)
lines(ecdf(p.values_noise.me1sum.ev[[1]])(alphas), ecdf(p.values_signal.me1sum.ev[[1]])(alphas))
lines(ecdf(p.values_noise.me1sum.ev[[2]])(alphas), ecdf(p.values_signal.me1sum.ev[[2]])(alphas), col = "blue")
lines(ecdf(p.values_noise.me1sum.ev[[3]])(alphas), ecdf(p.values_signal.me1sum.ev[[3]])(alphas), col = "violet")
lines(ecdf(p.values_noise.me1sum.ev[[4]])(alphas), ecdf(p.values_signal.me1sum.ev[[4]])(alphas), col = "green")
lines(ecdf(p.values_noise.me1sum.ev[[5]])(alphas), ecdf(p.values_signal.me1sum.ev[[5]])(alphas), col = "brown")

plot(punif)
lines(ecdf(p.values_noise.me1sum.fa[[1]])(alphas), ecdf(p.values_signal.me1sum.fa[[1]])(alphas))
lines(ecdf(p.values_noise.me1sum.fa[[2]])(alphas), ecdf(p.values_signal.me1sum.fa[[2]])(alphas), col = "blue")
lines(ecdf(p.values_noise.me1sum.fa[[3]])(alphas), ecdf(p.values_signal.me1sum.fa[[3]])(alphas), col = "violet")
lines(ecdf(p.values_noise.me1sum.fa[[4]])(alphas), ecdf(p.values_signal.me1sum.fa[[4]])(alphas), col = "green")
lines(ecdf(p.values_noise.me1sum.fa[[5]])(alphas), ecdf(p.values_signal.me1sum.fa[[5]])(alphas), col = "yellow")



lines(ecdf(p.values_noise.me1block.fa[[3]])(alphas), ecdf(p.values_signal.me1block.fa[[3]])(alphas), col = "green")


# cores <- detectCores()
# cl <- makeCluster(cores[1] - 1) #not to overload your computer
# registerDoParallel(cl)
# p <- foreach (iter(1), .combine = 'c', .export = c('toeplitz.mssa', 'Lcov', 'create.traj.mat', 'create.stacked.traj.mat', 'diag.avg', 
#                                                  'toeplitz.reconstruct', 'Toeplitz', 'eigen', 'new.hbhmat'), .packages = 'iterators') %dorng% {
#     f <- signal
#     s <- toeplitz.mssa(f, L = 20, D = 2, method = "block")
#     s$sigma
#     #@toeplitz.reconstruct(s, groups = list(cos=1:10))[[1]]
# }
# stopCluster(cl)









set.seed(Sys.time())
D <- 2
N <- 71
M <- 10000
signal <- ts(cbind(30 * cos(2*pi*(1:N)/12), 
                20 * cos(2*pi*(1:N)/12)))


Ls <- c(12, 24, 36, 48, 60)
MSEs <- cbind()
for (L in Ls) {
  MSE <- rbind()
  for (i in 1:100) {
    noise <- cbind(rnorm(N, sd = 5), rnorm(N, sd = 5))
    ts <- signal + noise
    s1 <- toeplitz.mssa(ts, L = L, D = 2, method = "sum")
    s2 <- toeplitz.mssa(ts, L = L, D = 2, method = "block")
    r1 <- toeplitz.reconstruct(s1, groups = list(cos = 1:2))[[1]]
    r2 <- toeplitz.reconstruct(s2, groups = list(cos = 1:2))[[1]]
    mse1 <- mean((signal - r1)^2)
    mse2 <- mean((signal - r2)^2)
    MSE <- rbind(MSE, cbind(mse1, mse2))
  }
  colMeans(MSE)
  MSEs <- cbind(MSEs, colMeans(MSE))
}
















