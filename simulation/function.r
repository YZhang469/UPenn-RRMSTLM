library(survival)
library(stabledist)
library(truncnorm)
library(dplyr)
library(data.table)
library(geepack)

generateData <- function(ph, 
                         n, K, 
                         a, b, # study entry
                         lambda0C, thetaC, # censoring time
                         lambda0A, c, # treatment eligibility
                         lambda0T # treatment time
                         ){
  
  # time-fixed covariate
  Za <- rbinom(n, 1, 0.5)
  
  if (ph){
    mu = 24; sigma = 1; rho = 0.8; # Vik
    gamma1 = -0.6; gamma2 = -0.3; # death time
    
    # time-varying covariate
    eta <- rnorm(n, mean = mu, sd = sigma)
    V <- rstable(n*K, alpha = rho, beta = 1, 
                 gamma = (a*cos(pi*rho/2))^(1/rho), delta = 0, pm = 1)
    Vmat <- matrix(V, nrow = n, ncol = K)
    Zb0 <- eta + apply(log(Vmat)/gamma2, 1, sum)
    Zbmat <- replicate(K, Zb0) - log(Vmat)/gamma2
    colnames(Zbmat) <- paste0("Zb", 1:K)
    # survival time
    V0 <- rstable(n, alpha = rho, beta = 1, 
                  gamma = (a*cos(pi*rho/2))^(1/rho), delta = 0, pm = 1)
    D <- a * (-log(runif(n, 0, 1)) * exp(-gamma1*Za-gamma2*Zb0) * V0^{-1/rho}) ^ (rho^2)
  }
  
  else if (!ph){
    mu = 36; sigma = 1; # Vik
    gamma1 = -0.6; gamma2 = -0.3; # death time
    
    # time-varying covariate
    eta <- rnorm(n, mean = mu, sd = sigma)
    V <- rtruncnorm(n*K, a = 0, b = Inf, mean = 0, sd = 1)
    Vmat <- matrix(V, nrow = n, ncol = K)
    Zb0 <- eta + apply(log(Vmat)/gamma2, 1, sum)
    Zbmat <- replicate(K, Zb0) - log(Vmat)/gamma2
    colnames(Zbmat) <- paste0("Zb", 1:K)
    # survival time
    V0 <- rtruncnorm(n, a = 0, b = Inf, mean = 0, sd = 2)
    D <- a * (-log(runif(n, 0, 1)) * exp(-gamma1*Za-gamma2*Zb0) * V0^{-1}) ^ 0.3
  }
  
  # study entry
  E = runif(n, min = 0, max = b)
  
  # censoring time
  C <- -log(runif(n, 0, 1))/lambda0C/exp(thetaC*Za)
  
  # treatment eligibility
  A <- -log(runif(n, 0, 1))/lambda0A/exp(c*V0)
  # treatment time
  thetaT1 = -0.1; thetaT2 = 0.1; 
  U <- runif(n, 0, 1)
  Trt <- rep(NA, n)
  for (i in 1:n){
    Zia <- Za[i]
    Ai <- A[i]
    f <- function(x){
      ft <- ifelse(x < Ai, 1, 0) * x + ifelse(x > Ai, 1, 0) * (Ai + exp(thetaT2)*(x-Ai))
      exp(-lambda0T * exp(thetaT1*Zia) * ft) - U[i]
    }
    Trt[i] <- try({uniroot(f, interval = c(1e-08, 1e+08), tol = 1e-6)$root}, silent = TRUE)
  }
  Trt <- suppressWarnings(as.numeric(Trt))
  Trt[is.na(Trt)|Trt > A] <- 1e+08
  
  dat <- cbind.data.frame("ID" = 1:n, "entry" = E, "Za" = Za, Zbmat, "X" = pmin(D, C, Trt), 
                          "deltaD" = ifelse(D < C & D < Trt, 1, 0), "deltaT" = ifelse(Trt < D & Trt < C, 1, 0), 
                          "A" = A)
  
  datA <- data.frame("ID" = rep(1:n, each = 2), "status" = rep(c("active", "inactive"), n), 
                     "Tstart" = rep(0, 2*n), "Tstop" = rep(dat$X, each = 2))
  datA$Tstop[c(T,F)] <- A
  datA$Tstart[c(F,T)] <- A
  datA <- datA %>% filter(Tstart < Tstop)
  
  return(list("datA" = datA, "dat" = dat))
}

trueBeta <- function(ph, L, K, interval, 
                     a, b, lambda0A, c, 
                     method = c("add", "multi")){
  n = 10000000
  
  # time-fixed covariate
  Za <- rbinom(n, 1, 0.5)
  
  if (ph){
    mu = 24; sigma = 1; rho = 0.8; # Vik
    gamma1 = -0.6; gamma2 = -0.3; # death time
    
    # time-varying covariate
    eta <- rnorm(n, mean = mu, sd = sigma)
    V <- rstable(n*K, alpha = rho, beta = 1, 
                 gamma = (a*cos(pi*rho/2))^(1/rho), delta = 0, pm = 1)
    Vmat <- matrix(V, nrow = n, ncol = K)
    Zb0 <- eta + apply(log(Vmat)/gamma2, 1, sum)
    Zbmat <- replicate(K, Zb0) - log(Vmat)/gamma2
    colnames(Zbmat) <- paste0("Zb", 1:K)
    
    # survival time
    V0 <- rstable(n, alpha = rho, beta = 1, 
                  gamma = (a*cos(pi*rho/2))^(1/rho), delta = 0, pm = 1)
    D <- a * (-log(runif(n, 0, 1)) * exp(-gamma1*Za-gamma2*Zb0) * V0^{-1/rho}) ^ (rho^2)
  }
  
  else if (!ph){
    mu = 36; sigma = 1; # Vik
    gamma1 = -0.6; gamma2 = -0.3; # death time
    
    # time-varying covariate
    eta <- rnorm(n, mean = mu, sd = sigma)
    V <- rtruncnorm(n*K, a = 0, b = Inf, mean = 0, sd = 1)
    Vmat <- matrix(V, nrow = n, ncol = K)
    Zb0 <- eta + apply(log(Vmat)/gamma2, 1, sum)
    Zbmat <- replicate(K, Zb0) - log(Vmat)/gamma2
    colnames(Zbmat) <- paste0("Zb", 1:K)
    # survival time
    V0 <- rtruncnorm(n, a = 0, b = Inf, mean = 0, sd = 2)
    D <- a * (-log(runif(n, 0, 1)) * exp(-gamma1*Za-gamma2*Zb0) * V0^{-1}) ^ 0.3
  }
  
  # assume patients enter the study between time 0 and b
  E = runif(n, min = 0, max = b)
  
  # assume patients are eligible for treatment up to A
  A = -log(runif(n, 0, 1))/lambda0A/exp(c*V0)
  
  dat <- data.frame("ID" = rep(1:n, K), "CS" = rep(1:K, each = n), "Za" = rep(Za, K), "Zb" = c(Zbmat), 
                    "Sik" = unlist(lapply(1:K, function(k){interval*k - E})), 
                    "Dik" = unlist(lapply(1:K, function(k){D - (interval*k - E)})), 
                    "Aik" = unlist(lapply(1:K, function(k){ifelse(interval*k - E < A, 1, 0)})))
  dat <- dat[dat$Sik >= 0 & dat$Dik >= 0 & dat$Aik == 1, ]
  dat$DLik <- pmin(dat$Dik, rep(L, nrow(dat)))

  if (method == "add"){
    Zabar <- tapply(dat$Za, dat$CS, mean)
    Zbbar <- tapply(dat$Zb, dat$CS, mean)
    dat.temp <- dat
    dat.temp$Za = dat$Za - Zabar[dat.temp$CS]
    dat.temp$Zb = dat$Zb - Zbbar[dat.temp$CS]
    beta <- coef(lm(DLik ~ Za + Zb - 1, dat = dat.temp))
  }
  else if (method == "multi"){
    dat.temp <- dat
    dat.temp$deltaX <- 1
    dat.temp$X <- 1
    dat.temp$W <- dat.temp$DLik
    mod.multi <- coxph(Surv(X, deltaX) ~ Za + Zb + offset(-log(DLik)) + strata(CS), data = dat.temp, 
                       weights = W, ties = "breslow", id = ID)
    beta <- coefficients(mod.multi)
  }
  
  return(beta)
}

dp <- function(dat, datA, L, K, interval,  # K = the number of cross-sections; interval = the time between consecutive cross-sections
               method = c("add", "multi")){
  
  # estimate the censoring model
  datC <- dat
  modC <- coxph(Surv(X, 1-(deltaD+deltaT)) ~ Za, data = datC, ties = "breslow")
  
  # estimate the treatment model
  datT <- merge(dat, datA, by = "ID")
  datT$weights <- ifelse(datT$status == "active", 1, 1e-08)
  dtimes <- sort(unique(with(datT, X[deltaT == 1])))
  datT.Extended <- survSplit(Surv(X, deltaT == 1) ~ ., datT, cut = dtimes)
  datT.Extended[, "Zt"] <- ifelse(datT.Extended$X > datT.Extended$A, 1, 0)
  modT <- coxph(Surv(tstart, X, event) ~ Za + Zt, data = datT.Extended, weights = weights, 
                ties = "breslow", timefix = FALSE)
  
  # create the stacked dataset
  dat.stacked <- data.frame()
  for (k in 1:K){
    # select subjects who are at-risk and treatment-eligible at Sik
    dat.temp <- dat
    dat.temp$Sik <- interval*k - dat.temp$entry
    dat.temp <- dat.temp %>%
      filter(Sik > 0 & Sik < X & Sik < A)
    if (nrow(dat.temp) > 0){
      dat.temp <- dat.temp %>%
        mutate(Xik = X - Sik) %>%
        rowwise() %>%
        mutate(Yik = min(Xik, L), deltaYik = ifelse(Xik <= L, deltaD, 1)) %>%
        select(-paste("Zb", (1:K)[-k], sep = ""))
      colnames(dat.temp)[grepl(paste("Zb", k, sep = ""), colnames(dat.temp))] <- "Zb"
      dat.stacked <- rbind.data.frame(dat.stacked, cbind.data.frame("CS" = k, dat.temp))
    }
    else{
      dat.stacked <- dat.stacked
    }
  }
  
  # calculate the inverse probability weights
  dat.temp1 = dat.temp2 = dat.temp3 = dat.stacked
  dat.temp1$X <- dat.temp1$Sik
  dat.temp2$X <- dat.temp2$Sik + dat.temp2$Yik
  dat.stacked$WC <- predict(modC, newdata = dat.temp1, type = "survival")/predict(modC, newdata = dat.temp2, type = "survival")
  
  dat.temp1$Zt <- ifelse(dat.temp1$A < dat.temp1$Sik, 1, 0)
  dat.temp1$tstart <- 0
  # dat.temp1$X <- dat.temp1$Sik
  dat.temp1$event <- dat.temp1$deltaT
  dat.temp2$Zt <- ifelse(dat.temp2$A < dat.temp2$Sik + dat.temp2$Yik, 1, 0)
  dat.temp2$tstart <- 0
  # dat.temp2$X <- dat.temp2$Sik + dat.temp2$Yik
  dat.temp2$event <- dat.temp2$deltaT
  dat.temp3$Zt <- ifelse(dat.temp3$A < dat.temp3$Yik, 1, 0)
  dat.temp3$tstart <- 0
  dat.temp3$X <- dat.temp3$Yik
  dat.temp3$event <- dat.temp3$deltaT

  dat.stacked$WT <- predict(modT, newdata = dat.temp1, type = "survival")/predict(modT, newdata = dat.temp2, type = "survival")
  dat.stacked$W <- dat.stacked$WC * dat.stacked$WT
  
  # point estimate and standard error of betahat
  dat.stacked$Wt <- dat.stacked$W * dat.stacked$deltaYik
  if (method == "add"){
    Zbar <- cbind.data.frame("Za" = sapply(split(dat.stacked, dat.stacked$CS), function(x){weighted.mean(x$Za, x$Wt)}), 
                             "Zb" = sapply(split(dat.stacked, dat.stacked$CS), function(x){weighted.mean(x$Zb, x$Wt)}))
    dat.temp <- dat.stacked
    dat.temp$Za = dat.stacked$Za - Zbar[dat.temp$CS, 1]
    dat.temp$Zb = dat.stacked$Zb - Zbar[dat.temp$CS, 2]
    betahat <- coef(lm(Yik ~ Za + Zb - 1, dat = dat.temp, weights = Wt))
    Zres <- dat.temp[, c("Za", "Zb")]
    A <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    for (i in 1:nrow(dat.stacked)){
      A = A + t(as.matrix(Zres[i, ])) %*% as.matrix(Zres[i, ]) * dat.stacked$Wt[i] / nrow(dat)
    }
    dat.temp$Yres <- dat.stacked$Yik - as.matrix(dat.stacked[, c("Za", "Zb")]) %*% betahat
    mu0k <- sapply(split(dat.temp, dat.temp$CS), function(x){weighted.mean(x$Yres, x$Wt)})
    B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    e <- diag(as.vector(dat.stacked$Wt * (dat.stacked$Yik - mu0k[dat.stacked$CS] - as.matrix(dat.stacked[, c("Za", "Zb")]) %*% betahat))) %*% as.matrix(Zres)
    E <- cbind.data.frame("ID" = dat.stacked$ID, e)
    Ei <- aggregate(. ~ ID, E, sum)
    for (i in 1:nrow(Ei)){
      ei <- t(Ei[i, -1])
      B <- B + ei %*% t(ei) / nrow(dat)
    }
    var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat)
  }
  else if (method == "multi"){
    dat.stacked$W.new <- dat.stacked$Wt * dat.stacked$Yik
    dat.temp <- dat.stacked
    dat.temp <- dat.temp[!(is.na(dat.temp$W.new)|dat.temp$W.new == 0), ]
    dat.temp$deltaX.new <- 1
    dat.temp$X.new <- 1
    mod.multi <- coxph(Surv(X.new, deltaX.new) ~ Za + Zb + offset(-log(Yik)) + strata(CS), data = dat.temp, 
                       weights = W.new, ties = "breslow", id = ID)
    betahat <- coefficients(mod.multi)
    Z <- dat.stacked[, c("Za", "Zb")]
    Zbar <- dat.stacked$Wt * exp(as.matrix(Z) %*% betahat)
    S0 <- tapply(Zbar, as.factor(dat.stacked$CS), sum)
    S1 <- aggregate(Z * Zbar, by = list(dat.stacked$CS), sum)[, -1]
    S2 <- aggregate(data.frame(t(apply(as.matrix(Z), 1, tcrossprod))) * Zbar, by = list(dat.stacked$CS), sum)[, -1]
    mu0k <- tapply(dat.stacked$W.new, as.factor(dat.stacked$CS), sum)/S0
    Sbar <- S1/S0
    Sk <- S2/S0 - t(apply(as.matrix(Sbar), 1, tcrossprod))
    A <- matrix(colSums(aggregate(Sk[dat.stacked$CS, ] * dat.stacked$W.new, by = list(dat.stacked$ID), sum)[, -1]) / nrow(dat), 
                nrow = length(betahat), ncol = length(betahat))
    B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    e <- dat.stacked$Wt * (Z - Sbar[dat.stacked$CS, ]) * (dat.stacked$Yik - mu0k[dat.stacked$CS] * as.vector(exp(as.matrix(Z) %*% betahat)))
    E <- cbind.data.frame("ID" = dat.stacked$ID, e)
    Ei <- aggregate(. ~ ID, E, sum)
    for (i in 1:nrow(Ei)){
      ei <- t(Ei[i, -1])
      B <- B + ei %*% t(ei) / nrow(dat)
    }
    var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat)
  }
  
  return(list("betahat" = betahat, "se" = sqrt(var)))
}

sim <- function(ph, n.sim, n, K, interval, L, 
                a, b, 
                lambda0C, 
                lambda0T, 
                method){ # vary L and percentage of censoring
  
  beta <- trueBeta(ph = ph, L = L, K = K, interval = interval, 
                   a = a, b = b, lambda0A = 0.002, c = 0.001, 
                   method = method)
  est <- data.frame(matrix(ncol = length(beta), nrow = n.sim))
  error <- data.frame(matrix(ncol = length(beta), nrow = n.sim))
  cover <- data.frame(matrix(ncol = length(beta), nrow = n.sim))
  for (iter in 1:n.sim){
    print(iter)
    dat.list <- generateData(ph = ph, n = n, K = K, 
                             a = a, b = b, lambda0C = lambda0C, thetaC = 0.001, 
                             lambda0A = 0.002, c = 0.001, lambda0T = lambda0T)
    dat <- dat.list$dat
    datA <- dat.list$datA
    res <- dp(dat, datA, L = L, K = K, interval = interval, 
              method = method)
    betahat <- res$betahat
    se <- res$se
    ci <- ifelse(beta > betahat-1.96*se & beta < betahat+1.96*se, 1, 0)
    est[iter, ] <- betahat
    error[iter, ] <- se
    cover[iter, ] <- ci
  }
  ## metrics
  # bias
  bias <- apply(est, 2, mean, na.rm = T) - beta
  # empirical standard deviation (ESD)
  esd <- sqrt(apply(apply(est, 1, function(x){(x - beta)^2}), 1, mean, na.rm = T))
  # average asymptotic standard error (ASE)
  ase <- apply(error, 2, mean, na.rm = T)
  # empirical coverage probabilities (CP): use ASE to create a confidence interval, and test if the true value lies within the interval
  cp <- apply(cover, 2, mean, na.rm = T)
  res <- cbind("true" = beta, "bias" = bias, "esd" = esd, "ase" = ase, "cp" = cp)
  
  return(res)
}
