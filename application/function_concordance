library(haven)
library(survival)
library(dplyr)
library(data.table)
library(geepack)

IOC <- function(predicted, time, event, dat, modC, method){ # method = "Harrell" or "Uno" (with IPTW)
  if (method == "Uno"){
    weights <- 1/predict(modC, newdata = dat, type = "survival")
  }
  num = 0
  denom = 0
  for (i in 1:nrow(dat)){
    if (method == "Harrell"){
      num <- num + sum(event[i] * ifelse(time[i] < time, 1, 0) * ifelse(predicted[i] < predicted, 1, 0))
      denom <- denom + sum(event[i] * ifelse(time[i] < time, 1, 0))
    }
    else if (method == "Uno"){
      dat.temp <- dat
      dat.temp$X <- time[i]
      weights.temp <- 1/predict(modC, newdata = dat.temp, type = "survival")
      num <- num + sum(event[i] * weights[i] * weights.temp * ifelse(time[i] < time, 1, 0) * ifelse(predicted[i] < predicted, 1, 0))
      denom <- denom + sum(event[i] * weights[i] * weights.temp * ifelse(time[i] < time, 1, 0))
    }
  }
  return(num/denom)
}

## load data
hist1 <- read_sas("history1_yuan_03aug2022.sas7bdat")
hist1$bili <- ifelse(hist1$bili == 0, log(0.01), log(hist1$bili))
hist1$creat <- log(pmin(hist1$creat+1,5))
hist1$inr <- log(hist1$inr)
study0 <- read_sas("study0_yuan_03aug2022.sas7bdat")
study0$bili_wl <- ifelse(study0$bili_wl == 0, log(0.01), log(study0$bili_wl))
study0$creat_wl <- log(pmin(study0$creat_wl+1, 5))
study0$inr_wl <- log(study0$inr_wl)

ZCnames <- c("inactive_wl", "age_wl5", "female", "diabetes", 
             "blood_a", "blood_ab", "blood_b", # reference: "blood_o"
             "diag_HCV", "diag_chol_cirr", "diag_ahn", "diag_met_dis", "diag_mal_neo", # reference: "diag_nonchol_cirr"
             "race_Black", "race_Hisp", "race_Asian", "race_oth", # reference: "race_White"
             "region_1", "region_2", "region_3", "region_4", "region_5", 
             "region_6", "region_7", "region_8", "region_9", "region_10", # reference: "region_11"
             "malig_wl", "dialysis_wl", "albumin_wl", "ascites_wl", "enceph_wl", 
             "height10", "weight5", "funcstat_wl", "working_wl", "dt_wl", 
             "creat_wl", "bili_wl", "inr_wl", "sodium_wl") # baseline covariates (study0)
ZTnames.TI <- c("age_wl5", "female", "diabetes", 
                "diag_HCV", "diag_chol_cirr", "diag_ahn", "diag_met_dis", "diag_mal_neo", 
                "race_Black", "race_Hisp", "race_Asian", "race_oth", 
                "region_1", "region_2", "region_3", "region_4", "region_5", 
                "region_6", "region_7", "region_8", "region_9", "region_10", 
                "malig_wl", "height10", "weight5", "funcstat_wl", "working_wl")
ZTnames.TV <- c("dialysis", "ascites", "enceph", "albumin", "meldall", "sodium") # time-varying covariates
Znames.TI <- c("age_wl5", "female", "diabetes", 
               "diag_HCV", "diag_chol_cirr", "diag_ahn", "diag_met_dis", "diag_mal_neo", 
               "race_Black", "race_Hisp", "race_Asian", "race_oth", 
               "region_1", "region_2", "region_3", "region_4", "region_5", 
               "region_6", "region_7", "region_8", "region_9", "region_10", 
               "malig_wl", "height10", "weight5", "funcstat_wl", "working_wl")
Znames.TV <- c("dialysis", "ascites", "enceph", "albumin", "bili", "creat", "inr", "sodium", "Sik") # time-varying covariates

CV.dat <- function(dat = study0, datA = hist1, IDname, TrtInelname, Ename, 
                   ZCnames, ZTnames.TI, ZTnames.TV, Znames.TI, Znames.TV, 
                   CSk, L, M = 5, stratified = "add"){
  
  ## create the stacked dataset without weights
  # impute missing covariate values as median
  names <- unique(c(ZCnames, ZTnames.TI, Znames.TI))
  for (j in 1:length(names)){
    if (sum(is.na(dat[, names[j]])) > 0){
      dat[is.na(dat[[names[j]]]), names[j]] <- median(dat[[names[j]]], na.rm = TRUE)
    }
  }
  
  dat.stacked <- data.frame()
  K <- length(CSk)
  for (k in 1:K){
    dat.temp <- merge(dat[, c(IDname, unique(c(ZCnames, ZTnames.TI, Znames.TI, Ename)), "X")], datA, by = IDname)
    dat.temp$Sik <- as.numeric(CSk[k] - dat.temp[, Ename])
    dat.temp <- dat.temp %>% 
      filter(Sik >= 0 & Sik <= X & Sik >= t1 & Sik < t2 & inactive == 0)
    if (nrow(dat.temp) > 0){
      dat.temp <- dat.temp %>%
        mutate(Xik = X - Sik) %>%
        rowwise() %>%
        mutate(Yik = min(Xik, L), deltaYik = ifelse(Xik <= L, dead, 1))
      dat.stacked <- rbind.data.frame(dat.stacked, cbind.data.frame("CS" = k, dat.temp))
    }
    else{
      dat.stacked <- dat.stacked
    }
  }
  
  ## split data into M folds based on subjects
  dat <- dat[sample(nrow(dat)), ]
  dat$index <- rep(1:M, length.out = nrow(dat))
  datA <- merge(dat[, c(IDname, "index")], datA, by = IDname)
  dat.stacked <- merge(dat[, c(IDname, "index")], dat.stacked, by = IDname)
  
  res <- rep(NA, M)
  
  for (m in 1:M){
    
    ## separate the training and test datasets
    train <- dat[dat$index != m, ]
    trainA <- datA[datA$index != m, ]
    train.stacked <- dat.stacked[dat.stacked$index != m, ]
    
    # estimate censoring model
    datC <- train
    modC <- coxph(as.formula(paste("Surv(X, 1-(dead+LT)) ~ ", paste(ZCnames, collapse = " + "), sep = "")), 
                  data = datC, ties = "breslow")
    
    # estimate treatment model
    datT <- merge(train[, c(IDname, ZTnames.TI, "X")], trainA, by = IDname)
    datT$weights <- ifelse(datT[, TrtInelname] == 0, 1, 1e-08)
    modT <- coxph(as.formula(paste("Surv(t1, t2, LT) ~ ", paste(c(ZTnames.TI, ZTnames.TV), collapse = " + "), sep = "")), 
                  data = datT, weights = weights, method = "breslow")
    
    dat.temp1 = dat.temp2 = train.stacked
    dat.temp1$X <- dat.temp1$Sik
    dat.temp2$X <- dat.temp2$Sik + dat.temp2$Yik
    train.stacked$WC <- predict(modC, newdata = dat.temp1, type = "survival")/predict(modC, newdata = dat.temp2, type = "survival")
    # some treatment weights are exceptionally large: there are strict rules for liver allocation, but centers sometimes "bend" them, resulting in the occurrence of highly improbably transplants
    dat.temp1$t1 <- 0
    dat.temp1$t2 <- ifelse(dat.temp1$Sik == 0, 1, dat.temp1$Sik)
    dat.temp2$t1 <- 0
    dat.temp2$t2 <- ifelse(dat.temp2$Sik + dat.temp2$Yik == 0, 1, dat.temp2$Sik + dat.temp2$Yik)
    train.stacked$WT <- predict(modT, newdata = dat.temp1, type = "survival")/predict(modT, newdata = dat.temp2, type = "survival")
    
    train.stacked$W <- train.stacked$WC * train.stacked$WT
    train.stacked$W[train.stacked$W > 100] <- 100
    
    train.stacked$Wt <- train.stacked$W * train.stacked$deltaYik
    if (stratified == "add"){
      if (length(c(Znames.TI, Znames.TV)) > 1){
        Zbar <- t(sapply(split(train.stacked, train.stacked$CS), function(x){colSums(x[, c(Znames.TI, Znames.TV)] * x$Wt)/sum(x$Wt)}))
      }
      else if (length(c(Znames.TI, Znames.TV)) == 1){
        Zbar <- matrix(sapply(split(train.stacked, train.stacked$CS), function(x){sum(x[, c(Znames.TI, Znames.TV)] * x$Wt)/sum(x$Wt)}))
      }
      dat.temp <- train.stacked
      dat.temp[, c(Znames.TI, Znames.TV)] = train.stacked[, c(Znames.TI, Znames.TV)] - Zbar[dat.temp$CS, ]
      betahat <- coef(lm(as.formula(paste("Yik ~ ", paste(c(Znames.TI, Znames.TV), collapse = " + "), " - 1", sep = "")), 
                         dat = dat.temp, weights = Wt))
    }
    else if (stratified == "multi"){
      train.stacked$W.new <- train.stacked$Wt * train.stacked$Yik
      dat.temp <- train.stacked
      dat.temp <- dat.temp[!(is.na(dat.temp$W.new)|dat.temp$W.new == 0), ]
      dat.temp$deltaX.new <- 1
      dat.temp$X.new <- 1
      mod.multi <- coxph(as.formula(paste("Surv(X.new, deltaX.new) ~ ", paste(c(Znames.TI, Znames.TV), collapse = " + "), " + offset(-log(Yik)) + strata(CS)", sep = "")), 
                         data = dat.temp, weights = W.new, ties = "breslow", id = WL_ID_CODE)
      betahat <- coefficients(mod.multi)
    }
    
    test.stacked <- dat.stacked[dat.stacked$index == m, ]
    Z.test <- test.stacked[, c(Znames.TI, Znames.TV)]
    
    if(stratified == "add"){
      predicted <- as.matrix(Z.test) %*% betahat
    }
    else if (stratified == "multi"){
      predicted <- exp(as.matrix(Z.test) %*% betahat)
    }
    res[m] <- IOC(predicted, test.stacked$Yik, test.stacked$deltaYik, test.stacked, modC, "Harrell")
  }
  
  return(mean(res))
  
}
