library(haven)
library(survival)
library(dplyr)

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
             "bili_wl", "creat_wl", "inr_wl", "sodium_wl") # baseline covariates (study0)
ZTnames.TI = Znames.TI = c("age_wl5", "female", "diabetes", 
                           "diag_HCV", "diag_chol_cirr", "diag_ahn", "diag_met_dis", "diag_mal_neo", 
                           "race_Black", "race_Hisp", "race_Asian", "race_oth", 
                           "region_1", "region_2", "region_3", "region_4", "region_5", 
                           "region_6", "region_7", "region_8", "region_9", "region_10", 
                           "malig_wl", "height10", "weight5", "funcstat_wl", "working_wl")
ZTnames.TV <- c("dialysis", "ascites", "enceph", "albumin", "meldall", "sodium") # time-varying covariates
Znames.TV <- c("dialysis", "ascites", "enceph", "albumin", "bili", "creat", "inr", "sodium", "Sik") # time-varying covariate

generateStackedData <- function(dat, datA, IDname, TrtInelname, Ename, 
                                ZCnames, ZTnames.TI, ZTnames.TV, Znames.TI, Znames.TV, 
                                CSk, L){
  
  # impute missing covariate values as median
  names <- unique(c(ZCnames, ZTnames.TI, Znames.TI))
  for (j in 1:length(names)){
    if (sum(is.na(dat[, names[j]])) > 0){
      dat[is.na(dat[[names[j]]]), names[j]] <- median(dat[[names[j]]], na.rm = TRUE)
    }
  }
  
  # estimate censoring model
  datC <- dat
  modC <- coxph(as.formula(paste("Surv(X, 1-(dead+LT)) ~ ", paste(ZCnames, collapse = " + "), sep = "")), data = datC, ties = "breslow")
  summary(modC)
  # estimate treatment model
  datT <- merge(dat[, c(IDname, ZTnames.TI, "X")], datA, by = IDname)
  datT$weights <- ifelse(datT[, TrtInelname] == 0, 1, 1e-08)
  modT <- coxph(as.formula(paste("Surv(t1, t2, LT) ~ ", paste(c(ZTnames.TI, ZTnames.TV), collapse = " + "), sep = "")), 
                data = datT, weights = weights, method = "breslow")
  # summary(modT)
  # write.csv(summary(modT)$coefficients, "modT.csv", row.names = TRUE)
  
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
  
  dat.temp1 = dat.temp2 = dat.stacked
  dat.temp1$X <- dat.temp1$Sik
  dat.temp2$X <- dat.temp2$Sik + dat.temp2$Yik
  dat.stacked$WC <- predict(modC, newdata = dat.temp1, type = "survival")/predict(modC, newdata = dat.temp2, type = "survival")
  dat.temp3 = dat.stacked
  dat.temp3$t1 <- dat.temp3$Sik
  dat.temp3$t2 <- dat.temp3$Sik + dat.temp3$Yik
  # some treatment weights are exceptionally large: there are strict rules for liver allocation, but centers sometimes "bend" them, resulting in the occurrence of highly improbably transplants
  dat.stacked$WT <- 1/predict(modT, newdata = dat.temp3, type = "survival")
  dat.stacked$W <- dat.stacked$WC * dat.stacked$WT
  dat.stacked$W[dat.stacked$W > 100] <- 100
  
  return(dat.stacked)
}

# point estimates and inference of model parameters
analyzeStackedData <- function(dat, dat.stacked, Znames.TI, Znames.TV, method){
  dat.stacked$Wt <- dat.stacked$W * dat.stacked$deltaYik
  if (method == "add"){
    Zbar <- t(sapply(split(dat.stacked, dat.stacked$CS), function(x){colSums(x[, c(Znames.TI, Znames.TV)] * x$Wt)/sum(x$Wt)}))
    dat.temp <- dat.stacked
    dat.temp[, c(Znames.TI, Znames.TV)] = dat.stacked[, c(Znames.TI, Znames.TV)] - Zbar[dat.temp$CS, ]
    betahat <- coef(lm(as.formula(paste("Yik ~ ", paste(c(Znames.TI, Znames.TV), collapse = " + "), " - 1", sep = "")), 
                       dat = dat.temp, weights = Wt))
    Zres <- dat.temp[, c(Znames.TI, Znames.TV)]
    A <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    for (i in 1:nrow(dat.stacked)){
      A = A + t(as.matrix(Zres[i, ])) %*% as.matrix(Zres[i, ]) * dat.stacked$Wt[i] / nrow(dat)
    }
    dat.temp$Yres <- dat.stacked$Yik - as.matrix(dat.stacked[, c(Znames.TI, Znames.TV)]) %*% betahat
    mu0k <- sapply(split(dat.temp, dat.temp$CS), function(x){weighted.mean(x$Yres, x$Wt)})
    B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    e <- Zres * as.vector(dat.stacked$Wt * (dat.stacked$Yik - mu0k[dat.stacked$CS] - as.matrix(dat.stacked[, c(Znames.TI, Znames.TV)]) %*% betahat))
    E <- cbind.data.frame("ID" = dat.stacked$WL_ID_CODE, e)
    Ei <- aggregate(. ~ ID, E, sum)
    for (i in 1:nrow(Ei)){
      ei <- t(Ei[i, -1])
      B <- B + ei %*% t(ei) / nrow(dat)
    }
    var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat)
  }
  else if (stratified == "multi"){
    dat.stacked$W.new <- dat.stacked$Wt * dat.stacked$Yik
    dat.temp <- dat.stacked
    dat.temp <- dat.temp[!(is.na(dat.temp$W.new)|dat.temp$W.new == 0), ]
    dat.temp$deltaX.new <- 1
    dat.temp$X.new <- 1
    mod.multi <- coxph(as.formula(paste("Surv(X.new, deltaX.new) ~ ", paste(c(Znames.TI, Znames.TV), collapse = " + "), " + offset(-log(Yik)) + strata(CS)", sep = "")), 
                       data = dat.temp, weights = W.new, ties = "breslow", id = WL_ID_CODE)
    betahat <- coefficients(mod.multi)
    Z <- dat.stacked[, c(Znames.TI, Znames.TV)]
    Zbar <- dat.stacked$Wt * exp(as.matrix(Z) %*% betahat)
    S0 <- tapply(Zbar, as.factor(dat.stacked$CS), sum)
    S1 <- aggregate(Z * Zbar, by = list(dat.stacked$CS), sum)[, -1]
    S2 <- aggregate(data.frame(t(apply(as.matrix(Z), 1, tcrossprod))) * Zbar, by = list(dat.stacked$CS), sum)[, -1]
    mu0k <- tapply(dat.stacked$W.new, as.factor(dat.stacked$CS), sum)/S0
    Sbar <- S1/S0
    Sk <- S2/S0 - t(apply(as.matrix(Sbar), 1, tcrossprod))
    A <- matrix(colSums(aggregate(Sk[dat.stacked$CS, ] * dat.stacked$W.new, by = list(dat.stacked$WL_ID_CODE), sum)[, -1]) / nrow(dat), 
                nrow = length(betahat), ncol = length(betahat))
    B <- matrix(0, ncol = length(betahat), nrow = length(betahat))
    e <- dat.stacked$Wt * (Z - Sbar[dat.stacked$CS, ]) * (dat.stacked$Yik - mu0k[dat.stacked$CS] * as.vector(exp(as.matrix(Z) %*% betahat)))
    E <- cbind.data.frame("ID" = dat.stacked$WL_ID_CODE, e)
    Ei <- aggregate(. ~ ID, E, sum)
    for (i in 1:nrow(Ei)){
      ei <- t(Ei[i, -1])
      B <- B + ei %*% t(ei) / nrow(dat)
    }
    var <- diag(solve(A) %*% B %*% t(solve(A))) / nrow(dat)
  }
  
  p <- pnorm(-abs(betahat/sqrt(var)), mean = 0, sd = 1, lower.tail = TRUE) * 2
  
  return(cbind.data.frame("betahat" = betahat, "se" = sqrt(var), "p" = p))
}

analyzeData <- function(dat = study0, datA = hist1, 
                        IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                        ZCnames, ZTnames.TI, ZTnames.TV, Znames.TI, Znames.TV, 
                        CSk, L, method){
  
  dat.stacked <- generateStackedData(dat, datA, IDname, TrtInelname, Ename, 
                                     ZCnames, ZTnames.TI, ZTnames.TV, Znames.TI, Znames.TV, 
                                     CSk, L)
  res <- analyzeStackedData(dat, dat.stacked, Znames.TI, Znames.TV, method)
  
  return(res)
}

sum.dts <- summary(c(hist1$dt1, hist1$dt2))

# monthly
day <- c("1/31/", "2/28/", "3/31/", "4/30/", "5/31/", "6/30/", "7/31/", "8/31/", "9/30/", "10/31/", "11/30/", "12/31/")
year <- 2010:(2010+floor((sum.dts["Max."]-sum.dts["Min."]+1)/365)-1-1)
CSk <- as.Date(do.call(paste0, expand.grid(day, year)), format = "%m/%d/%Y") - 
  as.Date("01/01/1960", format = "%m/%d/%Y")
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 30, method = "add")
write.csv(res.add, "analysis_add_month_1m.csv")
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 90, method = "add")
write.csv(res.add, "analysis_add_month_3m.csv")
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 180, method = "add")
write.csv(res.add, "analysis_add_month_6m.csv")
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 1*365, method = "add")
write.csv(res.add, "analysis_add_month_1y.csv")

# weekly
# CSk <- sum.dts["Min."] + (1:floor((sum.dts["Max."]-sum.dts["Min."]-365)/7)) * 7
CSk <- sum.dts["Min."] + (1:520) * 7
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 30, method = "add")
write.csv(res.add, "analysis_add_week_1m.csv")
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 90, method = "add")
write.csv(res.add, "analysis_add_week_3m.csv")
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 180, method = "add")
write.csv(res.add, "analysis_add_week_6m.csv")
res.add <- analyzeData(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 1*365, method = "add")
write.csv(res.add, "analysis_add_week_1y.csv")

# quarterly
day <- c("3/31/", "6/30/", "9/30/", "12/31/")
year <- 2010:(2010+floor((sum.dts["Max."]-sum.dts["Min."]+1)/365)-1-1)
CSk <- as.Date(do.call(paste0, expand.grid(day, year)), format = "%m/%d/%Y") - 
  as.Date("01/01/1960", format = "%m/%d/%Y")
res.add <- analyzeData(dat = study0, datA = hist1, 
                       IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
                       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, 
                       Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
                       CSk = CSk, L = 90, method = "add")
write.csv(res.add, "analysis_add_quarter.csv")
