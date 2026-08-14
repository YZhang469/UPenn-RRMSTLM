source("function_concordance.R")

sum.entry <- summary(study0$dt_wl)
sum.dts <- summary(c(hist1$dt1, hist1$dt2))

# monthly
day <- c("1/31/", "2/28/", "3/31/", "4/30/", "5/31/", "6/30/", "7/31/", "8/31/", "9/30/", "10/31/", "11/30/", "12/31/")
year <- 2010:(2010+floor((sum.dts["Max."]-sum.dts["Min."]+1)/365)-1-1)
CSk <- as.Date(do.call(paste0, expand.grid(day, year)), format = "%m/%d/%Y") - 
  as.Date("01/01/1960", format = "%m/%d/%Y")
CV.dat(dat = study0, datA = hist1, IDname = "WL_ID_CODE", TrtInelname = "inactive", Ename = "dt_wl", 
       ZCnames = ZCnames, ZTnames.TI = ZTnames.TI, ZTnames.TV = ZTnames.TV, 
       Znames.TI = Znames.TI, Znames.TV = Znames.TV, 
       CSk = CSk, L = 90, M = 2, stratified = "add")
