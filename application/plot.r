library(haven)
library(ggplot2)
library(ggpubr)
library(scales)

hist1 <- read_sas("history1_yuan_03aug2022.sas7bdat")
study0 <- read_sas("study0_yuan_03aug2022.sas7bdat")
sum.dts <- summary(c(hist1$dt1, hist1$dt2))

##########################
########## Plot ##########
##########################

# effect estimates and confidence interval
res <- read.csv("analysis_add_month_3m.csv", row.names = 1)
res <- res[-grep("region", rownames(res)), ]
res.TI <- res[1:17, ]
res.TV <- res[18:26, ]

res <- res.TI
res$covariate <- c("Age/5 (years)", "Sex (female)", "Diabetes (yes)", 
                   "Diagnosis (hepatitis C)", "(cholestatic cirrhosis)", "(acute hepatic nephropathy)", 
                   "(metastatic disease)", "(malignant neoplasm)", 
                   "Race (Black)", "(Hispanic)", "(Asian)", "(other)", 
                   "Malignancy", "Height/10 (cm)", "Weight/5 (kg)", 
                   "Functional status", "Working for income")
res$covariate <- factor(res$covariate, levels = unique(res$covariate))
res$signif <- ifelse(res$p < 0.001, "***", ifelse(res$p < 0.01, "**", ifelse(res$p < 0.05, "*", "")))
png(file = "effect_size_3m_TI.png", width = 12, height = 15, units = "cm", res = 1200)
ggplot(data = res, aes(x = betahat, y = covariate)) + 
  geom_point() + 
  geom_errorbar(aes(xmin = betahat - 1.96 * se, xmax = betahat + 1.96 * se), width = 0) + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.7) + 
  geom_text(aes(label = signif, x = betahat), size = 4, vjust = -0.2) + 
  scale_y_discrete(limits = rev(levels(res$covariate))) + 
  labs(x = expression(widehat(bold(beta))), y = "") + # tag = "Significance\n* p < 0.05\n** p < 0.01\n*** p < 0.001"
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0, "lines"), 
        # plot.tag.position = c(1.025, 0.12), 
        # plot.tag = element_text(family = "serif", lineheight = 1.25, hjust = 0, size = 10), 
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5, linetype = "solid"), 
        panel.grid.major = element_line(linetype = "dotted", color = "grey"), 
        panel.grid.minor = element_line(linetype = "dotted", color = "grey"))
dev.off()
