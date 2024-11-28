source("function.r")

set.seed(141)

res4 <- sim(ph = TRUE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
            a = 1, b = 50, lambda0C = 0.0012, lambda0T = 0.0024, 
            method = "multi")

write.csv(res4, "sim_ph_multi_lowC.csv")

set.seed(142)

res8 <- sim(ph = TRUE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
            a = 1, b = 50, lambda0C = 0.004, lambda0T = 0.0077, 
            method = "multi")

write.csv(res8, "sim_ph_multi_highC.csv")
