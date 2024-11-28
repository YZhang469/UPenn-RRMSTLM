source("function.r")

set.seed(131)

res3 <- sim(ph = TRUE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
            a = 1, b = 50, lambda0C = 0.0012, lambda0T = 0.0024, 
            method = "add")

write.csv(res3, "sim_ph_add_2500.csv")

set.seed(132)

res7 <- sim(ph = TRUE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
            a = 1, b = 50, lambda0C = 0.004, lambda0T = 0.0077, 
            method = "add")

write.csv(res7, "sim_ph_add_highC_2500.csv")
