source("function.r")

set.seed(231)

res11 <- sim(ph = FALSE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
             a = 1, b = 50, lambda0C = 0.0014, lambda0T = 0.0028, 
             method = "add")

write.csv(res11, "sim_nph_add_2500.csv")

set.seed(232)

res15 <- sim(ph = FALSE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
             a = 1, b = 50, lambda0C = 0.004, lambda0T = 0.0075, 
             method = "add")

write.csv(res15, "sim_nph_add_highC_2500.csv")
