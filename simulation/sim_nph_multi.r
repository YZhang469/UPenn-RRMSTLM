source("function.r")

set.seed(241)

res12 <- sim(ph = FALSE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
             a = 1, b = 50, lambda0C = 0.0014, lambda0T = 0.0028, 
             method = "multi")

write.csv(res12, "sim_nph_multi_lowC.csv")

set.seed(242)

res16 <- sim(ph = FALSE, n.sim = 1000, n = 2500, K = 5, interval = 10, L = 50, 
             a = 1, b = 50, lambda0C = 0.004, lambda0T = 0.0075, 
             method = "multi")

write.csv(res16, "sim_nph_multi_highC.csv")
