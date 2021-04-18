data <- PhDArticles
system.time(best1 <- search_all(data, join_unobserved = TRUE))

all <- all_orders(data = data, join_unobserved = TRUE)
bics <- sapply(all, BIC)
system.time(best2 <- search_dynamic(data, join_unobserved = TRUE))

system.time(approx_best <- search_greedy(data, join_unobserved = TRUE))
best1
best2
approx_best
range(bics)

BIC(best2, approx_best)

data <-readRDS("FallEld.rds")
data <- coronary
system.time(approx_best <- search_greedy(data,alg = simple_marginal, join_unobserved = TRUE))
system.time(best <- search_all(data,alg = simple_marginal, join_unobserved = FALSE))

plot(approx_best)
plot(best)

BIC(best, approx_best)
