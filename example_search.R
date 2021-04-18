data <- chds
system.time(best1 <- search_all(data, join_unobserved = TRUE))

system.time(best2 <- search_dynamic(data, join_unobserved = TRUE))

system.time(approx_best <- search_greedy(data, join_unobserved = TRUE))
best1
best2
approx_best

BIC(best1, best2, approx_best)
