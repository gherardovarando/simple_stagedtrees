load("chds.RData")
source("functions.R")
source("search_order.R")

alls <- search_all_order(chds, alg = simple_total_bhc,
                                           score = BIC, start = function(data, order){full(data, order, join_unobserved = FALSE, lambda = 0)})

bics <- sapply(alls, BIC)
plot(bics)
alls[1]
alls[11]
alls[14]
alls[24]
total_chds <- simple_total_bhc(full(chds, order = c("Social","Economic","Events","Admission"), join_unobserved = F, lambda =0))

total_chds
BIC(total_chds, best_simple_total_chds)

ceg.plot(total_chds)
ceg.plot(best_simple_total_chds)
