load("chds.RData")
source("functions.R")
source("search_order.R")

alls <- all_orders(PhDArticles, alg = simple_total_bhc,
                        start = function(data, order){full(data, order, join_unobserved = FALSE, lambda = 0)})

bics <- sapply(alls, BIC)
plot(bics)
best <- alls[[which.min(bics)]]
plot(best)
as_parentslist(best)
ceg.plot(best)
BIC(best)
best4s <- alls[which(bics < 2830)]
length(best4s)
best4s
lapply(best4s, as_parentslist)

######## we check bn 
bn.hc <- hc(chds)
bn.pc <- pc.stable(chds, alpha = 0.01, test = "x2")
plot(bn.hc)
plot(bn.pc)
plot(bnlearn::cpdag(bn.hc))
plot(bnlearn::cpdag(bn.pc))
########

alls_bhc <- all_orders(coronary, alg = stages_hc,
                   start = function(data, order){full(data, order, join_unobserved = TRUE, lambda = 0)})

bics_hc <- sapply(alls_hc, BIC)
plot(bics_hc)

best_hc <- alls_hc[[which.min(bics_hc)]]
plot(best_hc)
as_parentslist(best_hc)
plot(alls_hc[[18]])


##########

alls_cs <- all_orders(chds, alg = stages_csbhc,
                      start = function(data, order){full(data, order, join_unobserved = FALSE, lambda = 0)})
bics_cs <- sapply(alls_cs, BIC)
plot(bics_cs)

best_cs <- alls_hc[[which.min(bics_cs)]]
plot(best_cs)
as_parentslist(best_cs)

bests_cs <- alls_hc[bics_cs < 2827]
length(bests_cs)

lapply(bests_cs, logLik)

plot(bests_cs[[4]])


#### 

approx_best <- search_greedy(chds)
alls <- all_orders(chds, alg = stages_bhc)

bics <- sapply(alls, BIC)
plot(bics)
points(BIC(approx_best), col = "red")
best <- alls[[which.min(bics)]]
BIC(best)
BIC(approx_best)
