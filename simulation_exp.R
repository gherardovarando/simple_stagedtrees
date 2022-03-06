library("stagedtrees")
library("bnlearn")
library("pbapply")

source("functions.R")

################ define experiment function ################
experiment_fixorder <- function(data, lambda = 0, verbose = FALSE){
  model0 <- full(data, join_unobserved = FALSE, lambda = lambda)
  marginal <- simple_marginal(model0)
  total <- simple_total_bhc(model0)
  sevt_bhc <- stages_bhc(stndnaming(stages_bhc(model0, ignore = FALSE)))
  simplified <- simplify(sevt_bhc)
  return(list(simple_marginal = marginal, 
              simple_total = total,
              sevt_bhc = sevt_bhc,
              simplified = simplified))
}

dir.create("sim_results", showWarnings = FALSE)
pboptions(type = "txt")
############  simulation experiments 
M <- 100  # number of repetitions
Ntest <- 1000 # testing sample size
Ns <- c(25, 50, seq(from = 100, to = 1000, by = 100), 
        1250, 1500, 1750, 2000)
for (n in c(6)){ #number of variables
for (q in c(0.2, 0.5, 0.7)){# parameter that control number of stages/positions
for (N in Ns){
  message(q, "-", N)
  results <- t(pbreplicate(M, {
    true <- random_simple_sevt(n, q)
    data <- sample_from(true, nsim = N)
    res <- experiment_fixorder(data, lambda = 1)
    #test <- sample_from(true, nsim = Ntest)
    sapply(res, function(mod){
      htree <- hamming_stages(mod, true, return_tree = TRUE)
      sum(sapply(htree, function(x) {
        mean(as.numeric(x), na.rm = TRUE)
      }))
    })
  }, cl = 5))
  saveRDS(results, file = paste0("sim_results/", n, "_",
                                 N, "_", q, ".rds"))
}
}
}
