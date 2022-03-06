library("stagedtrees")
library("bnlearn")
library("igraph")
library("arules")
library("deal")
library("xtable")

source("functions.R")

################ define experiment function ################
experiment <- function(data, lambda = 0, verbose = FALSE){
  if (verbose) print(colnames(data))
  if (verbose) print(ncol(data))
  time.dag <- system.time(dag_hc <- bnlearn::hc(data))
  if (verbose) message("dag done")
  dag_fitted <- bn.fit(dag_hc, data = data)
  if (verbose) message("dag fit")
  sevt_from_dag <- sevt_fit(as_sevt(dag_fitted), data = data, 
                            lambda = lambda)
  if (verbose) message("sevt from dag done")
  order <- bnlearn::node.ordering(dag_hc)
  model0 <- full(data, order, join_unobserved = FALSE, lambda = lambda)
  if (verbose) message("full staged tree done")
  time.simple_marginal <- system.time(marginal <- simple_marginal(model0))
  if (verbose) message("marginal done")
  time.simple_total <- system.time(total <- simple_total_bhc(model0))
  if (verbose) message("total done")
  time.sevt_bhc <- system.time(sevt_bhc <- stages_bhc(model0, ignore = FALSE))
  if (verbose) message("bhc done")
  time.simplify <- system.time(simplified <- simplify(sevt_bhc))
  if (verbose) message("simplify done")
  time.greedy_marginal <- system.time(greedy_marginal <- search_greedy(data = data,
                                                                       alg = simple_marginal, 
                                                                       join_unobserved = FALSE, 
                                                                       lambda = lambda))
  if (verbose) message("greedy done")
  if (ncol(data)<7){
    time.all_marginal <- system.time(all_marginal <- search_all(data, alg = simple_marginal, 
                                                                join_unobserved = FALSE, lambda = lambda))
    if (verbose) message("all marginal done")
    time.all_total <- system.time(all_total <- search_all(data, alg = simple_total_bhc, 
                                                          join_unobserved = FALSE, lambda = lambda))
    if (verbose) message("all total done")
  }else{
    time.all_marginal <- time.all_total <- NA
    all_marginal <- all_total <- NA
  }
  return(list(models = list(dag = sevt_from_dag, 
                            simple_marginal = marginal, 
                            simple_total = total,
                            sevt_bhc = sevt_bhc,
                            simplified = simplified,
                            greedy_marginal = greedy_marginal,
                            all_marginal = all_marginal,
                            all_total = all_total), 
              time = list(
                dag = time.dag,
                simple_marginal = time.simple_marginal,
                simple_total = time.simple_total,
                sevt_bhc = time.sevt_bhc,
                simplified = time.simplify,
                greedy_marginal = time.greedy_marginal,
                all_marginal = time.all_marginal,
                all_total = time.all_total
              ),
              dag = dag_hc,
              data = data))
}

## change to TRUE to save results
save <- TRUE
save.latex <- TRUE

########## load data ############
titanic.df <- as.data.frame(Titanic)
titanic.df <- titanic.df[rep(row.names(titanic.df), titanic.df$Freq), 1:4]

load("data/cachexia_data.RData")
load("data/chds.RData")

load("data/mathmarks.RData")
mathmarks <- discretizeDF(mathmarks)

data("ksl")
ksl[,1] <- discretize(ksl[,1],breaks=2)
ksl[,2] <- discretize(ksl[,2],breaks=2)
ksl[,3] <- as.factor(ksl[,3])
ksl[,4] <- discretize(ksl[,4],breaks=2)

datasets <- list(
  asia = bnlearn::asia,
  cachexia = cbind(discretizeDF(cachexia_data[,-7]),MC = as.factor(cachexia_data$MC)),
  chds = chds,
  coronary = bnlearn::coronary,
  fall = readRDS("data/FallEld.rds"),
  ksl = ksl,
  mathmarks = mathmarks,
  phd = stagedtrees::PhDArticles,
  titanic = titanic.df
)


######## run experiments
results <- lapply(datasets, experiment, lambda = 0)


######## collect results and print tables
prettynames <- c("BN", "Marg.", "Total", "BHC", "Simpl.", "Greedy Marg.", "All Marg.", "All Total")


table.BIC <- t(as.data.frame(lapply(results, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    BIC(m) 
}))))


if (save) saveRDS(table.BIC, file = "tableBIC.rds")
if  (save.latex){
  colnames(table.BIC) <- prettynames
  tmp <- xtable(table.BIC, 
                align = c("l", rep("r", 8)),
                label = "table:bic",
         caption = "BIC for models learned over nine datasets: 
                    Bayesian networks (BN), 
                    simple staged trees with marginal algorithm (Marg.), 
                    simple staged trees with total algorithm (Total), 
                    generic staged trees with backward hill-climbing (BHC), 
                    simple staged trees derived via simplification (Simpl.), 
                    simple staged tree with variable ordering (Greedy Marg.), 
                    simple staged tree considering all orders 
                    and marginal algorithm (All-Marg.), 
                    simple staged tree considering all orders and total algorithm (All-Tot.).")
  print(tmp, booktabs = TRUE, file = "bic.table.tex", scalebox = 0.77)
}

table.time <- t(as.data.frame(lapply(results, function(x) sapply(x$time, function(m) {
  if (length(m) == 1) return(NA) else 
    m[3]
}))))

if (save) saveRDS(table.time, file = "tabletime.rds")
if  (save.latex){
  colnames(table.time) <- prettynames
  tmp <- xtable(table.time,
                align = c("l", rep("r", 8)),
                digits = 3,
                label = "table:time",
                caption = "Time (in seconds) to learn the models reported 
                in Table \\ref{table:bic}. 
                Computations carried out with an intel CORE i7 vPro 8th Gen.")
  print(tmp, booktabs = TRUE, file = "time.table.tex", scalebox = 0.77)
}

table.positions <- t(as.data.frame(lapply(results, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    num_positions(m)
}))))
if (save) saveRDS(table.positions, file = "tablepositions.rds")
if  (save.latex){
  colnames(table.positions) <- prettynames
  tmp <- xtable(table.positions, 
                digits = 0,
                align = c("l", rep("r", 8)),
                label = "table:pos",
                caption = "Number of positions for the staged trees reported in 
                Table \\ref{table:bic}. (TODO in R) For BHC the number in parenthesis is the ratio 
                between the number of stages and the number of positions.")
  print(tmp, booktabs = TRUE, file = "positions.table.tex", scalebox = 0.77)
}

table.stages <- t(as.data.frame(lapply(results, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    num_stages(m)
}))))
if (save) saveRDS(table.stages, file = "tablestages.rds")
if  (save.latex){
  colnames(table.stages) <- prettynames
  tmp <- xtable(table.stages, 
                digits = 0,
                align = c("l", rep("r", 8)),
                label = "table:stg",
                caption = "Number of stages for the staged trees reported in 
                Table \\ref{table:bic}. For BHC the number in parenthesis is the ratio 
                between the number of stages and the number of positions.")
  print(tmp, booktabs = TRUE, file = "stages.table.tex", scalebox = 0.77)
}

########## train-test

results_1 <- lapply(datasets, experiment, lambda = 1)

table.logLik <- t(as.data.frame(lapply(results_1, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    sum(stagedtrees::prob(m, x$test, log = TRUE)) 
}))))
if (save) saveRDS(table.logLik, file = "tableloglik.rds")
if  (save.latex){
  colnames(table.loglik) <- prettynames
  tmp <- xtable(table.loglik,
                align = c("l", rep("r", 8)),
                digits = 3,
                label = "table:loglik",
                caption = "Predictive log-likelihood for the models reported 
                           in Table \\ref{table:bic}. The results are based on 
                           a single random split of the data into 50\\% training 
                           and 50\\% test")
  print(tmp, booktabs = TRUE, file = "time.table.tex", scalebox = 0.77)
}



## CHDS PLOTS
plot(results$chds$dag)
plot(results$chds$models$simple_total)
ceg.plot(results$chds$models$simple_total)

## CORONARY PLOTS (TO DO)
plot(results$coronary$dag)
plot(results$coronary$models$simple_total)
ceg.plot(results$coronary$models$simple_total)







