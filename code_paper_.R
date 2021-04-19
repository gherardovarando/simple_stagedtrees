library(stagedtrees)
library(bnlearn)
library(igraph)
library(magrittr)
library(arules)
library(deal)
source("functions.R")
source("search_order.R")

titanic.df <- as.data.frame(Titanic)
titanic.df <- titanic.df[rep(row.names(titanic.df), titanic.df$Freq), 1:4]

load("cachexia_data.RData")
load("chds.RData")

load("mathmarks.RData")
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
  fall = readRDS("FallEld.rds"),
  ksl = ksl,
  mathmarks = mathmarks,
  phd = stagedtrees::PhDArticles,
  titanic = titanic.df
)


experiment <- function(data, lambda = 0, r_train = 1, seed = 0){
  set.seed(seed)
  ix <- sample(seq(nrow(data)), size = nrow(data)*r_train)
  train <- data[ix,]
  test <-  data[-ix,]
  print(colnames(data))
  time.dag <- system.time(dag_hc <- bnlearn::hc(train))
  dag_fitted <- bn.fit(dag_hc, data = train)
  sevt_from_dag <- sevt_fit(as_sevt(dag_fitted), data = train, lambda = lambda)
  order <- bnlearn::node.ordering(dag_hc)
  time.simple_marginal <- system.time(marginal <- simple_marginal(full(train, order, join_unobserved = FALSE, lambda = lambda)))
  time.simple_total <- system.time(total <- simple_total_bhc(full(train, order, join_unobserved = FALSE, lambda = lambda)))
  time.sevt_bhc <- system.time(sevt_bhc <- stages_bhc(stndnaming(stages_bhc(full(train, order, join_unobserved = TRUE, lambda = lambda)), ignore = FALSE)))
  time.simplify <- system.time(simplified <- simplify(sevt_bhc))
  time.greedy_marginal <- system.time(greedy_marginal <- search_greedy(data = train,
                                                                       alg = simple_marginal, join_unobserved = FALSE, lambda = lambda))
  if (ncol(data)<7){
    time.all_marginal <- system.time(all_marginal <- search_all(train, alg = simple_marginal, join_unobserved = FALSE, lambda = lambda)) 
    time.all_total <- system.time(all_total <- search_all(train, alg = simple_total_bhc, join_unobserved = FALSE, lambda = lambda))
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
              train = train,
              test = test))
}


######## comparisons 

results <- lapply(datasets, experiment, lambda = 0, r_train = 1)

table.BIC <- t(as.data.frame(lapply(results, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    BIC(m) 
}))))
View(table.BIC)
saveRDS(table.BIC, file = "tableBIC.rds")

table.time <- t(as.data.frame(lapply(results, function(x) sapply(x$time, function(m) {
  if (length(m) == 1) return(NA) else 
    m[3]
}))))

saveRDS(table.time, file = "tabletime.rds")

table.positions <- t(as.data.frame(lapply(results, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    num_positions(m)
}))))
saveRDS(table.positions, file = "tablepositions.rds")
View(table.positions)

table.stages <- t(as.data.frame(lapply(results, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    num_stages(m)
}))))
saveRDS(table.stages, file = "tablestages.rds")
View(table.stages)

########## train-test

results_1 <- lapply(datasets, experiment, lambda = 1, r_train = 0.5)

table.logLik <- t(as.data.frame(lapply(results_1, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    sum(stagedtrees::prob(m, x$test, log = TRUE)) 
}))))
saveRDS(table.logLik, file = "tableloglik.rds")
View(table.logLik)


## CHDS PLOTS
plot(results$chds$dag)
plot(results$chds$models$simple_total)
ceg.plot(results$chds$models$simple_total)

## CORONARY PLOTS (TO DO)
plot(results$coronary$dag)
plot(results$coronary$models$simple_total)
ceg.plot(results$coronary$models$simple_total)

