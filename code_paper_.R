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




######## comparisons 

results <- lapply(datasets, function(data){
  print(colnames(data))
  time.dag <- system.time(dag_hc <- bnlearn::hc(data))
  dag_fitted <- bn.fit(dag_hc, data = data)
  sevt_from_dag <- sevt_fit(as_sevt(dag_fitted), data = data, lambda =0)
  order <- bnlearn::node.ordering(dag_hc)
  time.simple_marginal <- system.time(marginal <- simple_marginal(full(data, order, join_unobserved = TRUE)))
  time.simple_total <- system.time(total <- simple_total_bhc(full(data, order, join_unobserved = TRUE)))
  time.sevt_bhc <- system.time(sevt_bhc <- stages_bhc(full(data, order, join_unobserved = TRUE)))
  time.simplify <- system.time(simplified <- simplify(sevt_bhc))
  time.greedy_marginal <- system.time(greedy_marginal <- search_greedy(data = data,
                                                                       alg = simple_marginal))
  if (ncol(data)<7){
    time.all_marginal <- system.time(all_marginal <- search_all(data, alg = simple_marginal)) 
    time.all_total <- system.time(all_total <- search_all(data, alg = simple_total_bhc))
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
              dag = dag_hc))
})

table.BIC <- t(as.data.frame(lapply(results, function(x) sapply(x$models, function(m) {
  if (length(m) == 1) return(NA) else 
    BIC(m) 
}))))

saveRDS(table.BIC, file = "tableBIC.rds")


table.time <- t(as.data.frame(lapply(results, function(x) sapply(x$time, function(m) {
  if (length(m) == 1) return(NA) else 
    m[3]
}))))

saveRDS(table.time, file = "tabletime.rds")


## CHDS PLOTS
plot(results$chds$dag)
plot(results$chds$models$simple_total)
ceg.plot(results$chds$models$simple_total)

## CORONARY PLOTS (TO DO)
plot(dag_cor)
plot(total_cor)
ceg.plot(total_cor)
