
library(igraph)
library(stagedtrees)
library(magrittr)


### 
nested_bhc <- function(data, lambda = 1, order = colnames(data),...){
  model <- full(data, lambda = lambda, order = order, ...)
  for (i in 2:length(order)){
    v <- order[i]
    lv <- length(model$tree[[order[i-1]]])
    os <- model$stages[[v]] == model$name_unobserved[1]
    if (i > 2) model$stages[[v]] <- 
        vapply(model$stages[[order[i-1]]], function(s){
          paste0(s, 1:lv)
        }, FUN.VALUE = rep("1", lv))[TRUE]
    model$stages[[v]][os] <- model$name_unobserved[1]
    print(model$stages[[v]])
    model <- sevt_fit(model, lambda = lambda)
    model <- stages_bhc(model, scope = v) 
  }
  model <- stndnaming(model)
  return(model)
}


## change to path of FallEld.rds 
fall <- readRDS(file = "FallEld.rds")
order <- c("A", "R", "T", "answer")

nested <- nested_bhc(data = fall, order = order) %>% stndnaming(uniq = TRUE) 

png(file = "fall-plot-nested.png", width = 1000, height = 2000, res = 200)
plot(nested)
dev.off()

ceg <- ceg(nested)
Adj <- ceg2adjmat(ceg)
g <- graph_from_adjacency_matrix(Adj)
plot(g, layout = layout_as_tree, edge.arrow.size = 0.2)
