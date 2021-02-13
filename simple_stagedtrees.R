library(igraph)
library(stagedtrees)
library(magrittr)


### 
simple_forward <- function(data, lambda = 1, order = colnames(data),
                       search = c("bhc", "fbhc", "bj", 'hclust', 'kmeans'),
                       ...){
  alg <- switch(search[1],
    bhc = stages_bhc,
    fbhc =  stages_fbhc,
    bj = stages_bj,
    hclust = stages_hclust,
    kmeans = stages_kmeans,
    NULL
  )
  model <- sevt(data, order = order)
  for (i in 2:length(order)){
    v <- order[i]
    lv <- length(model$tree[[order[i-1]]])
    os <- model$stages[[v]] == model$name_unobserved[1]
    if (i > 2) model$stages[[v]] <- 
        vapply(model$stages[[order[i-1]]], function(s){
          paste0(s, 1:lv)
        }, FUN.VALUE = rep("1", lv))[TRUE]
    model$stages[[v]][os] <- model$name_unobserved[1]
    #print(model$stages[[v]])
    model <- sevt_fit(model, data = data, lambda = lambda)
    model <- alg(model, scope = v, ...) 
  }
  model <- stndnaming(model)
  return(model)
}


## change to path of FallEld.rds 
fall_eld <- readRDS(file = "FallEld.rds")
order <- c("A", "R", "T", "answer")

model <- simple_forward(data = fall_eld, lambda = 1, order = order, search = "kmeans", k = 2) %>% stndnaming(uniq = TRUE) 

#png(file = "fall-plot-nested.png", width = 1000, height = 2000, res = 200)
plot(model)
#dev.off()

ceg <- ceg(model)
Adj <- ceg2adjmat(ceg)
g <- graph_from_adjacency_matrix(Adj)
plot(g, layout = layout_as_tree, edge.arrow.size = 0.2)
