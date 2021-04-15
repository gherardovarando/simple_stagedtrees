library(stagedtrees)
library(igraph)

simplify <- function(model){
  tmp_ceg <- stagedtrees::ceg(model)
  model$stages <- tmp_ceg$positions[-1]
  return(stagedtrees::sevt_fit(model))
}

simple_marginal <- function(data, lambda = 1, order = NULL,
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
  if (is.null(order)){
    if (is.table(data)){
      order <- names(dimnames(data))
    }else if(is.data.frame(data)){
      order <- colnames(data)
    }
  }
  model <- sevt(data, order = order, full = TRUE)
  order <- names(model$tree)
  for (i in 2:length(order)){
    v <- order[i]
    lv <- length(model$tree[[order[i-1]]])
    #os <- model$stages[[v]] == model$name_unobserved[1]
    if (i > 2) model$stages[[v]] <- 
      vapply(model$stages[[order[i-1]]], function(s){
        paste0(s, 1:lv)
      }, FUN.VALUE = rep("1", lv))[TRUE]
    #model$stages[[v]][os] <- model$name_unobserved[1]
    model <- sevt_fit(model, data = data, lambda = lambda)
    model <- alg(model, scope = v, ...) 
  }
  model <- stndnaming(model)
  return(model)
}

join_positions <- function(model, v, s1, s2){
  i <- which(v == names(model$tree))
  order <- names(model$tree)
  model <- join_stages(model, v, s1, s2)
  if (i == length(model$tree)) return(model)
  for (j in (i+1):length(model$tree)){
    w <- names(model$tree)[j]
    lv <- length(model$tree[[j-1]])
    model$stages[[w]] <- 
      vapply(model$stages[[order[j-1]]], function(s){
        paste0(s, 1:lv)
      }, FUN.VALUE = rep("1", lv))[TRUE]
  }
  return(sevt_fit(model))
}

simple_total_bhc <- function (object, score = function(x) {
  return(-BIC(x))
}, max_iter = Inf) {
  now_score <- score(object)
  scope <- names(object$tree)[-1]
  for (v in scope) {
    r <- 1
    iter <- 0
    done <- FALSE
    while (!done && iter < max_iter) {
      iter <- iter + 1
      temp <- object
      temp_score <- now_score
      done <- TRUE
      stages <- unique(object$stages[[v]])
      if (length(stages) > 1) {
        for (i in 2:length(stages)) {
          s1 <- stages[i]
          for (j in 1:(i - 1)) {
            s2 <- stages[j]
            try <- join_positions(object, v, s1, s2)
            try_score <- score(try)
            if (try_score >= temp_score) {
              temp <- try
              temp_score <- try_score
              s1a <- s1
              s2a <- s2
              done <- FALSE
            }
          }
        }
      }
      object <- temp
      now_score <- temp_score
    }
  }
  object$call <- sys.call()
  object$score <- list(value = now_score, f = score)
  return(object)
}

ceg.plot <- function(tree){
  strut <- ceg(tree)
  A <- ceg2adjmat(strut)
  gr <- graph_from_adjacency_matrix(A)
  plot.igraph(gr, layout = layout_as_tree, edge.arrow.size = 0.2)
}

num_positions <- function(tree){
  ceg <- ceg(tree)
  value <- 0
  for(i in 1:length(ceg$positions)){
    value <- value + length(unique(ceg$positions[[i]]))
  }
  value
}
num_stages <- function(tree){
  value <- 1
  for(i in 1:length(tree$stages)){
    value <- value + length(unique(tree$stages[[i]]))
  }
  value
}
