library(stagedtrees)
library(igraph)

simplify <- function(model){
  tmp_ceg <- stagedtrees::ceg(model)
  model$stages <- tmp_ceg$positions[-1]
  return(stagedtrees::sevt_fit(model))
}

simple_marginal <- function(object,
                            alg = stages_bhc,
                            scope = NULL,
                            ...){
  
  if (is.null(scope)){
    scope <- names(object$tree)[-1]
  }

  order <- names(object$tree)
  for (v in scope){
    i <- seq_along(names(object$tree))[names(object$tree) == v]
    lv <- length(object$tree[[order[i-1]]])
    os <- object$stages[[v]] %in% object$name_unobserved
    if (i > 2) object$stages[[v]] <- 
      vapply(object$stages[[order[i-1]]], function(s){
        paste0(s, 1:lv)
      }, FUN.VALUE = rep("1", lv))[TRUE]
    object$stages[[v]][os] <- object$name_unobserved[1]
    object <- sevt_fit(object, lambda = object$lambda)
    object <- alg(object, scope = v, ...) 
  }
  object <- stndnaming(object)
  return(object)
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

random_sevt <- function(n, k = 2){
  tree <- sapply(paste0("X", seq(n)), function(x) c("0","1"), 
                 USE.NAMES = TRUE, simplify = FALSE)
  model <- sevt(tree)
  model$stages <- lapply(model$stages, FUN = function(stages){
    paste0(sample(1:k, size = length(stages), replace = TRUE))
  })
  model$prob <- list()
  model$prob <- lapply(model$stages, function(stages){
    sapply(unique(stages), FUN = function(s){
      p <- runif(2) 
      p <- p / sum(p)
      names(p) <- c("0", "1")
      attr(p, "n") <- 1
      return(p)
    }, simplify = FALSE, USE.NAMES = TRUE)
  })
  p <- runif(2) 
  p <- p / sum(p)
  names(p) <- c("0", "1")
  attr(p, "n") <- 1
  model$prob <- c(list("X1" = list("1" = p)), model$prob)
  return(model)
}

random_simple_sevt <- function(n, q = 0.5){
  tree <- sapply(paste0("X", seq(n)), function(x) c("0","1"), 
                 USE.NAMES = TRUE, simplify = FALSE)
  model <- sevt(tree)
  model$stages$X2 <- sample(c("1", "2"), replace = TRUE)
  for (i in seq_along(tree)[-1]){
    v <- names(tree)[i]
    lv <- length(tree[[i-1]])
    if (i > 2) model$stages[[v]] <- 
      vapply(model$stages[[names(tree)[i-1]]], function(s){
        paste0(s, 1:lv)
      }, FUN.VALUE = rep("1", lv))[TRUE]
    for (s in unique(model$stages[[v]])){
      if (runif(1) < q){
        model$stages[[v]][model$stages[[v]] == s] <- sample(unique(model$stages[[v]]), size = 1) 
      }
    }
  }
  model$prob <- list()
  model$prob <- lapply(model$stages, function(stages){
    sapply(unique(stages), FUN = function(s){
      p <- c(runif(1, max = 1),runif(1, max = 5)) 
      p <- sample(p) / sum(p)
      names(p) <- c("0", "1")
      attr(p, "n") <- 1
      return(p)
    }, simplify = FALSE, USE.NAMES = TRUE)
  })
  p <- runif(2) 
  p <- p / sum(p)
  names(p) <- c("0", "1")
  attr(p, "n") <- 1
  model$prob <- c(list("X1" = list("1" = p)), model$prob)
  return(model)
}
