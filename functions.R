library(stagedtrees)
library(igraph)

#' Simplify a stagedtree 
#' 
#' Function to simplify a staged tree model. 
#' @param model an object of class \code{sevt}
#' @return an object of class \code{sevt} 
#'         representing the simplified model. 
#'         The returned model will be fitted if the input \code{model} was. 
#' @details The \code{simplify} function will produce the corresponding simple
#'          staged tree, that is a staged tree where stages and positions are 
#'          equivalent. 
#'          To do so the function \code{ceg} is used to compute positions, and 
#'          then the stages' vectors are replaced with the positions' vectors. 
#'          The model is the re-fitted if the input was a fitted staged tree. 
#'          
#'          Despite the name, the simplified staged tree has always a number 
#'          of stages greater or equal to the initial staged tree, thus it is 
#'          a more complex statistical model. 
simplify <- function(model){
  tmp_ceg <- stagedtrees::ceg(model)
  model$stages <- tmp_ceg$positions[-1]
  ## if model was fitted then refit it
  if (stagedtrees:::is_fitted_sevt(model)){
    model <- stagedtrees::sevt_fit(model)
  }
  return(model)
}

#' Marginal search of simple staged trees
#' 
#' @param object an object of class \code{sevt}.
#' @param alg an algorithm to search stages, should accept the \code{scope} argument. 
#' @param scope NULL.
#' @param ... additional parameters passed to \code{alg}.
#' @return an object of class \code{sevt}, the simple staged tree resulting 
#'         from the search. 
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


#' join positions in a staged tree model
#' 
#' @param model an object of class \code{sevt}.
#' @param v the name of a variable in the model.
#' @param s1 stage to join
#' @param s2 stage to join
#' @details this functions works similarly to the \code{join_stages}
#' function in the \code{stagedtrees} package, but it also joins 
#' downstream stages to make nodes with stages \code{s1,s2} in the same 
#' position. This function works properly only when downstream variables 
#' from \code{v} have full stages vectors. 
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

#' Total BHC search of simple staged trees
#' 
#' @param object an object of class \code{sevt}.
#' @param score a scroe to be maximized in the search.
#' @return an object of class \code{sevt}, the simple staged tree resulting 
#'         from the search. 
simple_total_bhc <- function (object, 
                              score = function(x) {return(-BIC(x))}, 
                              max_iter = Inf) {
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

#' obtain number of positions 
num_positions <- function(tree){
  ceg <- ceg(tree)
  value <- 0
  for(i in 1:length(ceg$positions)){
    value <- value + length(unique(ceg$positions[[i]]))
  }
  value
}

#' obtain number of stages
num_stages <- function(tree){
  value <- 1
  for(i in 1:length(tree$stages)){
    value <- value + length(unique(tree$stages[[i]]))
  }
  value
}

#' generate a random staged tree
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

#' generate a random simple staged tree
random_simple_sevt <- function(n, q = 0.5){
  tree <- sapply(paste0("X", seq(n)), function(x) c("a","b"), 
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
      p <- runif(2) 
      p <- p / sum(p)
      names(p) <- c("a", "b")
      attr(p, "n") <- 1
      return(p)
    }, simplify = FALSE, USE.NAMES = TRUE)
  })
  p <- runif(2) 
  p <- p / sum(p)
  names(p) <- c("a", "b")
  attr(p, "n") <- 1
  model$prob <- c(list("X1" = list("1" = p)), model$prob)
  return(model)
}


search_all <- function(data, alg = stages_bhc, 
                       search_score = BIC,
                       join_unobserved = TRUE, lambda = 0, ...){
  if (is.data.frame(data)){
    vv <- colnames(data)
  }
  if (is.table(data)){
    vv <- names(dimnames(data))
  }
  allres <- combinat::permn(vv, fun = function(order){
    alg(full(data, order =  order, 
             join_unobserved = join_unobserved, lambda = lambda))
  })    
  scores <- sapply(allres, search_score)
  return(allres[[which.min(scores)]])
}


