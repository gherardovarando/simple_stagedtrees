all_orders <- function(data, alg = stages_bhc, start = full){
  if (is.data.frame(data)){
    vv <- colnames(data)
  }
  if (is.table(data)){
    vv <- names(dimnames(data))
  }
  allres <- combinat::permn(vv, fun = function(order){
    alg(start(data, order =  order))
  })    
  scores <- sapply(allres, BIC)
  return(allres)
}


sevt_add <- function(object, var, data){
  if (is.data.frame(data)) {
    data <- table(data[, c(names(object$tree), var)])
  }
  path <- names(object$tree)
  tt <- apply(data, MARGIN = c(path, var), sum)
  ctable <- ftable(tt, col.vars = var, row.vars = path)
  tmp <- sevt(data, full = TRUE, order = c(path, var))
  object$tree <- tmp$tree
  object$stages[[var]] <- tmp$stages[[var]]
  object$ctables[[var]] <- ctable
  object$prob[[var]] <- lapply(seq_along(object$stages[[var]]), function(ix){
    tt <- ctable[ix, ]
    names(tt) <- object$tree[[var]]
    n <- sum(tt)
    tt <- (tt + object$lambda)
    tt <- tt / sum(tt)
    tt[is.nan(tt)] <- NA
    attr(tt, "n") <- n
    return(tt)
  })
  names(object$prob[[var]]) <- object$stages[[var]]
  object$ll <- NULL ## force recompute log-likelihood
  object$ll <- logLik(object)
  return(object)
}

### should use github version in branch main
search_greedy <- function(data, alg = stages_bhc, score = BIC){
  if (is.data.frame(data)){
    vs <- colnames(data)
  }
  if (is.table(data)){
    vs <- names(dimnames(data))
  }
  ## initialize best
  best <- full(data, order = vs[1])
  ## check all other possible first variable
  if (length(vs) < 2) return(best)
  for (v in vs){
    tmp <- full(data, order = v)
    #print(score(tmp))
    if (score(tmp) < score(best)){
      best <- tmp
    }
  }
  object <- best
  ## add the best one by one 
  svs <- vs[!(vs %in% names(object$tree))]
  for (i in seq_along(vs)[-1]){
    #done <- FALSE
    best <- alg(sevt_add(object, svs[1], data), scope = svs[1])
    for (v in svs[-1]){
      tmp <- alg(sevt_add(object, v, data), scope = v)
      if (score(tmp) < score(best)){
        best <- tmp
        #done <- TRUE
      }  
    }
    #if (!done) break
    object <- best
    svs <- vs[!(vs %in% names(object$tree))]     
  }
  return(object)
}

