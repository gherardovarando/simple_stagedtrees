search_all_order <- function(data, alg = stages_bhc, start = full, score = BIC){
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
  #return(allres[[which.min(scores)]]()
  return(allres)
}

