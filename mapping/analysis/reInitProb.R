reInitProb = function(arena){
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  arena@specs = lapply(arena@specs, function(x){
    x@lpobj = sybil::sysBiolAlg(x@model, algorithm="fba")
    x@fbasol <- sybil::optimizeProb(x@lpobj)
    names(x@fbasol$fluxes) = x@model@react_id
    return(x)
  })
  return(arena)
}