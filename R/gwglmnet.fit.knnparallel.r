gwglmnet.fit.knnparallel = function(x, y, family, coords, fit.loc, oracle, D, verbose, varselect.method, prior.weights, tuning, predict, simulation, gweight, target, beta1, beta2, tol.loc, longlat=FALSE, N, interact, resid.type) {
  if (!is.null(fit.loc)) {
    coords.unique = unique(fit.loc)
  } else {
    coords.unique = unique(coords)
  }
  n = dim(coords.unique)[1]

  gwglmnet.object = list()

  prior.weights = drop(prior.weights)
  max.weights = rep(1, length(prior.weights))
  total.weight = sum(max.weights * prior.weights)
  if (is.null(tol.loc)) {tol.loc = target / 1000}
  
  models = foreach(i=1:n, .packages=c('SGL'), .errorhandling='remove') %dopar% {
    loc = coords.unique[i,]
    dist = drop(D[i,])

    opt = optimize(gwglmnet.knn, lower=beta1, upper=beta2, 
      maximum=FALSE, tol=tol.loc, coords=coords, loc=loc,
      gweight=gweight, verbose=verbose, dist=dist, total.weight=total.weight,
      prior.weights=prior.weights, target=target)
    bandwidth = opt$minimum
    kernweights = drop(gweight(dist, bandwidth))
    
    if (is.null(oracle)) {
      m = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, N=N, varselect.method=varselect.method, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=kernweights, prior.weights=prior.weights, interact=interact)
    } else {
      m = gwselect.fit.oracle(x=x, y=y, family=family, coords=coords, loc=loc, N=N, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=kernweights, prior.weights=prior.weights, interact=interact, oracle=oracle[[i]])
    }
    
    if (verbose) {
      cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bandwidth,3), "; s=", m[['s']], "; sigma2=", round(tail(m[['sigma2']],1),3), "; nonzero=", paste(m[['nonzero']], collapse=","), "; weightsum=", round(m[['weightsum']],3), ".\n", sep=''))
    }
    return(m)
  }

  gwglmnet.object[['models']] = models

	if (tuning) { }
  else if (predict) { }
  else if (simulation) { }
  else {gwglmnet.object[['coords']] = coords}

  class(gwglmnet.object) = 'gwglmnet.object'
  return(gwglmnet.object)
}
