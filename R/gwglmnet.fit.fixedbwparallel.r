gwglmnet.fit.fixedbwparallel = function(x, y, family, coords, fit.loc=NULL, bw, D=NULL, oracle, verbose=FALSE, mode.select, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat, simulation, tuning, predict, interact, N, resid.type) {
  if (!is.null(fit.loc)) { coords.unique = fit.loc }
  else { coords.unique = unique(coords) }
  n = dim(coords.unique)[1]

  gwglmnet.object = list()
  models = list()

  if (is.null(gwr.weights)) {
    gwr.weights = gweight(D, bw)    
  }       

  gweights = list()
  for (j in 1:nrow(gwr.weights)) {
      gweights[[j]] = as.vector(gwr.weights[j,])
  }      

  models = foreach(i=1:n, .packages=c('SGL'), .errorhandling='remove') %dopar% {
    #Fit one location's model here
    loc = coords.unique[i,]
    gw = drop(gweights[[i]])
	
  	if (is.null(oracle)) {
      m = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, N=N, mode.select=mode.select, gwr.weights=gw, prior.weights=prior.weights, predict=predict, tuning=tuning, simulation=simulation, verbose=verbose, interact=interact)
  	} else {
      m = gwselect.fit.oracle(x=x, y=y, family=family, coords=coords, loc=loc, N=N, mode.select=mode.select, gwr.weights=gw, prior.weights=prior.weights, predict=predict, tuning=tuning, simulation=simulation, verbose=verbose, interact=interact, oracle=oracle[[i]])
    }
      
    if (verbose) {
      cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bw,3), "; s=", m[['s']], "; sigma2=", round(tail(m[['sigma2']],1),3), "; nonzero=", paste(m[['nonzero']], collapse=","), "; weightsum=", round(m[['weightsum']],3), ".\n", sep=''))
    }
    return(m)
  }

  gwglmnet.object[['models']] = models
    
	if (tuning) { }
  else if (predict) { }
  else if (simulation) { }
  else { gwglmnet.object[['coords']] = coords }
	
  class(gwglmnet.object) = 'gwglmnet.object'
  return(gwglmnet.object)
}
