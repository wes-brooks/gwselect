gwglmnet.fit.fixedbw = function(x, y, family, coords, fit.loc=NULL, oracle, bw, D=NULL, verbose=FALSE, mode.select, gwr.weights=NULL, indx=NULL, prior.weights=NULL, tuning=FALSE, predict=FALSE, gweight=NULL, longlat=FALSE, alpha, simulation, interact, N, shrunk.fit, resid.type) {
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
    gweights[[j]] = drop(gwr.weights[j,])
  }

  for (i in 1:n) {
    #Fit one location's model here
    loc = coords.unique[i,]
    gw = gweights[[i]]

  	if (is.null(oracle)) {
      models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, verbose=verbose, N=N, mode.select=mode.select, gwr.weights=gw, prior.weights=prior.weights, predict=predict, tuning=tuning, simulation=simulation, interact=interact)
  	} else {
      models[[i]] = gwselect.fit.oracle(x=x, y=y, family=family, coords=coords, loc=loc, verbose=verbose, N=N, mode.select=mode.select, gwr.weights=gw, prior.weights=prior.weights, predict=predict, tuning=tuning, simulation=simulation, interact=interact, oracle=oracle[[i]])
    }

    if (verbose) {
      cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bw, 3), "; loss=", round(tail(models[[i]][['tunelist']][['ssr-loc']][[resid.type]],1),3), "; s=", models[[i]][['s']], "; sigma2=", round(tail(models[[i]][['sigma2']],1),3), "; nonzero=", paste(models[[i]][['nonzero']], collapse=","), "; weightsum=", round(models[[i]][['weightsum']],3), ".\n", sep=''))
  	}
  }

  gwglmnet.object[['models']] = models
	
	if (tuning) { }
  else if (predict) { }
  else if (simulation) { }
  else { gwglmnet.object[['coords']] = coords }

  class(gwglmnet.object) = 'gwglmnet.object'
  return(gwglmnet.object)
}
