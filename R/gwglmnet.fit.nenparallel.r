gwglmnet.fit.nenparallel = function(x, y, family, coords, fit.loc, D, varselect.method, verbose, prior.weights, gweight, target, beta1, beta2, tol.loc, longlat=FALSE, tuning, simulation, predict, N, oracle, interact, resid.type) {
  coords.unique = unique(coords)
  n = dim(coords.unique)[1]
  gwglmnet.object = list()
  models = list()

  if (verbose) {cat(paste('beta1:', round(beta1,3), ', beta2:', round(beta2,3), '\n', sep=''))}
  if (is.null(tol.loc)) {tol.loc = target / 1000}
  
  models = foreach(i=1:n, .packages=c('SGL'), .errorhandling='remove') %dopar% {
    loc = coords.unique[i,]
    dist = D[i,]

    #Compute the bandwidth to achieve total local weighted squared error equal to the target:
    opt = optimize(gwglmnet.ssr, lower=beta1, upper=beta2, 
      maximum=FALSE, tol=tol.loc, x=x, y=y, coords=coords, loc=loc, resid.type=resid.type,
      gweight=gweight, verbose=verbose, dist=dist, varselect.method=varselect.method, family=family, oracle=oracle,
      prior.weights=prior.weights, target=target, interact=interact)
    bandwidth = opt$minimum
    kernweights = drop(gweight(dist, bandwdth))

	  if (verbose) {
    	cat(paste("For i=", i, ", target:", round(target,3), ", bw=", round(bandwidth,3), ", tolerance=", round(target/1000,3), ", miss=", round(opt$objective,3), ".\n", sep=''))
    }

    #Fit the local model:
    if (is.null(oracle)) {
      m = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, N=N, varselect.method=varselect.method, verbose=verbose, gwr.weights=kernweights, prior.weights=prior.weights, tuning=tuning, predict=predict, simulation=simulation, interact=interact)
    } else {
      m = gwselect.fit.oracle(x=x, y=y, family=family, coords=coords, loc=loc, N=N, verbose=verbose, gwr.weights=kernweights, prior.weights=prior.weights, tuning=tuning, predict=predict, simulation=simulation, interact=interact, oracle=oracle[[i]])
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
