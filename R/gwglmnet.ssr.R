gwglmnet.ssr = function(bw, x, y, family, coords, loc, dist, verbose, prior.weights, gweight, target, mode.select, interact, oracle, resid.type) {
  #Calculate the local weights:
  kernweights = drop(gweight(dist, bw))
  
  if (is.null(oracle)) {
    gwglmnet.object = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, mode.select=mode.select, interact=interact, predict=TRUE, tuning=FALSE, simulation=FALSE, verbose=verbose, gwr.weights=kernweights, prior.weights=prior.weights)
  } else {
    gwglmnet.object = gwselect.fit.oracle(x=x, y=y, family=family, coords=coords, loc=loc, mode.select=mode.select, interact=interact, predict=TRUE, tuning=FALSE, simulation=FALSE, verbose=verbose, gwr.weights=kernweights, prior.weights=prior.weights, oracle=oracle[[i]])
  }
  
  loss = gwglmnet.object[['ssr']][[resid.type]]
  if (verbose) {cat(paste('loc:(', paste(round(loc,3), collapse=","), '), target: ', round(target,3), ', bw:', round(bw,3), ', ssr:', round(loss,3), ', miss:', round(abs(loss-target),3), '\n', sep=""))}
  return(abs(loss-target))
}
