gwglmnet.ssr = function(bw, x, y, family, coords, loc, dist, verbose, prior.weights, gweight, target, mode.select, interact, alpha, oracle, shrunk.fit, resid.type) {
    if (is.null(oracle)) {
        gwglmnet.object = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bw, dist=dist, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, mode.select=mode.select, interact=interact, alpha=alpha, shrunk.fit=shrunk.fit, predict=TRUE, tuning=FALSE, simulation=FALSE)
    }
    else {
        gwglmnet.object = gwselect.fit.oracle(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=NULL, interact=interact, prior.weights=prior.weights, gweight=gweight)
    }
    
    loss = gwglmnet.object[['ssr']][[resid.type]]
    if (verbose) {cat(paste('loc:(', paste(round(loc,3), collapse=","), '), target: ', round(target,3), ', bw:', round(bw,3), ', ssr:', round(loss,3), ', miss:', round(abs(loss-target),3), '\n', sep=""))}
    return(abs(loss-target))
}
