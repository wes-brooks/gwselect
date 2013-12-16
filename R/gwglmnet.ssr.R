gwglmnet.ssr = function(bw, x, y, family, coords, loc, dist, s, verbose, prior.weights, gweight, target, adapt, mode.select, precondition, interact, alpha, oracle, shrunk.fit, AICc) {
    if (is.null(oracle)) {
        gwglmnet.object = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bw, dist=dist, s=s, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, mode.select=mode.select, precondition=precondition, interact=interact, alpha=alpha, shrunk.fit=shrunk.fit, AICc=AICc)
    }
    else {
        gwglmnet.object = gwselect.fit.oracle(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=NULL, interact=interact, prior.weights=prior.weights, gweight=gweight, AICc=AICc)
    }

    loss = gwglmnet.object[['loss']]
    if (verbose) {cat(paste('loc:(', paste(loc, collapse=","), '), target: ', target, ', bw:', bw, ', ssr:', loss, ', miss:', abs(loss-target), '\n', sep=""))}
    return(abs(loss-target))
}
