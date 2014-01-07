gwglmnet.fit.nen = function(x, y, family, coords, indx, fit.loc, D, mode.select, verbose, prior.weights, gweight, target, beta1, beta2, tol.loc, longlat=FALSE, precondition=FALSE, tuning, simulation, predict, N, oracle, interact, alpha, shrunk.fit, resid.type) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwglmnet.object = list()
    models = list()

    if (verbose) {cat(paste('beta1:', round(beta1,3), ', beta2:', round(beta2,3), '\n', sep=''))}
    if (is.null(tol.loc)) {tol.loc = target / 1000}
    
    for (i in 1:n) {
        loc = coords.unique[i,]
        dist = D[i,]

        opt = optimize(gwglmnet.ssr, lower=beta1, upper=beta2, maximum=FALSE, tol=tol.loc,
            x=x, y=y, coords=coords, loc=loc, alpha=alpha, resid.type=resid.type,
            gweight=gweight, verbose=verbose, dist=dist, mode.select=mode.select, family=family, oracle=oracle,
            prior.weights=prior.weights, target=target, precondition=precondition, interact=interact, shrunk.fit=shrunk.fit)
        bandwidth = opt$minimum

        models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bandwidth, dist=dist, mode.select=mode.select, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, tuning=tuning, simulation=simulation, predict=predict, precondition=precondition, N=N, interact=interact, shrunk.fit=shrunk.fit, alpha=alpha)

		if (verbose) {
        	cat(paste("For i=", i, ", target:", round(target,3), ", bw=", round(bandwidth,3), ", tolerance=", round(target/1000,3), ", miss=", round(opt$objective,3), ".\n", sep=''))
        }
    }

    gwglmnet.object[['models']] = models
    
	if (tuning) {
    } else if (predict) {
    } else if (simulation) {
    } else {
	    gwglmnet.object[['coords']] = coords
	}

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
