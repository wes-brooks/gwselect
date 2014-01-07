gwglmnet.fit.knn = function(x, y, family, coords, fit.loc, oracle, D, verbose, mode.select, prior.weights, tuning, predict, simulation, indx, gweight, target, beta1, beta2, tol.loc, longlat=FALSE, precondition=FALSE, N, interact, alpha, shrunk.fit) {
    if (!is.null(fit.loc)) {
        coords.unique = unique(fit.loc)
    } else {
        coords.unique = unique(coords)
    }
    n = dim(coords.unique)[1]

    gwglmnet.object = list()
    models = list()

    prior.weights = drop(prior.weights)
    max.weights = rep(1, length(prior.weights))
    total.weight = sum(max.weights * prior.weights)
    if (is.null(tol.loc)) {tol.loc = target / 1000}

    for (i in 1:n) {
        loc = coords.unique[i,]
        dist = drop(D[i,])

        opt = optimize(gwglmnet.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=tol.loc, coords=coords, loc=loc, indx=indx,
            gweight=gweight, verbose=verbose, dist=dist, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum

		if (is.null(oracle)) {
	        models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, indx=indx, bw=bandwidth, dist=dist, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, precondition=precondition, N=N, interact=interact, alpha=alpha, shrunk.fit=shrunk.fit)
        } else {
            models[[i]] = gwselect.fit.oracle(x=x, y=y, bw=bandwidth, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, dist=dist, prior.weights=prior.weights, gweight=gweight, interact=interact)
        }
        
        if (verbose) {
	        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bandwidth, 3), "; loss=", round(models[[i]][['tunelist']][['ssr-loc']][[resid.type]],3), "; s=", models[[i]][['s']], "; sigma2=", round(tail(models[[i]][['sigma2']],1),3), "; nonzero=", paste(models[[i]][['nonzero']], collapse=","), "; weightsum=", round(models[[i]][['weightsum']],3), ".\n", sep=''))
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
