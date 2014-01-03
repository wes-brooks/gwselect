gwglmnet.fit.nen = function(x, y, family, coords, indx, fit.loc, D, s, mode.select, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE, tuning, simulation, predict, N, oracle, interact, alpha, shrunk.fit, bw.select, resid.type) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwglmnet.object = list()
    models = list()

    if (verbose) {cat(paste('beta1:', round(beta1,3), ', beta2:', round(beta2,3), '\n', sep=''))}

    for (i in 1:n) {
        loc = coords.unique[i,]
        dist = D[i,]

        opt = optimize(gwglmnet.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, coords=coords, loc=loc, s=s, alpha=alpha, bw.select=bw.select, resid.type=resid.type,
            gweight=gweight, verbose=verbose, dist=dist, mode.select=mode.select, adapt=adapt, family=family, oracle=oracle,
            prior.weights=prior.weights, target=target, precondition=precondition, interact=interact, shrunk.fit=shrunk.fit)
        bandwidth = opt$minimum

        models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bandwidth, dist=dist, s=s, mode.select=mode.select, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, tuning=tuning, simulation=simulation, predict=predict, precondition=precondition, N=N, interact=interact, shrunk.fit=shrunk.fit, alpha=alpha, bw.select=bw.select, resid.type=resid.type)

		if (verbose) {
        	cat(paste("For i=", i, ", target:", round(target,3), ", bw=", round(bandwidth,3), ", tolerance=", round(target/1000,3), ", miss=", round(opt$objective,3), ".\n", sep=''))
        }
    }

    gwglmnet.object[['models']] = models
    gwglmnet.object[['coords']] = coords
    gwglmnet.object[['s.range']] = s

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
