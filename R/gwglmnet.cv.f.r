gwglmnet.cv.f = function(formula, data, weights, indx, family, bw, coords, gweight, env, oracle, mode.select, verbose, longlat, tol.loc, bw.method, N, parallel, precondition, interact, alpha, shrunk.fit, bw.select, resid.type) {    
    #Generate the model with the given bandwidth:
    cat(paste("starting bw:", round(bw, 3), '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, family=family,
        weights=weights, tuning=TRUE, indx=indx, coords=coords, gweight=gweight,
        oracle=oracle, bw=bw, N=N, mode.select=mode.select, verbose=verbose,
        longlat=longlat, bw.method=bw.method, parallel=parallel,
        precondition=precondition, interact=interact, alpha=alpha, tol.loc=tol.loc,
        shrunk.fit=shrunk.fit, resid.type=resid.type)
    
    if (bw.select=='AICc') {
        trH = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {tail(x[['tunelist']][['trace.local']],1)}))
       	loss = nrow(data) * (log(mean(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['tunelist']][['ssr-loc']][[resid.type]]}))) + 1 + (2*(trH+1))/(nrow(data)-trH-2) + log(2*pi))
    } else if (bw.select=='GCV') {
        trH = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {tail(x[['tunelist']][['trace.local']],1)})) 
        loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['tunelist']][['ssr-loc']][[resid.type]]})) / (nrow(data)-trH)**2
    } else if (bw.select=='BICg') {
        loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {
            s2 = x[['tunelist']][['s2']]
            if (family=='gaussian') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]])/s2 + log(s2) }
            else if (family=='binomial') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]]) }
            else if (family=='poisson') { ll = min(x[['tunelist']][['ssr-loc']][[resid.type]])/s2 }
            df = v[['tunelist']][['df']]
            return(ll + log(x[['tunelist']][['n']]) * df / x[['tunelist']][['n']])
            }))
        loss = loss + sum(sapply(gwglmnet.model[['model']][['models']], function(x) {min(x[['tunelist']][['ssr-loc']][[resid.type]])}))
        #"Simplistic" BIC - based on eq4.22 from the Fotheringham et al. book:
        #loss = nrow(data) * (log(mean(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['ssr.local']]}))) + 1 + log(2*pi)) + trH * log(nrow(data))/2
    }

    res = mget('trace', env=env, ifnotfound=list(matrix(NA, nrow=0, ncol=2)))
    res$trace = rbind(res$trace, c(bw, loss))
    assign('trace', res$trace, env=env)

    cat(paste('Bandwidth: ', round(bw, 3), '. Loss: ', signif(loss, 5), '\n', sep=''))
    return(loss)
}
