gwglmnet.cv.f = function(formula, data, weights, indx, family, bw, coords, gweight, env, oracle, mode.select, verbose, adapt, longlat, s, tol, method, N, parallel, precondition, interact, alpha, shrunk.fit, bw.select resid.type) {    
    #Generate the model with the given bandwidth:
    cat(paste("starting bw:", round(bw, 3), '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, family=family,
        weights=weights, tuning=TRUE, indx=indx, coords=coords, gweight=gweight,
        oracle=oracle, bw=bw, N=N, mode.select=mode.select, verbose=verbose,
        longlat=longlat, adapt=adapt, s=s, method=method, parallel=parallel,
        precondition=precondition, interact=interact, alpha=alpha,
        shrunk.fit=shrunk.fit, bw.select=bw.select, resid.type=resid.type)

    if (bw.select=='AICc') {
        trH = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {tail(x[['tunelist']][['trace.local']],1)})) 
       	loss = nrow(data) * (log(mean(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['ssr.local']]}))) + 1 + (2*(trH+1))/(nrow(data)-trH-2) + log(2*pi))
    } else if (bw.select=='GCV') {
        trH = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {tail(x[['tunelist']][['trace.local']],1)})) 
        if (resid.type=='deviance') {
            loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['tunelist']]})) / (nrow(data)-trH)**2
        } else if (resid.type=='pearson') {
            loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['tunelist']]})) / (nrow(data)-trH)**2
        }
    } else if (bw.select=='BICg') {
        loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {min(x[['loss.local']])}))
        #"Simplistic" BIC - based on eq4.22 from the Fotheringham et al. book:
        #loss = nrow(data) * (log(mean(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['ssr.local']]}))) + 1 + log(2*pi)) + trH * log(nrow(data))/2
    }

    res = mget('trace', env=env, ifnotfound=list(matrix(NA, nrow=0, ncol=2)))
    res$trace = rbind(res$trace, c(bw, loss))
    assign('trace', res$trace, env=env)

    cat(paste('Bandwidth: ', round(bw, 3), '. Loss: ', signif(loss, 5), '\n', sep=''))
    return(loss)
}
