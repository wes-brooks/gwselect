gwglmnet.sel = function(formula, data=list(), family, range=NULL, weights=NULL, coords, oracle=NULL, indx=NULL, gweight=gwr.Gauss, bw.method="dist", mode.select='AIC', verbose=FALSE, longlat=FALSE, tol.loc=.Machine$double.eps^0.25, tol.bw=.Machine$double.eps^0.25, parallel=FALSE, alpha=1, precondition=FALSE, interact=FALSE, shrunk.fit=TRUE, bw.select=c('AICc', 'GCV', 'BICg'), resid.type=c('deviance', 'pearson')) {

    if (is.null(longlat) || !is.logical(longlat)) 
        longlat <- FALSE
    if (missing(coords)) 
        stop("Observation coordinates have to be given")
        
    mf <- match.call(expand.dots = FALSE)
    #m <- match(c("formula", "data", "weights"), names(mf), 0)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    #weights <- as.vector(model.extract(mf, "weights"))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) 
        weights <- rep(as.numeric(1), dp.n)
    if (any(is.na(weights))) 
        stop("NAs in weights")
    if (any(weights < 0)) 
        stop("negative weights")
    y <- model.extract(mf, "response")

        
        
    if (!is.null(range)) {
        beta1 = min(range)
        beta2 = max(range)
    } else {
        if (bw.method == "dist") {
            bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
            difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
            if (any(!is.finite(difmin))) 
                difmin[which(!is.finite(difmin))] <- 0
            beta1 <- difmin/1000
            beta2 <- 10*difmin
        } else if (bw.method == 'knn') {
            beta1 <- 0
            beta2 <- 1
        } else if (bw.method == 'nen') {
            if (family=='binomial') {beta2 = sum(weights/(mean(y)*(1-mean(y))) * (y-mean(y))**2)}
            else if (family=='poisson') {beta2 = sum(weights/(mean(y)) * (y-mean(y))**2)}
            else if (family=='gaussian') {beta2 = sum(weights * (y-mean(y))**2)}
            beta1 = beta2/1000
        }
    }

    #Create a new environment, in which we will store the likelihood trace from bandwidth selection.
    oo = new.env()
    opt <- optimize(gwglmnet.cv.f, interval=c(beta1, beta2), tol=tol.bw, maximum=FALSE,
        formula=formula, indx=indx, coords=coords, env=oo, oracle=oracle, family=family, mode.select=mode.select,
        gweight=gweight, verbose=verbose, longlat=longlat, data=data, bw.method=bw.method, alpha=alpha, shrunk.fit=shrunk.fit,
        weights=weights, tol.loc=tol.loc, parallel=parallel, precondition=precondition, N=1, interact=interact,
        resid.type=resid.type, bw.select=bw.select)
    trace = oo$trace
    rm(oo)

    bdwt <- opt$minimum
    res <- bdwt
    return(list(bw=res, trace=trace, bw.select=bw.select, resid.type=resid.type))
}
