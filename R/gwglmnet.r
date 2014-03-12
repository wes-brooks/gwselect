gwglmnet <- function(formula, data, family, weights=NULL, coords, fit.loc=NULL, tuning=FALSE, predict=FALSE, simulation=FALSE, oracle=NULL, gweight, bw=NULL, varselect.method=c('AIC','BIC','AICc'), verbose=FALSE, longlat, tol.loc=NULL, N=1, bw.type=c('dist','knn','nen'), parallel=FALSE, D=NULL, interact=FALSE, resid.type=c('deviance','pearson')) {
    if (is(data, "Spatial")) {
        if (!missing(coords)) 
            warning("data is Spatial* object, ignoring coords argument")
        coords <- coordinates(data)
        if ((is.null(longlat) || !is.logical(longlat)) && !is.na(is.projected(data)) && 
            !is.projected(data)) {
            longlat <- TRUE
        }
        else longlat <- FALSE
        data <- as(data, "data.frame")
    }
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
    x <- model.matrix(mt, mf)

    #Get the matrices of distances and weights
    if (is.null(D)) {
        n = dim(coords)[1]
        if (longlat) {
            D = as.matrix(earth.dist(coords),n,n)
        } else {
            Xmat = matrix(rep(coords[,1], times=n), n, n)
            Ymat = matrix(rep(coords[,2], times=n), n, n)
            D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
        }
    }

    res = list()
    resid.type = match.arg(resid.type)
    bw.type = match.arg(bw.type)
    varselect.method = match.arg(varselect.method)

    if (bw.type=='dist') {
        weight.matrix = gweight(D, bw)
        if (parallel) {
            res[['model']] = gwglmnet.fit.fixedbwparallel(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, N=N, coords=coords, oracle=oracle, varselect.method=varselect.method, bw=bw, fit.loc=fit.loc, gwr.weights=weight.matrix, verbose=verbose, interact=interact, resid.type=resid.type)
        } else {
            res[['model']] = gwglmnet.fit.fixedbw(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, N=N, coords=coords, oracle=oracle, varselect.method=varselect.method, bw=bw, fit.loc=fit.loc, gwr.weights=weight.matrix, verbose=verbose, interact=interact, resid.type=resid.type)
        }
    } else {        
        bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0
        beta1 = difmin/300
        beta2 = 10*difmin

        if (bw.type=='nen') {
            if (parallel) {
                res[['model']] = gwglmnet.fit.nenparallel(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, coords=coords, oracle=oracle, fit.loc=fit.loc, N=N, D=D, longlat=longlat, varselect.method=varselect.method, verbose=verbose, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol.loc=tol.loc, interact=interact, resid.type=resid.type)
            } else {
                res[['model']] = gwglmnet.fit.nen(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, coords=coords, oracle=oracle, fit.loc=fit.loc, N=N, D=D, longlat=longlat, varselect.method=varselect.method, verbose=verbose, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol.loc=tol.loc, interact=interact, resid.type=resid.type)
            }
        } else if (bw.type=='knn') {
            if (parallel) {
                res[['model']] = gwglmnet.fit.knnparallel(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, coords=coords, oracle=oracle, fit.loc=fit.loc, N=N, D=D, longlat=longlat, varselect.method=varselect.method, verbose=verbose, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol.loc=tol.loc, interact=interact, resid.type=resid.type)
            } else {
                res[['model']] = gwglmnet.fit.knn(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, coords=coords, oracle=oracle, fit.loc=fit.loc, N=N, D=D, longlat=longlat, varselect.method=varselect.method, verbose=verbose, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol.loc=tol.loc, interact=interact, resid.type=resid.type)
            }
        }
    }

    if (!tuning) {
        res[['data']] = data
        res[['response']] = as.character(formula[[2]])
        res[['family']] = family
        res[['weights']] = weights
        res[['coords']] = coords
        res[['fit.locs']] = fit.loc
        res[['longlat']] = longlat
        res[['gweight']] = gweight
        res[['bw']] = bw
        res[['bw.type']] = bw.type
        res[['varselect.method']] = varselect.method
        res[['interact']] = interact
    }
    class(res) = "gwselect"

    res
}
