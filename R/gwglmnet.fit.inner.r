gwglmnet.fit.inner = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, event=NULL, family, mode.select, tuning, predict, simulation, verbose, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, interact, N=1, alpha, tau=3, shrunk.fit) {
    if (!is.null(indx)) {
        colocated = which(round(coords[indx,1],5)==round(as.numeric(loc[1]),5) & round(coords[indx,2],5) == round(as.numeric(loc[2]),5))
    }
    else {
        colocated = which(round(coords[,1],5) == round(as.numeric(loc[1]),5) & round(coords[,2],5) == round(as.numeric(loc[2]),5))        
    }
    reps = length(colocated)

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)     
    } 
    gwr.weights = drop(gwr.weights)  

    if (!is.null(indx)) {
        gwr.weights = gwr.weights[indx]
    }

	  #Allow for the adaptive elastic net penalty:
	  if (substring(as.character(alpha), 1, 1) == 'a') {
		    cormat = abs(cor(x))
		    diag(cormat) = NA
        alpha = 1 - max(cormat, na.rm=TRUE)
    }

    #For interaction on location:
    groups = 0:(ncol(x)-1)
    oldnames = colnames(x)
    if (interact) {
        newnames = vector()
        for (l in 1:length(oldnames)) {
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[1], sep=""))
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[2], sep=""))
        }
        interacted = matrix(ncol=2*ncol(x), nrow=nrow(x))
        for (k in 1:ncol(x)) {
            interacted[,2*(k-1)+1] = x[,k]*(coords[,1]-loc[1])
            interacted[,2*k] = x[,k]*(coords[,2]-loc[2])
            groups = c(groups, groups[k], groups[k])
        }
        x.interacted = cbind(x, interacted)
        colnames(x.interacted) = c(oldnames, newnames)
    }

    if (mode.select=='CV') { 
        xx = as.matrix(x[-colocated,])
        if (interact) {xx.interacted = as.matrix(x.interacted[-colocated,])}
        yy = as.matrix(y[-colocated])
        w <- prior.weights[-colocated] * gwr.weights[-colocated]  
    } else {
        xx = as.matrix(x)
        if (interact) {xx.interacted = as.matrix(x.interacted)}
        yy = as.matrix(y)
        w <- prior.weights * gwr.weights
    }
    
    if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 

    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
    
    xx = xx[weighted,]
    if (interact) {xx.interacted = xx.interacted[weighted,]}
    yy = as.matrix(yy[weighted])
    w = w[weighted]

    int.list = list()
    coef.list = list()
    coef.unshrunk.list=list()    
    coef.unshrunk.interacted.list=list()
    tunelist = list()

    for (i in 1:N) {
        #Final permutation is the original ordering of the data:
        if (i==N) {            
            permutation = 1:n.weighted
        } else {
            tot.w = sum(w)
            permutation = vector()
            while (sum(w[permutation]) < tot.w) {
                permutation = c(permutation, sample(1:n.weighted, size=1))
            }

            permutation = permutation[1:which.min(tot.w - cumsum(w[permutation]))]
        }

        colocated = which(gwr.weights[weighted][permutation]==1)
        xxx = xx[permutation,]
        if (interact) {xxx = xx.interacted[permutation,]}
        yyy = yy[permutation]
        
        model = SGL(data=list(x=xxx, y=yyy), weights=w[permutation], index=groups, standardize=FALSE, alpha=0, nlam=100, min.frac=0.0001, adaptive=TRUE)
print(model[['beta']])
        nsteps = length(model$lambda) + 1       
        vars = apply(as.matrix(model[['beta']]), 2, function(x) {which(x!=0)})
print(vars)
        df = sapply(vars, length)
print(df)

        if (sum(w) > ncol(x)) {
            #Extract the fitted values for each lambda:
            coefs = t(as.matrix(model[['beta']]))
            fitted = model[['results']][['fitted']] #predict(model, newx=predx, type="response")  

            s2 = sum(w[permutation]*(model[['results']][['residuals']][,ncol(fitted)])**2) / (sum(w[permutation]) - df) 
            
            #Compute the loss (varies by family)
            #loss = model[[mode.select]]
            if (mode.select == 'AIC') {penalty = 2*df}
            if (mode.select == 'AICc') {penalty = 2*df + 2*df*(df+1)/(sum(w[permutation]) - df - 1)}
            if (mode.select == 'BIC') {penalty = sum(w[permutation])*df}
print(apply(model[['results']][['residuals']], 2, function(x) sum(w[permutation] * x**2)))
#print(model[['residuals']]**2)
#print(w[permutation])
print(penalty)
            #loss = (model[['results']][['residuals']][colocated,])**2 + penalty*w[permutation][colocated] / sum(w[permutation])
loss = apply(model[['results']][['residuals']], 2, function(x) sum(w[permutation]*x**2)) + penalty

            #Pick the lambda that minimizes the loss:
            k = which.min(loss)
print(loss)
print(k)
            fitted = fitted[,k]
            localfit = fitted[colocated]
            df = df[k]
            if (k > 1) {
                varset = vars[[k]]
                
                #modeldata = data.frame(y=yy[permutation], xxx[,varset])
                #m = glm(y~., data=modeldata, weights=w[permutation], family=family)
                #working.weights = as.vector(m$weights)
                #result = tryCatch({
                #    Xh = diag(sqrt(working.weights)) %*% as.matrix(cbind(rep(1,length(permutation)), xxx[,varset]))
                #    H = Xh %*% solve(t(Xh) %*% Xh) %*% t(Xh)
                #    Hii = H[colocated,colocated]
                #}, error = function(e) {
                #    Hii = nrow(x) - 2
                #})
                #if (!shrunk.fit) {
                #    fitted = m$fitted
                #    localfit = fitted[colocated]
                #    df = length(varset) + 1
                #    s2 = sum((w[permutation]*m$residuals**2) / (sum(w[permutation]) - df))
                #}
                #coefs.unshrunk = rep(0, ncol(x) + 1)
                #coefs.unshrunk[c(1, varset + 1)] = coef(m)
                #s2.unshrunk = sum(w[permutation]*m$residuals**2)/sum(w[permutation])
            } else {
                #modeldata = data.frame(y=yy[permutation], xxx)
                #m = glm(y~1, data=modeldata, weights=w[permutation], family=family)
                
                #coefs.unshrunk = rep(0, ncol(xx) + 1)
                #coefs.unshrunk[1] = sum(fity * w[permutation]) / sum(w[permutation])
                #s2.unshrunk = sum((sqrt(w[permutation])*fity)**2)/sum(w[permutation])

                #Hii = 1 / sum(w[permutation])
            }
            
            if (length(colocated)>0) {
cat(paste("colocated obs: ", length(colocated), "\n", sep=""))
cat(paste("family: ", family, "\n", sep=""))
                tunelist[['ssr-loc']] = list()
                tunelist[['ssr']] = list()
                
                #Pearson residuals:
                if (family=='gaussian') {
                    tunelist[['ssr-loc']][['pearson']] = sum((w[permutation]*(fitted - yyy)**2)[colocated])
                    tunelist[['ssr']][['pearson']] = sum(w[permutation]*(fitted - yyy)**2)
                } else if (family=='poisson') {
                    tunelist[['ssr-loc']][['pearson']] = sum((w[permutation]*(yyy - fitted)**2/fitted)[colocated])
                    tunelist[['ssr']][['pearson']] = sum(w[permutation]*(fitted - yyy)**2/fitted)
                } else if (family=='binomial') {
                    tunelist[['ssr-loc']][['pearson']] = sum((w[permutation]*(yyy - fitted)**2/(fitted*(1-fitted)))[colocated])
                    tunelist[['ssr']][['pearson']] = sum(w[permutation]*(fitted - yyy)**2/(fitted*(1-fitted)))
                }

                #Deviance residuals:
                if (family=='gaussian') {
                    tunelist[['ssr-loc']][['deviance']] = sum((w[permutation]*(fitted - yyy)**2)[colocated])
                    tunelist[['ssr']][['deviance']] = sum(w[permutation]*(fitted - yyy)**2)
                } else if (family=='poisson') {
                    tunelist[['ssr-loc']][['deviance']] = sum((2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - (yyy-fitted)))[colocated])
                    tunelist[['ssr']][['deviance']] = sum(2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - (yyy-fitted)))
                } else if (family=='binomial') {
                    tunelist[['ssr-loc']][['deviance']] = sum((2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - ylogy(1-yyy) + (1-yyy)*log(1-fitted)))[colocated])
                    tunelist[['ssr']][['deviance']] = sum(2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - ylogy(1-yyy) + (1-yyy)*log(1-fitted)))
                }

                if (family=='gaussian') {
                    tunelist[['s2']] = s2
                } else if (family=='poisson') { 
                    tunelist[['s2']] = summary(m)$dispersion
                } else if (family=='binomial') {
                    tunelist[['s2']] = 1
                }
                tunelist[['n']] = sum(w[permutation])
                #tunelist[['trace.local']] = Hii
                tunelist[['df']] = df
                tunelist[['df-local']] = df*w[permutation][colocated] / sum(w[permutation])
            } else {
                loss.local = NA
            }                   
        } else {
            fitted = rep(meany, length(permutation))
            s2 = 0
            loss = Inf
            loss.local = c(Inf)   
            localfit = meany
        }
    
        #Get the tuning parameter to minimize the loss:
        s.optimal = which.min(loss)
        
        #We have all we need for the tuning stage.
        if (!tuning) {
            #Get the coefficients:
            coefs = coefs[s.optimal,]
            coefs = Matrix(coefs, ncol=1)
            #rownames(coefs) = c("(Intercept)", colnames(x))   

            #coefs = coefs * c(1, adapt.weight) * c(1, 1/normx)
            #if (length(coefs)>1) {coefs[1] = mean(sqrt(w[permutation])*fity) - sum(coefs[2:length(coefs)] * drop(sqrt(w[permutation]) %*% xxx) / nrow(xxx))}
    
            #coefs.unshrunk = Matrix(coefs.unshrunk[1:(ncol(x)+1)], ncol=1)
            #rownames(coefs.unshrunk) = c("(Intercept)", oldnames)  
    
            #coef.unshrunk.list[[i]] = coefs.unshrunk
            #coef.list[[i]] = coefs
        }
    }
    
    if (tuning) {
        return(list(tunelist=tunelist, s=s.optimal, sigma2=s2, nonzero=colnames(xxx)[vars[[s.optimal]]], weightsum=sum(w), loss=loss, alpha=alpha))
    } else if (predict) {
        return(list(tunelist=tunelist, coef=coefs, weightsum=sum(w), s=s.optimal, sigma2=s2, nonzero=colnames(xxx)[vars[[s.optimal]]]))
    } else if (simulation) {
        #return(list(tunelist=tunelist, coef=coefs, coeflist=coef.list, s=s.optimal, bw=bw, sigma2=s2, coef.unshrunk=coefs.unshrunk, s2.unshrunk=s2.unshrunk, coef.unshrunk.list=coef.unshrunk.list, fitted=localfit, alpha=alpha, nonzero=colnames(x)[vars[[s.optimal]]], actual=predy[colocated], weightsum=sum(w), loss=loss))
      return(list(tunelist=tunelist, coef=coefs, coeflist=coef.list, s=s.optimal, bw=bw, sigma2=s2, fitted=localfit, alpha=alpha, nonzero=colnames(xxx)[vars[[s.optimal]]], actual=yyy[colocated], weightsum=sum(w), loss=loss))
    } else {
        return(list(model=model, loss=loss, coef=coefs, coeflist=coef.list, nonzero=colnames(xxx)[vars[[s.optimal]]], s=s.optimal, loc=loc, bw=bw, df=df, loss.local=loss, sigma2=s2, sum.weights=sum(w), N=N, fitted=localfit, alpha=alpha, weightsum=sum(w)))
    }
}
