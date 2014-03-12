gwglmnet.fit.inner = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, event=NULL, family, mode.select, tuning, predict, simulation, verbose, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, interact, N=1, alpha, shrunk.fit) {
  #Find which observations were made at the model location  
  if (!is.null(indx)) {colocated = which(round(coords[indx,1],5)==round(as.numeric(loc[1]),5) & round(coords[indx,2],5) == round(as.numeric(loc[2]),5))}
  else {colocated = which(round(coords[,1],5) == round(as.numeric(loc[1]),5) & round(coords[,2],5) == round(as.numeric(loc[2]),5))}

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

  #Establish groups for the group lasso
  vargroup = 1:ncol(x)
  
  #Compute the covariate-by-location interactions
  raw.names = colnames(x)
  if (interact) {
    interact.names = vector()
    for (l in 1:length(raw.names)) {
      interact.names = c(interact.names, paste(raw.names[l], ":", colnames(coords)[1], sep=""))
      interact.names = c(interact.names, paste(raw.names[l], ":", colnames(coords)[2], sep=""))
    }
    interacted = matrix(ncol=2*ncol(x), nrow=nrow(x))
    for (k in 1:ncol(x)) {
      interacted[,2*(k-1)+1] = x[,k]*(coords[,1]-loc[1])
      interacted[,2*k] = x[,k]*(coords[,2]-loc[2])
      vargroup = c(vargroup, vargroup[k], vargroup[k])
    }
    x.interacted = cbind(x, interacted)
    colnames(x.interacted) = c(raw.names, interact.names)
  }

  #Combine prior weights and kernel weights
  w <- prior.weights * gwr.weights
  if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 
  weighted = which(w>0)
  n.weighted = length(weighted)
  
  #Limit our attention to the observations with nonzero weight
  xx = as.matrix(x[weighted,])
  if (interact) {xx.interacted = as.matrix(x.interacted[weighted,])}
  yy = as.matrix(y[weighted])
  w = w[weighted]

  #Instantiate objects to store our output
  int.list = list()
  coef.list = list()
  coef.unshrunk.list=list()    
  coef.unshrunk.interacted.list=list()
  tunelist = list()

  #Loop through N bootstrap replicates
  for (i in 1:N) {
      
    if (i==N) {
      #Final permutation is the original ordering of the data:
      permutation = 1:n.weighted
    } else {
      #Generate permutations to preserve the total local weight:
      tot.w = sum(w)
      permutation = vector()
      while (sum(w[permutation]) < tot.w) {
        permutation = c(permutation, sample(1:n.weighted, size=1))
      }

      permutation = permutation[1:which.min(tot.w - cumsum(w[permutation]))]
    }

    #Sample the data according to the permutation:
    colocated = which(gwr.weights[weighted][permutation]==1)
    xxx = xx[permutation,]
    if (interact) {xxx = xx.interacted[permutation,]}
    yyy = yy[permutation]
    sumw = sum(w[permutation])
    
    #Use the adaptive group lasso to produce a local model:
    model = SGL(data=list(x=xxx, y=yyy), weights=w[permutation], index=vargroup, standardize=FALSE, alpha=0, nlam=100, min.frac=0.0001, adaptive=TRUE)

    
    vars = apply(as.matrix(model[['beta']]), 2, function(x) {which(x!=0)})
    df = sapply(vars, length)

    if (sumw > ncol(x)) {
      #Extract the fitted values for each lambda:
      fitted = model[['results']][['fitted']] #predict(model, newx=predx, type="response")
      s2 = sum(w[permutation]*(model[['results']][['residuals']][,ncol(fitted)])**2) / (sumw - df) 
      
      #Compute the loss (varies by family)
      #loss = model[[mode.select]]
      if (mode.select == 'AIC') {penalty = 2*df}
      if (mode.select == 'AICc') {penalty = 2*df + 2*df*(df+1)/(sumw - df - 1)}
      if (mode.select == 'BIC') {penalty = sumw*df}


#Assuming scale from the largest model:
#loss = sumw * log(s2) + apply(model[['results']][['residuals']], 2, function(x) sum(w[permutation]*x**2))/s2 + penalty

#Estimating scale in penalty formula:
loss = sumw * (log(apply(model[['results']][['residuals']], 2, function(x) sum(w[permutation]*x**2))) - log(sumw)) + penalty

#Estimating the loss only at the modeling location (not the total local loss)
#loss = (model[['results']][['residuals']][colocated,])**2 + penalty*w[permutation][colocated] / sumw

      #Pick the lambda that minimizes the loss:
      k = which.min(loss)
      fitted = fitted[,k]
      localfit = fitted[colocated]
      df = df[k]
      if (k > 1) {
        varset = vars[[k]]
      } else {
        varset = NULL
      }
      
      if (length(colocated)>0) {
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

        #Compute the dispersion parameter:
        if (family=='gaussian') { tunelist[['s2']] = s2 }
        else if (family=='poisson') { tunelist[['s2']] = summary(m)$dispersion }
        else if (family=='binomial') { tunelist[['s2']] = 1 }
        
        #Prepare some outputs for the bandwidth-finding scheme:
        tunelist[['n']] = sumw
        tunelist[['df']] = df
        tunelist[['df-local']] = df*w[permutation][colocated] / sumw
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
    
    #Get the coefficients:
    coefs = t(rbind(model[['intercept']], model[['beta']])[k,]
    coefs = Matrix(coefs, ncol=1)
    rownames(coefs) = c("(Intercept)", colnames(xxx))

    coef.list[[i]] = coefs
  }
  
  if (tuning) {
    return(list(tunelist=tunelist, s=k, sigma2=s2, nonzero=colnames(xxx)[vars[[k]]], weightsum=sum(w), loss=loss, alpha=alpha))
  } else if (predict) {
    return(list(tunelist=tunelist, coef=coefs, weightsum=sum(w), s=k, sigma2=s2, nonzero=colnames(xxx)[vars[[k]]]))
  } else if (simulation) {
    #return(list(tunelist=tunelist, coef=coefs, coeflist=coef.list, s=k, bw=bw, sigma2=s2, coef.unshrunk=coefs.unshrunk, s2.unshrunk=s2.unshrunk, coef.unshrunk.list=coef.unshrunk.list, fitted=localfit, alpha=alpha, nonzero=colnames(x)[vars[[k]]], actual=predy[colocated], weightsum=sum(w), loss=loss))
    return(list(tunelist=tunelist, coef=coefs, coeflist=coef.list, s=k, bw=bw, sigma2=s2, fitted=localfit, alpha=alpha, nonzero=colnames(xxx)[vars[[k]]], actual=yyy[colocated], weightsum=sum(w), loss=loss))
  } else {
    return(list(model=model, loss=loss, coef=coefs, coeflist=coef.list, nonzero=colnames(xxx)[vars[[k]]], s=k, loc=loc, bw=bw, df=df, loss.local=loss, sigma2=s2, sum.weights=sum(w), N=N, fitted=localfit, alpha=alpha, weightsum=sum(w)))
  }
}
