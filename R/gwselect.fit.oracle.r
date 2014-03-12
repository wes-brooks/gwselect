gwselect.fit.oracle = function(x, y, coords, loc, family='gaussian', oracle, tuning, predict, simulation, verbose, interact, gwr.weights, prior.weights, longlat, N) {

  colocated = which(round(coords[,1],5) == round(as.numeric(loc[1]),5) & round(coords[,2],5) == round(as.numeric(loc[2]),5))        

  #Establish the oracular data frame (possibly with location interactions)
  xx = as.matrix(x[,oracle])   
  colnames(xx) = c(oracle)
  if (interact && ncol(xx)>0) {
    newnames = vector()
    for (l in 1:length(oracle)) {
      newnames = c(newnames, paste(oracle[l], ":", colnames(coords)[1], sep=""))
      newnames = c(newnames, paste(oracle[l], ":", colnames(coords)[2], sep=""))
    }

    interacted = matrix(ncol=2*ncol(xx), nrow=nrow(xx))
    for (k in 1:ncol(xx)) {
      interacted[,2*(k-1)+1] = xx[,k] * (coords[,1]-loc[1,1])
      interacted[,2*k] = xx[,k] * (coords[,2]-loc[1,2])
    }
    xx = cbind(xx, interacted)
    colnames(xx) = c(oracle, newnames)
  }

  ssr.local = NA
  df = length(oracle) + 1
  
  yy = as.matrix(y)    
  w <- prior.weights * gwr.weights
  
  if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 

  n <- nrow(xx)
  weighted = which(w>0)
  n.weighted = length(weighted)
  
  xx = xx[weighted,]
  yy = as.matrix(yy[weighted])
  w = w[weighted]
  colocated = which(gwr.weights[weighted]==1)
  fitdata = data.frame(y=yy, xx)
  localdata = data.frame(fitdata[colocated,])
  colnames(fitdata) = colnames(localdata) = c("y", colnames(xx))
  
  int.list = list()
  coef.list = list()
  tunelist = list()
  
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
    
    permuted = data.frame(fitdata[permutation,])
    colnames(permuted) = colnames(fitdata)

    model = glm(y~., data=permuted, weights=w[permutation], family=family)
    colocated = which(gwr.weights[weighted][permutation]==1)
    yyy = yy[permutation]
            
    #Get the coefficients:
    coefs = rep(0, ncol(x)+1)
    names(coefs) = c("(Intercept)", colnames(x))
    coefs[c("(Intercept)", oracle)] = coef(model)[c("(Intercept)", oracle)]
    coefs = Matrix(coefs, ncol=1)
    rownames(coefs) = c("(Intercept)", colnames(x))
    coef.list[[i]] = coefs
  
  	if (i==N) { 
    	if (sum(w) > dim(permuted)[2]) {  		
  			s2 = sum(w[permutation] * model$residuals**2)/(sum(w) - df)
        
				fitted = predict(model, newdata=permuted, type='response')
		    Xh = diag(sqrt(w[permutation])) %*% as.matrix(cbind(rep(1,length(permutation)), xx))
        H = Xh %*% solve(t(Xh) %*% Xh) %*% t(Xh)
        Hii = sum(H[colocated,colocated])

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

          if (family=='gaussian') { tunelist[['s2']] = s2 }
          else if (family=='poisson') { tunelist[['s2']] = summary(model)$dispersion }
          else if (family=='binomial') { tunelist[['s2']] = 1 }
          
          tunelist[['n']] = sum(w[permutation])
          tunelist[['trace.local']] = Hii
          tunelist[['df']] = df
          tunelist[['df-local']] = df*w[permutation][colocated] / sum(w[permutation])
				} else {
					loss.local = NA
				}
		  } else {
  			loss.local = Inf
  			fitted = NA
  			s2 = NA
  			coef.list = NA
  			coefs = NA
		  }
    }
  }
  
  #Return the results
  if (tuning) {
    return(list(tunelist=tunelist, s=NULL, sigma2=s2, nonzero=oracle, weightsum=sum(w)))
  } else if (predict) {
    return(list(tunelist=tunelist, coef=coefs))
  } else if (simulation) {
    return(list(tunelist=tunelist, coef=coefs, coeflist=coef.list, bw=bw, sigma2=s2, fitted=fitted[colocated][1], nonzero=oracle, weightsum=sum(w), s=NULL))
  } else {
    return(list(model=model, coef=coefs, coeflist=coef.list, loc=loc, bw=bw, tunelist=tunelist, sigma2=s2, nonzero=oracle, sum.weights=sum(w), N=N))
  }
}
