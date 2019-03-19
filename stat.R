#Goodness of Fit Measures Used
#KS, MSE of distirbution, Log-likelihood
#
#Ari Purwanto Sarwo Prasojo & Puguh Prasetyoputra (2019)
#Research Center for Population, Indonesian Institute of Sciences
#________________________________________________________________


#Calculate KS distance ----
ks.stat <- function(x,distr,...){
  
#________________________________________________________________
#ks.stat
#Calculating ks distance for lognormal and weibull distribution
#
#/Usage
#ks.stat(x,distr,...)
#/Arguments
#x              : univariate continuous random sample
#distr          : choices of distribution ("lognormal" or "weibull")
#...          	: additional arguments (parameters of distribution)
#
#/Value
#ks            	: ks distance
#________________________________________________________________
  
  #assign distr pars
  pars <- list(...)
  for(par in names(pars)){
    assign(par, pars[[par]])
  }
  #sort data
  x <- sort(x)
  n <- length(x)
  #calculating empirical c.d.f
  ecdfx <- ecdf(x)
  Sx <- ecdfx(x)
  Sx1 <- c(0, Sx[-n])
  #calculating theoretical c.d.f
  if(distr == "lognormal"){
    Fx <- plnorm(x,meanlog,sdlog)
  }else if(distr == "weibull"){
    Fx <- pweibull(x,shape,scale)
  }
  #max distance between emprical & theoretical c.d.f
  ks <- max(pmax(abs(Fx-Sx1), abs(Sx-Fx)))
  return(ks)
}


#Calculate MSE of distribution fitting ----
mse.stat <- function(x,distr,...){
  
#________________________________________________________________
#mse.stat
#Calculating mean square error (mse) of distribution fitting
#
#/Usage
#mse.stat(x,distr,...)
#/Arguments
#x              : univariate continuous random sample
#distr          : choices of distribution ("lognormal" or "weibull")
#...          	: additional arguments (parameters of distribution)
#
#/Value
#mse          	: mse of distribution fitting
#________________________________________________________________
  
  #assign distr pars
  pars <- list(...)
  for(par in names(pars)){
    assign(par, pars[[par]])
  }
  #sort data
  x <- sort(x)
  n <- length(x)
  #calculating empirical c.d.f
  ecdfx <- ecdf(x)
  Sx <- ecdfx(x)
  #calculating theoretical c.d.f
  if(distr == "lognormal"){
    Fx <- plnorm(x,meanlog,sdlog)
  }else if(distr == "weibull"){
    Fx <- pweibull(x,shape,scale)
  }
  mse <- sum((Fx-Sx)^2)/n
  return(mse)
}


#Calculate Log-likelihood ----
loglik <- function(x,distr,theta){
  
#________________________________________________________________
#loglik
#Calculating log-likelihood of lognormal or weibull distribution
#
#/Usage
#loglik(x,distr,theta)
#/Arguments
#x              : univariate continuous random sample
#distr          : choices of distribution ("lognormal" or "weibull")
#theta        	: parameters of distirbution (meanlog,sdlog) or (shape,scale)
#
#/Value
#llh          	: log-likelihood value
#________________________________________________________________  
  
  if(distr == "lognormal"){
    llh <- sum(dlnorm(x,meanlog = theta[1],sdlog = theta[2],log = TRUE))
  }else if(distr == "weibull"){
    llh <- sum(dweibull(x,shape = theta[1],scale = theta[2],log = TRUE))
  }
  return(llh)
}