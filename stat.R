ks.stat <- function(x,distr,...){
  #assign dist pars
  pars <- list(...)
  for(par in names(pars)){
    assign(par, pars[[par]])
  }
  #sort data
  x <- sort(x)
  n <- length(x)
  #calculating empirical cum distr
  #Sx <- cumsum(table(x))/n
  xpu <- seq(1, n)/n
  xpl <- seq(0, n-1)/n  
  #calculating theoretical cum distr
  if(distr == "lognormal"){
    Fx <- plnorm(x,meanlog,sdlog)
  }else if(distr == "weibull"){
    Fx <- pweibull(x,shape,scale)
  }
  #distance between emprical & theoretical
  #SxFx <- abs(Sx-Fx)
  #Sx1Fx <- abs(c(0,Sx[-n])-Fx)
  #KS statistics
  #ks <- max(SxFx,Sx1Fx)
  #ks <- max(abs(Fx-xpl), abs(xpu-Fx))
  ks <- max(pmax(abs(Fx-xpl), abs(xpu-Fx)))
  return(ks)
}

mse.stat <- function(x,distr,...){
  #assign dist pars
  pars <- list(...)
  for(par in names(pars)){
    assign(par, pars[[par]])
  }
  #sort data
  x <- sort(x)
  n <- length(x)
  #calculating empirical cum distr
  ecdf <- seq(1,n)/n
  #calculating theoretical cum distr
  if(distr == "lognormal"){
    tcdf <- plnorm(x,meanlog,sdlog)
  }else if(distr == "weibull"){
    tcdf <- pweibull(x,shape,scale)
  }
  mse <- sum((tcdf-ecdf)^2)/n
  return(mse)
}


loglik <- function(x,distr,theta){
  if(distr == "lognormal"){
    llh <- sum(dlnorm(x,meanlog = theta[1],sdlog = theta[2],log = TRUE))
  }else if(distr == "weibull"){
    llh <- sum(dweibull(x,shape = theta[1],scale = theta[2],log = TRUE))
  }
  return(llh)
}