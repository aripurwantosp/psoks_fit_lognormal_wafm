#Fit Distribution using minimum KS or MSE fit
#For lognormal and weibull distribution
#
#Ari Purwanto Sarwo Prasojo & Puguh Prasetyoputra (2019)
#Research Center for Population, Indonesian Institute of Sciences
#________________________________________________________________


pso.fitdist <- function(dt,distr,stat,limit,cf=c(2,2),max.iter=200,n.swarm=20){
  
#________________________________________________________________
#pso.fitdist
#Lognormal, Weibull parameter estimation through minimizing KS or MSE using PSO
#Require stat.r source
#
#/Usage
#pso.fitdist(dt,distr,stat,limit,cf,max.iter,n.swarm)
#/Arguments
#dt             : univariate continuous random sample
#distr          : choices of distribution ("lognormal" or "weibull")
#stat			      : stat function to minimize
#limit        	: searching limit c(lower, upper)
#cf				      : social factor (default (c(2,2))
#max.iter     	: maximum iteration (default 100)
#n.swarm      	: particle size
#
#/Value
#solution     	: estimators
#stat           : string indicate that stistics which minimize
#stat.value    	: statistic value, according stat
#loglik    	    : log-likelihood value
#iter         	: iteration
#stat.trace    	: behaviour's record
#________________________________________________________________

  #/fitness function
  fn.fit <- function(theta){
    if(distr == "lognormal"){
      do.call(stat,list(x = dt, distr = distr, meanlog = theta[1], sdlog = theta[2]))
    }else if(distr == "weibull"){
      do.call(stat,list(x = dt, distr = distr, shape = theta[1], scale = theta[2]))
    }
  }
					 
	#/PSO parameter
  #//inertia weight, SA Type (Orkcu et al., 2015)
	omegamin <- .4
	omegamax <- .9
	it <- 1:max.iter
	omega <- omegamin+((omegamax-omegamin)*.95^(it-1))

	#/Initialization
	#//problem dimension
	p <- 2
	#//particle's position
	x <- matrix(runif(n.swarm*p,limit[1],limit[2]),nrow=n.swarm)
	#//particle's velocity
	v <- matrix(0,nrow=n.swarm,ncol=p)
  #//evaluate particle
	eval <- apply(x,1,FUN = fn.fit)
	it <- 0
	Pbest <- x
	besteval <- eval #personal best
	Gbest <- x[which.min(besteval),] #global best
	fneval.step <- min(besteval)
	#//initialization of looping condition
	loop <- TRUE

	#/iteration process
	while(loop){
		it <- it+1
		#//particle update
		r1 <- runif(p); r2 <- runif(p)
		for(j in 1:n.swarm){
		  x_prev <- x[j,]
			v[j,] <- omega[it]*v[j,]+cf[1]*r1*(Pbest[j,]-x[j,])+
					cf[2]*r2*(Gbest-x[j,])
			x[j,] <- x[j,]+v[j,]
				
			temp <- x[j,]<limit[1]
			if(any(temp)){
				x[j,temp] <- x_prev[temp]
			}

			temp <- x[j,]>limit[2]
			if(any(temp)){
				x[j,temp] <- x_prev[temp]
			}
			
			#///evaluate particle
			eval[j] <- fn.fit(x[j,])
		}
		
		Pbest.idx <- which(eval<besteval)
		besteval[Pbest.idx] <- eval[Pbest.idx]
		Pbest[Pbest.idx,] <- x[Pbest.idx,]
		Gbest.idx <- which.min(besteval)
		Gbest <- Pbest[Gbest.idx,]
		fneval <- besteval[Gbest.idx]
		fneval.step <- c(fneval.step,fneval)
		loop <- it != max.iter
	}
	
	#/labelling result & calculate log likelihood
  if(distr == "lognormal"){
	  names(Gbest) <- c("meanlog","sdlog")
	}else if (distr == "weibull"){
	  names(Gbest) <- c("shape","scale")
	}
	llh <- loglik(dt,distr,Gbest) #log-likelihood function from stat.r source
	
	#/se of estimators
	#se <- sqrt(diag(solve(pracma::hessian(f=fnloglik,x0=Gbest,x=dt))))
	#names(se) <- names(Gbest)

	#/output
	names(fneval.step) <- c(0,1:max.iter)
	output <- list(solution=Gbest,stat=stat,stat.value=fneval,
	               loglik=llh,iter=it,stat.trace=fneval.step)
	return(output)
}