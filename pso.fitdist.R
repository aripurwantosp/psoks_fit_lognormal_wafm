pso.fitdist <- function(dt,distr,stat,limit,cf=c(2,2),max.iter=100,n.swarm=30){
##----------------------------------------------------------------------
##kspso.optim
#
##Deskripsi :
#Macro ini digunakan untuk optimasi fungsi objektif (kasus maksimisasi)
#
##Cara menggunakan :
#ketik optimpso(param,fn,...,limit,inertia,cf,max.iter,n.swarm), dimana:
#param			: inisial parameter (isikan NA sebanyak ukuran dimensi)
#fn         	: fungsi objektif (fitness)
#...         	: argumen tambahan pada fungsi objektif
#limit        	: batas-batas pencarian c(batas bawah, batas atas)
#cf				: faktor kognitif dan sosial (default (c(2,2))
#max.iter     	: iterasi maksimum (default 100)
#n.swarm      	: ukuran populasi
#
##Output
#solution     	: solusi optimasi
#fitness      	: nilai objektif
#iteration    	: banyaknya iterasi
#fitness.step 	: record nilai objektif setiap iterasi
##----------------------------------------------------------------------

	
#   #ks as fitness function
# 	ks.fit <- function(param){
#     if(distr == "lognormal"){
# 	    ks.stat(dt,distr,meanlog = param[1],sdlog = param[2])
# 	  }else if (distr == "weibull"){
# 	    ks.stat(dt,distr,shape = param[1],sdlog = param[2])
# 	  }
# 	}
	
	# #fitness function
	# if(method == "ks"){
	#   fn.fit <- function(param){
	#     if(distr == "lognormal"){
	#       ks.stat(dt,distr,meanlog = param[1],sdlog = param[2])
	#     }else if (distr == "weibull"){
	#       ks.stat(dt,distr,shape = param[1],sdlog = param[2])
	#     }
	#   }
	# }else if(method == "mse"){
	#   function(param){
	#     if(distr == "lognormal"){
	#       mse.stat(dt,distr,meanlog = param[1],sdlog = param[2])
	#     }else if (distr == "weibull"){
	#       mse.stat(dt,distr,shape = param[1],sdlog = param[2])
	#     }
	#   }
	# }
  
  #fitness function
  fn.fit <- function(theta){
    if(distr == "lognormal"){
      do.call(stat,list(x = dt, distr = distr, meanlog = theta[1], sdlog = theta[2]))
    }else if(distr == "weibull"){
      do.call(stat,list(x = dt, distr = distr, shape = theta[1], scale = theta[2]))
    }
  }
					 
	#PSO parameter
	omegamin <- .4
	omegamax <- .9

	#Inisialisasi
	#dimensi pencarian
	p <- 2
	#posisi partikel
	x <- matrix(runif(n.swarm*p,limit[1],limit[2]),nrow=n.swarm)
	#kecepatan partikel
	v <- matrix(0,nrow=n.swarm,ncol=p)
	eval <- 0

	#bobot inersia
	it <- 1:max.iter
	omega <- omegamin+((omegamax-omegamin)*.95^(it-1))

	#eval <- apply(x,1,FUN = ks.fit)
	eval <- apply(x,1,FUN = fn.fit)
	it <- 0
	Pbest <- x
	besteval <- eval
	Gbest <- x[which.min(besteval),]
	fneval.step <- min(besteval)
	loop <- TRUE

	#Proses Iterasi
	while(loop){
		it <- it+1
		#update partikel
		r1 <- runif(p); r2 <- runif(p)
		for(j in 1:n.swarm){
		  x_prev <- x[j,]
			v[j,] <- omega[it]*v[j,]+cf[1]*r1*(Pbest[j,]-x[j,])+
					cf[2]*r2*(Gbest-x[j,])
			x[j,] <- x[j,]+v[j,]
				
			temp <- x[j,]<limit[1]
			if(any(temp)){
				x[j,temp] <- x_prev[temp]
				#v[j,temp] <- 0
			}

			temp <- x[j,]>limit[2]
			if(any(temp)){
				x[j,temp] <- x_prev[temp]
				#v[j,temp] <- 0
			}
			
			#evaluasi partikel
			#eval[j] <- ks.fit(x[j,])
			eval[j] <- fn.fit(x[j,])
		}
		
		index1 <- which(eval<besteval)
		besteval[index1] <- eval[index1]
		Pbest[index1,] <- x[index1,]
		index2 <- which.min(besteval)
		Gbest <- Pbest[index2,]
		fneval <- besteval[index2]
		fneval.step <- c(fneval.step,fneval)
		loop <- it!=max.iter
	}
	
	#labelling result & calc log likelihood
  if(distr == "lognormal"){
	  names(Gbest) <- c("meanlog","sdlog")
	  fnloglik <- function(x,theta){
	    sum(dlnorm(x,meanlog = theta[1],sdlog = theta[2],log = TRUE))
	  }
	}else if (distr == "weibull"){
	  names(Gbest) <- c("shape","scale")
	  fnloglik <- function(x,theta){
	    sum(dweibull(x,shape = theta[1],scale = theta[2],log = TRUE))
	  }
	}
	loglik <- fnloglik(dt,Gbest)
	#se <- sqrt(diag(solve(pracma::hessian(f=fnloglik,x0=Gbest,x=dt))))
	#names(se) <- names(Gbest)

	#output
	names(fneval.step) <- c(0,1:max.iter)
	output <- list(solution=Gbest,stat=stat,stat.value=fneval,
	               loglik=loglik,iter=it, ks.trace=fneval.step)
	return(output)
}