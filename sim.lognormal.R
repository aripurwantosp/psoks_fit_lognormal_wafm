#Simulation Study for Lognormal Distribution

#Load packages and sources ----
source("pso.fitdist.R")
source("stat.R")
library(fitdistrplus)
library(xlsx)

#Global ----
nsamp <- c(50,100,200,500)
rep <- 1000

#Lognormal distribution
lnorm_par <- rbind(c(0.5,1),
                   c(1,1),
                   c(1,.5),
                   c(1,2))
# lnorm_par <- rbind(c(0,1),
#                    c(1,1))
colnames(lnorm_par) <- c("mulog","sd")

mn <- 0
mx <- 3

#Graph ----
# x <- seq(0,8,length = 200)
# for(i in 1:nrow(lnorm_par)){
#   par <- lnorm_par[i,]
#   denx <- dlnorm(x, meanlog = par[1], sdlog = par[2])
#   if(i == 1){
#     plot(x,denx,type="l",xlim=c(0,8),ylim=c(0,.7))
#   }else{
#     lines(x,denx)
#   }
# }

value <- c("mulog","sdlog")
meas <- c("mse_dist","ks","mle")
meas_ag <- c("bias_mulog","bias_sdlog","mse_mulog","mse_sdlog")
method <- c("PSO_KS_","MLE_","MME_")
namecol <- c("meanlog0","sdlog0","n",
             as.vector(sapply(method, FUN = paste0, value)))
nameass <- as.vector(sapply(method, FUN = paste0, meas))
nameass_ag <- as.vector(sapply(method, FUN = paste0, meas_ag))

namecolfit <- paste0("step",0:200)
namerow <- paste0("sim",1:rep)
meansim <- varsim <- ass_mean <- assess <- rectime<- namerectime <- NULL

for(i in 1:nrow(lnorm_par)){
  par <- lnorm_par[i,]
  
  for(n in nsamp){
    
    estimate_all <- matrix(NA, nrow = rep, ncol = 9)
    assess_all <- matrix(NA, nrow = rep, ncol = 9)
    psoks_fitrec <-  matrix(NA, nrow = rep, ncol=201)
    cat("\n\nlognormal: ",par[1],"-",par[2],"nsample: ",n,"-> ")
    
    start <- Sys.time()
    
    for(k in 1:rep){
      cat(k, " ")
      
      xrnd <- rlnorm(n, meanlog = par[1], sdlog = par[2])
      
      #Estimation and individual assessment
      #pso-ks
      psoks_est <- pso.fitdist(dt = xrnd, dist = "lognormal", stat = "ks.stat",
                               limit = c(mn,mx), max.iter = 200, n.swarm = 20)
      psoks <- psoks_est$solution
      psoks_mse_dist <- mse.stat(xrnd, distr = "lognormal", meanlog = psoks[1], sdlog = psoks[2])
      psoks_ks <- ks.stat(xrnd, distr = "lognormal", meanlog = psoks[1], sdlog = psoks[2])
      psoks_loglik <- loglik(xrnd, distr = "lognormal", theta = psoks)
      
      #mle estimation
      mle_est <- fitdist(xrnd, distr = "lnorm", method = "mle")
      mle <- mle_est$estimate
      mle_mse_dist <- mse.stat(xrnd, distr = "lognormal", meanlog = mle[1], sdlog = mle[2])
      mle_ks <- ks.stat(xrnd, distr = "lognormal", meanlog = mle[1], sdlog = mle[2])
      mle_loglik <- loglik(xrnd, distr = "lognormal", theta = mle)
      
      #mme estimation
      mme_est <- fitdist(xrnd, distr = "lnorm", method = "mme")
      mme <- mme_est$estimate
      mme_mse_dist <- mse.stat(xrnd, distr = "lognormal", meanlog = mme[1], sdlog = mme[2])
      mme_ks <- ks.stat(xrnd, distr = "lognormal", meanlog = mme[1], sdlog = mme[2])
      mme_loglik <- loglik(xrnd, distr = "lognormal", theta = mme)
      
      estimate_all[k, ] <- c(par, n, psoks, mle, mme)
      psoks_fitrec[k, ] <- psoks_est$ks.trace
      assess_all[k, ] <- c(psoks_mse_dist, psoks_ks, psoks_loglik, #psoks
                           mle_mse_dist, mle_ks, mle_loglik, #mle
                           mme_mse_dist, mme_ks, mme_loglik) #mme 
      
    }
    
    #Record time
    end <- Sys.time()
    time <- end-start
    rectime <- c(rectime, time)
    namerectime <- c(namerectime,paste0("lognormal ", par[1],"-",par[2]," nsample: ",n))
    
    #Save to excel for each scenario
    colnames(estimate_all) <- namecol
    rownames(estimate_all) <- namerow
    
    #psoks_fitrec <- rbind(psoks_fitrec, colMeans(psoks_fitrec))
    colnames(psoks_fitrec) <- namecolfit
    rownames(psoks_fitrec) <- namerow
    
    colnames(assess_all) <- nameass
    rownames(assess_all) <- namerow
    
    wb <- createWorkbook()
    tmpsheet <- createSheet(wb,sheetName = "results")
    addDataFrame(estimate_all,tmpsheet)
    
    tmpsheet <- createSheet(wb,sheetName = "psoks_step")
    addDataFrame(psoks_fitrec,tmpsheet)
    
    tmpsheet <- createSheet(wb,sheetName = "assess")
    addDataFrame(assess_all,tmpsheet)
    
    saveWorkbook(wb,paste0("lognormal_",par[1],"-",par[2],"_n",n,".xlsx"))
    
    
    
    #Record of aggregate values
    meanscen <- colMeans(estimate_all)
    meansim <- rbind(meansim, meanscen)
    varsim <- rbind(varsim, apply(estimate_all, 2, FUN = var))
    ass_mean <- rbind(ass_mean, colMeans(assess_all))
    
    #Aggregate assessment
    #pso-ks
    ass_bias_psoks <- meanscen[4:5]-par
    ass_mse_psoks <- c(sum((estimate_all[,4]-par[1])^2),
                       sum((estimate_all[,5]-par[2])^2))/rep

    #mle
    ass_bias_mle <- meanscen[6:7]-par
    ass_mse_mle <- c(sum((estimate_all[,6]-par[1])^2),
                     sum((estimate_all[,7]-par[2])^2))/rep
    
    #mme
    ass_bias_mme <- meanscen[8:9]-par
    ass_mse_mme <- c(sum((estimate_all[,8]-par[1])^2),
                     sum((estimate_all[,9]-par[2])^2))/rep
    
    #assessment record
    assess <- rbind(assess, c(ass_bias_psoks, ass_mse_psoks, #psoks
                              ass_bias_mle, ass_mse_mle, #mle
                              ass_bias_mme, ass_mse_mme #mme
    ))
    
  }
}

#Saving aggregate results to excel sheet
rectime <- data.frame(namerectime,rectime)

colnames(meansim) <- namecol
rownames(meansim) <- namerectime

colnames(varsim) <- namecol
rownames(varsim) <- namerectime

colnames(ass_mean) <- nameass
rownames(ass_mean) <- namerectime

colnames(assess) <- nameass_ag
rownames(assess) <- namerectime


wb <- createWorkbook()
tmpsheet <- createSheet(wb,sheetName="times")
addDataFrame(rectime,tmpsheet)

tmpsheet <- createSheet(wb,sheetName="aggregate_mean")
addDataFrame(meansim,tmpsheet)

tmpsheet <- createSheet(wb,sheetName="aggregate_var")
addDataFrame(varsim,tmpsheet)

tmpsheet <- createSheet(wb,sheetName="aggregate_ass")
addDataFrame(ass_mean,tmpsheet)

tmpsheet <- createSheet(wb,sheetName="assessment")
addDataFrame(assess,tmpsheet)

saveWorkbook(wb,"lognormal_sim_aggregate.xlsx")