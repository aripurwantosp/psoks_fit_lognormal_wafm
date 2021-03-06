#PSO-KS Empirical Study Fromm IDHS 2017 Data
#Fitting Lognormal Distribution to Women Age at First Marriage Data
#
#Ari Purwanto Sarwo Prasojo & Puguh Prasetyoputra (2019)
#Research Center for Population, Indonesian Institute of Sciences
#________________________________________________________________

#Load package, source, data ----
source("macro/stat.R")
source("macro/pso.fitdist.R")
library(fitdistrplus)
library(readstata13)
library(dplyr)
library(ggplot2)

#/function to calculate lognormal properties (mode, mean, variance)
lnorm_prop <- function(theta){
  mode <- exp(theta[1]-theta[2]^2)
  mean <- exp(theta[1]+.5*theta[2]^2)
  var <- (exp(theta[2]^2)-1)*exp(2*theta[1]+theta[2]^2)
  centr <- c(mode,mean,var)
  names(centr) <- c("mode","mean","var")
  return(centr)
}

#/read data form stata format
evwomen <- read.dta13("sdki17-evwomen.dta")
# View(evwomen)
# str(evwomen)


# National Aggregate ------------------------------------------------------

#/data sample ----
all_afm <- filter(evwomen, V501 != "Never in union")
all_afm <- all_afm$V511

#/fit distribution ----
psoks_all_afm <- pso.fitdist(dt = all_afm, dist = "lognormal", stat = "ks.stat",
                             limit = c(0,100), max.iter = 200, n.swarm = 20)
#//pso behaviour
pso_beh <- data.frame(iter = 0:200, ks = psoks_all_afm$stat.trace)
ggplot(pso_beh, aes(iter, ks)) + geom_line(color = "brown1") + geom_point(color = "brown1") + 
  labs(title = "PSO-KS's Behaviour, Fitting Lognormal Distribution", 
       subtitle = "Women's age at first marriage, national aggregation",
       x = "Iteration", y = "KS Distance") +
  theme_bw()

mle_all_afm <- fitdist(all_afm, distr = "lnorm", method = "mle")
mme_all_afm <- fitdist(all_afm, distr = "lnorm", method = "mme")

est_all_afm <- rbind(psoks_all_afm$solution, mle_all_afm$estimate, mme_all_afm$estimate)
colnames(est_all_afm) <- c("meanlog","sdlog")
rownames(est_all_afm) <- c("PSOKS","MLE","MME")
cat("Parameter estimation, national\n")
cat("---------------------------------------------\n")
print(est_all_afm)

#/properties ----
all_prop <- rbind(lnorm_prop(psoks_all_afm$solution),
                     lnorm_prop(mle_all_afm$estimate),
                     lnorm_prop(mme_all_afm$estimate))
rownames(all_prop) <- c("PSOKS","MLE","MME")
cat("Properties, national\n")
cat("---------------------------------------------\n")
print(all_prop)

#/goodness of fit statistcs ----
#//mse of distribution fitting
mse_psoks_all_afm <- mse.stat(all_afm,"lognormal", meanlog = psoks_all_afm$solution[1],
                              sdlog = psoks_all_afm$solution[2])
mse_mle_all_afm <- mse.stat(all_afm,"lognormal", meanlog = mle_all_afm$estimate[1],
                            sdlog = mle_all_afm$estimate[2])
mse_mme_all_afm <- mse.stat(all_afm,"lognormal", meanlog = mme_all_afm$estimate[1],
                            sdlog = mme_all_afm$estimate[2])

#//ks distance
ks_psoks_all_afm <- ks.stat(all_afm,"lognormal", meanlog = psoks_all_afm$solution[1],
                              sdlog = psoks_all_afm$solution[2])
ks_mle_all_afm <- ks.stat(all_afm,"lognormal", meanlog = mle_all_afm$estimate[1],
                            sdlog = mle_all_afm$estimate[2])
ks_mme_all_afm <- ks.stat(all_afm,"lognormal", meanlog = mme_all_afm$estimate[1],
                            sdlog = mme_all_afm$estimate[2])

#//log-likelihood
loglik_psoks_all_afm <- loglik(all_afm, "lognormal", psoks_all_afm$solution)
loglik_mle_all_afm <- loglik(all_afm, "lognormal", mle_all_afm$estimate)
loglik_mme_all_afm <- loglik(all_afm, "lognormal", mme_all_afm$estimate)

stat_all_afm <- rbind(c(mse_psoks_all_afm, mse_mle_all_afm, mse_mme_all_afm), #mse
                      c(ks_psoks_all_afm, ks_mle_all_afm, ks_mme_all_afm), #ks
                      c(loglik_psoks_all_afm, loglik_mle_all_afm, loglik_mme_all_afm)) #loglik
colnames(stat_all_afm) <- c("PSOKS","MLE","MME")
rownames(stat_all_afm) <- c("mse dist","ks distance","loglik")
cat("Goodness of fit measure, national aggregation\n")
cat("---------------------------------------------\n")
print(stat_all_afm)


#/graph fitting ----
#//histogram & pdf
len_all_afm <- length(all_afm)
den_all_psoks <- dlnorm(all_afm, meanlog = psoks_all_afm$solution[1],
                          sdlog = psoks_all_afm$solution[2])
den_all_mle <- dlnorm(all_afm, meanlog = mle_all_afm$estimate[1],
                        sdlog = mle_all_afm$estimate[2])
den_all_mme <- dlnorm(all_afm, meanlog = mme_all_afm$estimate[1],
                        sdlog = mme_all_afm$estimate[2])
den_all_afm <- data.frame(afm = rep(all_afm,3),
                         Fit = c(rep("PSO-KS", len_all_afm), rep("MLE", len_all_afm),
                                 rep("MME", len_all_afm)),
                         Density = c(den_all_psoks, den_all_mle, den_all_mme))

all_hist <- hist(all_afm, freq = FALSE)
ggplot(den_all_afm, aes(x = afm)) +
  geom_histogram(aes(y = ..density..), breaks = all_hist$breaks,
                 fill = "gray70", color = "white") +
  xlim(min(all_afm),max(all_afm)) +
  geom_line(aes(x = afm, y = Density, color = Fit), size = .8) +
  labs(title = "Histogram & Lognormal PDF", 
       subtitle = "Women's age at first marriage, national",
       x = "Age", y = "Density") +
  theme_bw()

#//ecdf & cdf
cumden_all_psoks <- plnorm(all_afm, meanlog = psoks_all_afm$solution[1],
                        sdlog = psoks_all_afm$solution[2])
cumden_all_mle <- plnorm(all_afm, meanlog = mle_all_afm$estimate[1],
                      sdlog = mle_all_afm$estimate[2])
cumden_all_mme <- plnorm(all_afm, meanlog = mme_all_afm$estimate[1],
                      sdlog = mme_all_afm$estimate[2])
cumden_all_afm <- data.frame(afm = rep(all_afm,3),
                          Fit = c(rep("PSO-KS", len_all_afm), rep("MLE", len_all_afm),
                                  rep("MME", len_all_afm)),
                          Cumden = c(cumden_all_psoks, cumden_all_mle, cumden_all_mme))

ggplot(data = cumden_all_afm, aes(afm))+
  stat_ecdf(geom = "step") +
  geom_line(aes(x = afm, y = Cumden, color = Fit), size = .7) +
  labs(title = "Empirical & Lognormal CDF", 
       subtitle = "Women's age at first marriage, national",
       x = "Age", y = "Cumulative Density") +
  theme_bw()


# Diassagregation by Residence (Rural/Urban) ---------------------------------

#/data sample ----
rural_afm <- filter(evwomen, V501 != "Never in union" & V025 == "Rural")$V511
urban_afm <- filter(evwomen, V501 != "Never in union" & V025 == "Urban")$V511

#/fit distribution ----
psoks_rural_afm <- pso.fitdist(dt = rural_afm, dist = "lognormal", stat = "ks.stat",
                             limit = c(0,100), max.iter = 200, n.swarm = 20)
psoks_urban_afm <- pso.fitdist(dt = urban_afm, dist = "lognormal", stat = "ks.stat",
                               limit = c(0,100), max.iter = 200, n.swarm = 20)
#//pso behaviour
pso_beh_reg <- data.frame(iter = rep(0:200,2), Residence = c(rep("Rural", 201),rep("Urban", 201)),
                          ks = c(psoks_rural_afm$stat.trace, psoks_urban_afm$stat.trace))
ggplot(pso_beh_reg, aes(iter, ks, group = Residence)) + geom_line(aes(color = Residence)) +
  geom_point(aes(color = Residence)) + 
  labs(title = "PSO-KS's Behaviour, Fitting Lognormal Distribution", 
       subtitle = "Women's age at first marriage, national by residence",
       x = "Iteration", y = "KS Distance") +
  theme_bw()

mle_rural_afm <- fitdist(rural_afm, distr = "lnorm", method = "mle")
mle_urban_afm <- fitdist(urban_afm, distr = "lnorm", method = "mle")
mme_rural_afm <- fitdist(rural_afm, distr = "lnorm", method = "mme")
mme_urban_afm <- fitdist(urban_afm, distr = "lnorm", method = "mme")

est_rural_afm <- rbind(psoks_rural_afm$solution, mle_rural_afm$estimate, mme_rural_afm$estimate)
colnames(est_rural_afm) <- c("meanlog","sdlog")
rownames(est_rural_afm) <- c("PSOKS","MLE","MME")
cat("Parameter estimation, national-rural region\n")
cat("---------------------------------------------\n")
print(est_rural_afm)

est_urban_afm <- rbind(psoks_urban_afm$solution, mle_urban_afm$estimate, mme_urban_afm$estimate)
colnames(est_urban_afm) <- c("meanlog","sdlog")
rownames(est_urban_afm) <- c("PSOKS","MLE","MME")
cat("Parameter estimation, national-urban region\n")
cat("---------------------------------------------\n")
print(est_urban_afm)

#/properties ----
rural_prop <- rbind(lnorm_prop(psoks_rural_afm$solution),
                     lnorm_prop(mle_rural_afm$estimate),
                     lnorm_prop(mme_rural_afm$estimate))
rownames(rural_prop) <- c("PSOKS","MLE","MME")
cat("Properties, national-rural region\n")
cat("---------------------------------------------\n")
print(rural_prop)

urban_prop <- rbind(lnorm_prop(psoks_urban_afm$solution),
                       lnorm_prop(mle_urban_afm$estimate),
                       lnorm_prop(mme_urban_afm$estimate))
rownames(urban_prop) <- c("PSOKS","MLE","MME")
cat("Properties, national-urban region\n")
cat("---------------------------------------------\n")
print(urban_prop)

#/goodness of fit statistcs ----
#//mse of distribution fitting
mse_psoks_rural_afm <- mse.stat(rural_afm,"lognormal", meanlog = psoks_rural_afm$solution[1],
                              sdlog = psoks_rural_afm$solution[2])
mse_mle_rural_afm <- mse.stat(rural_afm,"lognormal", meanlog = mle_rural_afm$estimate[1],
                            sdlog = mle_rural_afm$estimate[2])
mse_mme_rural_afm <- mse.stat(rural_afm,"lognormal", meanlog = mme_rural_afm$estimate[1],
                            sdlog = mme_rural_afm$estimate[2])

mse_psoks_urban_afm <- mse.stat(urban_afm,"lognormal", meanlog = psoks_urban_afm$solution[1],
                                sdlog = psoks_urban_afm$solution[2])
mse_mle_urban_afm <- mse.stat(urban_afm,"lognormal", meanlog = mle_urban_afm$estimate[1],
                              sdlog = mle_urban_afm$estimate[2])
mse_mme_urban_afm <- mse.stat(urban_afm,"lognormal", meanlog = mme_urban_afm$estimate[1],
                              sdlog = mme_urban_afm$estimate[2])

#//ks distance
ks_psoks_rural_afm <- ks.stat(rural_afm,"lognormal", meanlog = psoks_rural_afm$solution[1],
                            sdlog = psoks_rural_afm$solution[2])
ks_mle_rural_afm <- ks.stat(rural_afm,"lognormal", meanlog = mle_rural_afm$estimate[1],
                          sdlog = mle_rural_afm$estimate[2])
ks_mme_rural_afm <- ks.stat(rural_afm,"lognormal", meanlog = mme_rural_afm$estimate[1],
                          sdlog = mme_rural_afm$estimate[2])

ks_psoks_urban_afm <- ks.stat(urban_afm,"lognormal", meanlog = psoks_urban_afm$solution[1],
                              sdlog = psoks_urban_afm$solution[2])
ks_mle_urban_afm <- ks.stat(urban_afm,"lognormal", meanlog = mle_urban_afm$estimate[1],
                            sdlog = mle_urban_afm$estimate[2])
ks_mme_urban_afm <- ks.stat(urban_afm,"lognormal", meanlog = mme_urban_afm$estimate[1],
                            sdlog = mme_urban_afm$estimate[2])

#//log-likelihood
loglik_psoks_rural_afm <- loglik(rural_afm, "lognormal", psoks_rural_afm$solution)
loglik_mle_rural_afm <- loglik(rural_afm, "lognormal", mle_rural_afm$estimate)
loglik_mme_rural_afm <- loglik(rural_afm, "lognormal", mme_rural_afm$estimate)

loglik_psoks_urban_afm <- loglik(urban_afm, "lognormal", psoks_urban_afm$solution)
loglik_mle_urban_afm <- loglik(urban_afm, "lognormal", mle_urban_afm$estimate)
loglik_mme_urban_afm <- loglik(urban_afm, "lognormal", mme_urban_afm$estimate)

stat_rural_afm <- rbind(c(mse_psoks_rural_afm, mse_mle_rural_afm, mse_mme_rural_afm), #mse
                      c(ks_psoks_rural_afm, ks_mle_rural_afm, ks_mme_rural_afm), #ks
                      c(loglik_psoks_rural_afm, loglik_mle_rural_afm, loglik_mme_rural_afm)) #loglik

stat_urban_afm <- rbind(c(mse_psoks_urban_afm, mse_mle_urban_afm, mse_mme_urban_afm), #mse
                        c(ks_psoks_urban_afm, ks_mle_urban_afm, ks_mme_urban_afm), #ks
                        c(loglik_psoks_urban_afm, loglik_mle_urban_afm, loglik_mme_urban_afm)) #loglik

colnames(stat_rural_afm) <- colnames(stat_urban_afm) <- c("PSOKS","MLE","MME")
rownames(stat_rural_afm) <- rownames(stat_urban_afm) <- c("mse dist","ks distance","loglik")
cat("Goodness of fit measure, national-rural\n")
cat("---------------------------------------------\n")
print(stat_rural_afm)

cat("Goodness of fit measure, national-urban\n")
cat("---------------------------------------------\n")
print(stat_urban_afm)


#/graph fitting ----
#//histogram & pdf
reg_hist <- hist(rural_afm)
len_rural_afm <- length(rural_afm)
len_urban_afm <- length(urban_afm)
den_rural_psoks <- dlnorm(rural_afm, meanlog = psoks_rural_afm$solution[1],
                          sdlog = psoks_rural_afm$solution[2])
den_rural_mle <- dlnorm(rural_afm, meanlog = mle_rural_afm$estimate[1],
                        sdlog = mle_rural_afm$estimate[2])
den_rural_mme <- dlnorm(rural_afm, meanlog = mme_rural_afm$estimate[1],
                        sdlog = mme_rural_afm$estimate[2])
den_urban_psoks <- dlnorm(urban_afm, meanlog = psoks_urban_afm$solution[1],
                          sdlog = psoks_urban_afm$solution[2])
den_urban_mle <- dlnorm(urban_afm, meanlog = mle_urban_afm$estimate[1],
                        sdlog = mle_urban_afm$estimate[2])
den_urban_mme <- dlnorm(urban_afm, meanlog = mme_urban_afm$estimate[1],
                        sdlog = mme_urban_afm$estimate[2])
reg_df <- data.frame(afm = c(rep(rural_afm, 3), rep(urban_afm, 3)),
                     Residence = c(rep("Rural", 3*len_rural_afm), rep("Urban", 3*len_urban_afm)),
                     Fit = c(rep("PSO-KS",len_rural_afm), rep("MLE",len_rural_afm),
                             rep("MME",len_rural_afm),
                             rep("PSO-KS",len_urban_afm), rep("MLE",len_urban_afm),
                             rep("MME",len_urban_afm)),
                     Density = c(den_rural_psoks, den_rural_mle, den_rural_mme,
                                 den_urban_psoks, den_urban_mle, den_urban_mme))

ggplot(data = reg_df)+
  geom_histogram(aes(x = afm, y = ..density..),
                 fill = "gray70", color = "white", breaks = reg_hist$breaks) +
  xlim(min(all_afm),max(all_afm)) +
  geom_line(aes(x = afm, y = Density, color = Fit), size = .8) +
  labs(title = "Histogram & Lognormal PDF", 
       subtitle = "Women's age at first marriage, national by residence",
       x = "Age", y = "Density") +
  theme_bw() +
  facet_grid(Residence ~.)


#//ecdf & cdf
cum_rural_psoks <- plnorm(rural_afm, meanlog = psoks_rural_afm$solution[1],
                          sdlog = psoks_rural_afm$solution[2])
cum_rural_mle <- plnorm(rural_afm, meanlog = mle_rural_afm$estimate[1],
                        sdlog = mle_rural_afm$estimate[2])
cum_rural_mme <- plnorm(rural_afm, meanlog = mme_rural_afm$estimate[1],
                        sdlog = mme_rural_afm$estimate[2])
cum_urban_psoks <- plnorm(urban_afm, meanlog = psoks_urban_afm$solution[1],
                          sdlog = psoks_urban_afm$solution[2])
cum_urban_mle <- plnorm(urban_afm, meanlog = mle_urban_afm$estimate[1],
                        sdlog = mle_urban_afm$estimate[2])
cum_urban_mme <- plnorm(urban_afm, meanlog = mme_urban_afm$estimate[1],
                        sdlog = mme_urban_afm$estimate[2])
reg_cumd <- data.frame(afm = c(rep(rural_afm, 3), rep(urban_afm, 3)),
                     Residence = c(rep("Rural", 3*len_rural_afm), rep("Urban", 3*len_urban_afm)),
                     Fit = c(rep("PSO-KS",len_rural_afm), rep("MLE",len_rural_afm),
                             rep("MME",len_rural_afm),
                             rep("PSO-KS",len_urban_afm), rep("MLE",len_urban_afm),
                             rep("MME",len_urban_afm)),
                     Cumden = c(cum_rural_psoks, cum_rural_mle, cum_rural_mme,
                                 cum_urban_psoks, cum_urban_mle, cum_urban_mme))

ggplot(data = reg_cumd, aes(afm))+
  stat_ecdf(geom = "step") +
  geom_line(aes(x = afm, y = Cumden, color = Fit), size = .7) +
  labs(title = "Empirical & Lognormal CDF", 
       subtitle = "Women's age at first marriage, national by residence",
       x = "Age", y = "Cumulative Density") +
  theme_bw() +
  facet_grid(Residence ~.)


# Properties --------------------------------------------------------------

prop <- function(theta){
  meanlog <- theta[1]
  sdlog <- theta[2]
  
  mode <- exp(meanlog-sdlog^2)
  mean <- exp(meanlog + .5*sdlog^2)
  vari <- (exp(sdlog^2)-1)*exp(2*meanlog + sdlog^2)
  proper <- c(mode, mean, vari)
  names(proper) <- c("mode", "mean", "variance")
  return(proper)
}

# t <- rlnorm(100, meanlog = 1, sdlog = .5)
# Mode(t)
# mean(t)
# var(t)
# fit <- fitdist(t,"lnorm","mle")
# prop(fit$estimate)


#National
nat_psoks <- c(2.994, 0.213)
prop(nat_psoks)
nat_mle <- c(2.997, 0.221)
prop(nat_mle)
nat_mme <- c(2.997, 0.223)
prop(nat_mme)

#Reside
#/rural
rur_psoks <- c(2.949, 0.201)
prop(rur_psoks)
rur_mle <- c(2.956, 0.220)
prop(rur_mle)
rur_mme <- c(2.955, 0.224)
prop(rur_mme)

#/urban
ur_psoks <- c(3.037, 0.209)
prop(ur_psoks)
ur_mle <- c(3.038, 0.215)
prop(ur_mle)
ur_mme <- c(3.038, 0.215)
prop(ur_mme)