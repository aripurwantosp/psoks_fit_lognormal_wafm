#Empirical Study Form IDHS 2017 Data
#Fitting Lognormal Distribution to Women Age at First Marriage Data

#Load package and source
source("stat.R")
source("pso.fitdist.R")
library(fitdistrplus)
library(readstata13)
library(dplyr)
library(ggplot2)

lnorm_meanvar <- function(theta){
  mode <- exp(theta[1]-theta[2]^2)
  mean <- exp(theta[1]+.5*theta[2]^2)
  var <- (exp(theta[2]^2)-1)*exp(2*theta[1]+theta[2]^2)
  centr <- c(mode,mean,var)
  names(centr) <- c("mode","mean","var")
  return(centr)
}

#Read data form stata format
evwomen <- read.dta13("sdki17-evwomen.dta")
View(evwomen)
str(evwomen)


# Aggregate ---------------------------------------------------------------

all_afm <- filter(evwomen, V501 != "Never in union")
all_afm <- all_afm$V511


#fitdist
psoks_all_afm <- pso.fitdist(dt = all_afm, dist = "lognormal", stat = "ks.stat",
                             limit = c(0,100), max.iter = 200, n.swarm = 20)
#/pso behaviour
pso_beh <- data.frame(iter = 0:200, ks = psoks_all_afm$ks.trace)
ggplot(pso_beh, aes(iter, ks)) + geom_line(color = "brown1") + geom_point(color = "brown1") + 
  labs(title = "PSO-KS's Behaviour, Fitting Lognormal Distribution", 
       subtitle = "Women age at first marriage, national aggregation",
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

#/central tendency
all_central <- rbind(lnorm_meanvar(psoks_all_afm$solution),
                     lnorm_meanvar(mle_all_afm$estimate),
                     lnorm_meanvar(mme_all_afm$estimate))
rownames(all_central) <- c("PSOKS","MLE","MME")
cat("Parameter estimation, national\n")
cat("---------------------------------------------\n")
print(all_central)

#goodness of fit statistcs
#/mse dsitr stat
mse_psoks_all_afm <- mse.stat(all_afm,"lognormal", meanlog = psoks_all_afm$solution[1],
                              sdlog = psoks_all_afm$solution[2])
mse_mle_all_afm <- mse.stat(all_afm,"lognormal", meanlog = mle_all_afm$estimate[1],
                            sdlog = mle_all_afm$estimate[2])
mse_mme_all_afm <- mse.stat(all_afm,"lognormal", meanlog = mme_all_afm$estimate[1],
                            sdlog = mme_all_afm$estimate[2])

#/ks stat
ks_psoks_all_afm <- ks.stat(all_afm,"lognormal", meanlog = psoks_all_afm$solution[1],
                              sdlog = psoks_all_afm$solution[2])
ks_mle_all_afm <- ks.stat(all_afm,"lognormal", meanlog = mle_all_afm$estimate[1],
                            sdlog = mle_all_afm$estimate[2])
ks_mme_all_afm <- ks.stat(all_afm,"lognormal", meanlog = mme_all_afm$estimate[1],
                            sdlog = mme_all_afm$estimate[2])

#/loglik
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


#graph fitting
#/histogram & pdf
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
       subtitle = "Women age at first marriage, national",
       x = "Age", y = "Density") +
  theme_bw()

#/ecdf & cdf
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
       subtitle = "Women age at first marriage, national",
       x = "Age", y = "Cumulative Density") +
  theme_bw()



ggplot(filter(evwomen, V501 != "Never in union"), aes(x = V511)) +
  stat_ecdf(geom = "step") +
  stat_function(aes(color = "PSO-KS"), fun = plnorm, size = .7,
                args = list(meanlog = psoks_all_afm$solution[1],
                            sdlog = psoks_all_afm$solution[2])) +
  stat_function(aes(color = "MLE"), fun = plnorm, size = .7,
                args = list(meanlog = mle_all_afm$estimate[1],
                            sdlog = mle_all_afm$estimate[2])) +
  stat_function(aes(color = "MME"), fun = plnorm, size = .7,
                args = list(meanlog = mme_all_afm$estimate[1],
                            sdlog = mme_all_afm$estimate[2])) +
  scale_colour_manual("Fit Lognormal", values = c("orange", "blue", "red"))+
  labs(title = "Empirical & Lognormal CDF", 
       subtitle = "Women age at first marriage, national",
       x = "Age", y = "Cumulative Density") +
  theme_minimal()


# Diassagregation by Region (Rural/Urban) ---------------------------------

rural_afm <- filter(evwomen, V501 != "Never in union" & V025 == "Rural")$V511
urban_afm <- filter(evwomen, V501 != "Never in union" & V025 == "Urban")$V511


#fitdist
psoks_rural_afm <- pso.fitdist(dt = rural_afm, dist = "lognormal", stat = "ks.stat",
                             limit = c(0,100), max.iter = 200, n.swarm = 20)
psoks_urban_afm <- pso.fitdist(dt = urban_afm, dist = "lognormal", stat = "ks.stat",
                               limit = c(0,100), max.iter = 200, n.swarm = 20)
#/pso behaviour
pso_beh_reg <- data.frame(iter = rep(0:200,2), region = c(rep("Rural", 201),rep("Urban", 201)),
                          ks = c(psoks_rural_afm$ks.trace, psoks_urban_afm$ks.trace))
ggplot(pso_beh_reg, aes(iter, ks, group = region)) + geom_line(aes(color = region)) +
  geom_point(aes(color = region)) + 
  labs(title = "PSO-KS's Behaviour, Fitting Lognormal Distribution", 
       subtitle = "Women age at first marriage, national by region",
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

#/central tendency
rural_central <- rbind(lnorm_meanvar(psoks_rural_afm$solution),
                     lnorm_meanvar(mle_rural_afm$estimate),
                     lnorm_meanvar(mme_rural_afm$estimate))
rownames(rural_central) <- c("PSOKS","MLE","MME")
cat("Parameter estimation, national-rural region\n")
cat("---------------------------------------------\n")
print(rural_central)

urban_central <- rbind(lnorm_meanvar(psoks_urban_afm$solution),
                       lnorm_meanvar(mle_urban_afm$estimate),
                       lnorm_meanvar(mme_urban_afm$estimate))
rownames(urban_central) <- c("PSOKS","MLE","MME")
cat("Parameter estimation, national-urban region\n")
cat("---------------------------------------------\n")
print(urban_central)


#graph fitting
#/histogram & pdf
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
                     Reside = c(rep("Rural", 3*len_rural_afm), rep("Urban", 3*len_urban_afm)),
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
       subtitle = "Women age at first marriage, national by reside",
       x = "Age", y = "Density") +
  theme_bw() +
  facet_grid(Reside ~.)


#/ecdf & cdf
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
                     Reside = c(rep("Rural", 3*len_rural_afm), rep("Urban", 3*len_urban_afm)),
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
       subtitle = "Women age at first marriage, national by reside",
       x = "Age", y = "Cumulative Density") +
  theme_bw() +
  facet_grid(Reside ~.)