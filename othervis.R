#Supplementary Visualization
#Lognormal variation, ECDF, PSO-KS Behavior
#
#Ari Purwanto Sarwo Prasojo & Puguh Prasetyoputra (2019)
#Research Center for Population, Indonesian Institute of Sciences
#________________________________________________________________

library(ggplot2)
library(xlsx)

#Lognormal p.d.f ----

n <- 200
t <- seq(0,10,length = n)
mu <- c(0, 0.5, 1, 1)
sd <- c(1, 1, 1, 1.5)
lnorm <-  NULL
lnormlab <- NULL

for(i in 1:length(mu)){
  lnorm <- c(lnorm, dlnorm(t, meanlog = mu[i], sdlog = sd[i]))
  lnormlab <- c(lnormlab, rep(paste0("meanlog = ",mu[i],"; sdlog = ",sd[i]),n))
}

t <- rep(t, 4)

df <- data.frame("T" = t, Density = lnorm, Parameters = lnormlab)

ggplot(df, aes(T, Density, group = Parameters, color = Parameters))+
  geom_line(size = .725)+
  theme_bw()


#Empirical and theoretical cdf of lognormal ----

cdf <- plnorm(t, meanlog = 0, sdlog = 1)
df <- data.frame("T" = t, c.d.f = cdf)

dftext <- data.frame(x = c(2, 5), y = c(0.880, 0.860), 
                     label = c("Theoretical", "Empirical"),
                     color = c("red", "black"))

ggplot(df, aes(x = T))+
  stat_ecdf(geom = "step") +
  stat_ecdf(geom = "point") +
  xlim(c(0, max(t))) +
  stat_function(fun = plnorm, size = .7, args = list(meanlog = 0, sdlog = 1)) +
  geom_text(aes(x,y, label = label), data = dftext) +
  ylab("c.d.f") +
  theme_bw()


#PSO-KS behaviour ----

visstep <- read.xlsx("sim2/vispsostep.xlsx",2)
visstep$n <- as.factor(visstep$n )

ggplot(visstep, aes(Iteration, KS, group = n, color = n)) +
  geom_point() +
  geom_line() +
  ylab("KS Distance") +
  theme_bw()