# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script includes an example of SEM/occupancy models:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(latex2exp)
library(reshape2)
library(jagsUI)

set.seed(1234)
n <- 300
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate our two covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
forest <- rnorm(n, 0, 1)
snow <- rnorm(n, 0, 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate bobcat occupancy probability (psiB)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha0 <- -1
alpha1 <- -3
psiB <- plogis(alpha0 + alpha1 * snow)
plot(psiB ~ snow) # bobcats aren't present in deep snow

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate bobcat occupancy (zB)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zB <- rbinom(n, 1, psiB)
table(zB)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate lynx occupancy probability (psiL)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta0 <- 0
beta1 <- 0.5
beta2 <- -2
psiL <- plogis(beta0 + beta1 * forest + beta2 * zB)
plot(psiL ~ forest) 
# these two lines of points represent lynx occupancy probabilities at sites
# with (upper) and without (lower) bocats

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate lynx occupancy (zL)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zL <- rbinom(n, 1, psiL)
table(zL)

table(zL,zB)
# 97 sites have neither species, only 19 have both


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate detection data for bobcat (yB) and lynx (yL)
# we'll sample each site for four weeks (K)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
K <- 4
yB <- rep(0, n)
yL <- rep(0, n)
p <- 0.5
for (i in 1:n){
  yB[i] <- rbinom(1, K * zB[i], p)
  yL[i] <- rbinom(1, K * zL[i], p)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(yB,zB)
table(yL,zL) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("cat_occupancy.jags")
cat("
      model {

      alpha0 ~ dlogis(0, 1)
      alpha1 ~ dnorm(0, 0.1)

      beta0 ~ dlogis(0, 1)
      beta1 ~ dnorm(0, 0.1)
      beta2 ~ dnorm(0, 0.1)
      
      pB ~ dbeta(1, 1)
      pL ~ dbeta(1, 1)


      for (i in 1:n){
      
        logit(psiB[i]) = alpha0 + alpha1 * s[i]
        logit(psiL[i]) = beta0 + beta1 * f[i] + beta2 * zB[i]
        # note we use zB as a covariate instead of psiB because we care 
        # about the effect of bobcat presence (zB)...
        # not the effect of bocat probability of presence (psiB)
        
        zB[i] ~ dbern(psiB[i])
        zL[i] ~ dbern(psiL[i])
        
        yB[i] ~ dbin(pB, zB[i] * K)
        yL[i] ~ dbin(pL, zL[i] * K)
        
      }

      
      }
      ",fill = TRUE)
sink()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bundle data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zB.dat <- rep(NA, n); zB.dat[yB > 0] <- 1
zL.dat <- rep(NA, n); zL.dat[yL > 0] <- 1


jags.data <- list(n = n, s = snow, f = forest,
                  zB = zB.dat, zL = zL.dat, K = K,
                  yB = yB, yL = yL)

# Initial values
inits <- function(){list(alpha = c(0,0))}  

# Parameters monitored
parameters <- c('alpha0','alpha1',
                'beta0','beta1','beta2',
                'pB','pL')

nc <- 4
nt <- 10
ni <- 25000
nb <- 10000

library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "cat_occupancy.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()
print(m, digits = 3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
#
# 1) Make a plot of the effect of bobcat occupancy on lynx occupancy
#
# 2) Work with a neighbor to think about a different model parameterization...
#      i.e., you might include a composite or latent variable rather
#            than effects among measured variable
#    
# 3) revise the code above the model to strengthen or weaken the effects
#    of bobcat on lynx. What happens to the raw 'contingency tables'
#    i.e., table(zB,zL)
#    Why?
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
