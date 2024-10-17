library(blavaan)
library(lavaan)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample size
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 200

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate (collinear) covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x1 <- rnorm(n,0,1)
x2 <- x1 + rnorm(n, 0, 1)
plot(x1 ~ x2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parameter estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta1 <- 0.5
beta2 <- 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this is how we build a composite covariate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c <- beta1*x1 + beta2*x2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compare expected values of y
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ey.comp <- c
Ey.beta <- beta1*x1 + beta2*x2
# these are the same numbers... they should be... that's a good thing
plot(Ey.comp ~ Ey.beta)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate data with measurement/process error
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y <- rnorm(n, c, 0.1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(y ~ c)
plot(y ~ x1)
plot(y ~ x2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a data.frame with data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- data.frame(y,x1,x2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run model in lavaan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mL <- sem('y ~ c
           c <~ x1 + 1*x2',
           data = dat)
summary(mL) # great

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stan model via blavaan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mS <- bsem('y ~ c
           c <~ x1 + 1*x2',
           data = dat,
           mcmcfile = T,
           target = 'stan')
summary(mS) # great!

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JAGS model via blavaan, rut roh!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mJ <- bsem('y ~ c
           c <~ x1 + 1*x2',
          data = dat,
          mcmcfile = T,
          target = 'jags')
summary(mJ) # huh?!

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# revised JAGS model (a cool work-around suggested by
# the blavaan creator, E. Merkle, to create an almost identical model)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mR <- bsem('c =~ NA * y
           c ~ x1 + 1*x2
           c ~~ 0.0001*c
           y ~~ NA*y',
           data = dat,
           mcmcfile = T,
           target = 'jags')
summary(mR) 
# cool?! Check out the JAGS file to see how this is constructed
# and compare it to the model below...


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# So how do we actually do this in JAGS?!
# see below...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("composite_custom.jags")
cat("
      model {


      alpha ~ dnorm(0, 1)

      beta[1] ~ dnorm(0, 1)
      beta[2] = 1
      beta[3] ~ dnorm(0, 0.1)
      
      sigma ~ dgamma(1,1)
      tau = 1/(sigma * sigma)


      for (i in 1:n){
      
        c[i] = beta[1]*x1[i] + beta[2]*x2[i]  # our composite covariate has no error
        y[i] ~ dnorm(alpha + beta[3] * c[i], tau)
        
      }

      
      }
      ",fill = TRUE)
sink()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bundle data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jags.data <- list(y = y, x1 = x1, x2 = x2, n = n)

# Initial values
inits <- function(){list(alpha = 0)}  

# Parameters monitored
parameters <- c('alpha','beta','c','sigma')

nc <- 4
nt <- 10
ni <- 25000
nb <- 10000

library(jagsUI)
Sys.time()
mC <- jags(jags.data, inits, parameters, "composite_custom.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()
print(mC, digits = 3)

plot(mC$q50$c ~ c)



