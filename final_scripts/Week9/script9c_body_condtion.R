library(blavaan)
library(jagsUI)
library(reshape2)
library(MASS)
library(latex2exp)


set.seed(1234)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# creating a covariance matrix using the SRS separation strategy following: 
# Barnard et al. (2000) Statistical Sinica & 
# Alvarez, Niemi, and Simpson (2014)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# variances of each measured trait
sigma <- diag(c(1,0.5,1,1))
# correlation matrix
R <- matrix(c(1,    0.8, -0.8,  0.8,
              0.8,  1,   -0.8,  0.8,
              -0.8, -0.8, 1,    -0.8,
              0.8,  0.8,  -0.8, 1), byrow = T, nrow = 4)
# variance-covariance matrix
Sigma <- sigma %*% R %*% sigma


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 200
mu <- rep(0,4)
y <- mvrnorm(n, mu, Sigma)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# normal v. bivariate normal distribution figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot a normal
hist(y[,1], main = '', xlab = 'y', breaks = 25, cex.lab = 2)

plot(y[,1] ~ y[,2], ylab = TeX("$y_1$"), xlab = TeX("$y_2"),
     cex.lab = 2, las = 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we're going to use these values (y) in our model, but will back-transform
# these z-standardized values to actual measurements
# m: mass (grams)
# f: fat (stays z-standardized, mean 0, sd = 1)
# c: cort (stays z-standardized, mean 0, sd = 1)
# e: eggs (number of eggs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mass <- rlnorm(n, 3 + 0.1 * y[,1], 0.01)
fat <- y[,2]
cort <- y[,3]
eggs <- round(exp(1.25 + 0.25 * y[,4]))

dat <- data.frame(mass = as.numeric(scale(mass)), fat, cort, eggs = log(eggs))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# panel plot data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,3))
plot(eggs ~ mass, las = 1)
plot(eggs ~ fat, las = 1)
plot(eggs ~ cort, las = 1)
plot(mass ~ fat, las = 1)
plot(cort ~ fat, las = 1)
plot(cort ~ mass, las = 1)







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's build two models
# model 1: latent variable model
# model 2: composite covariate model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("m1_latent.jags")
cat("
      model {

      for (j in 1:4){
        alpha[j] ~ dnorm(0,0.1)
        sigma[j] ~ dgamma(1,1)
        tau[j] = 1/(sigma[j] * sigma[j])
      }

      beta[1] = 1
      beta[2] ~ dnorm(0,0.1)
      beta[3] ~ dnorm(0,0.1)
      beta[4] ~ dnorm(0,0.1)

      for (i in 1:n){
      
        b[i] ~ dnorm(0, tau[4])
        
        # linear regression for fat and cort
        # note mass was z-standardized below
        m[i] ~ dnorm(alpha[1] + beta[1] * b[i], tau[1])
        f[i] ~ dnorm(alpha[2] + beta[2] * b[i], tau[2])
        c[i] ~ dnorm(alpha[3] + beta[3] * b[i], tau[3])
        
        # Poisson regression for number of eggs
        e[i] ~ dpois(exp(alpha[4] + beta[4] * b[i]))
        
      }

      }
      ",fill = TRUE)
sink()


sink("m2_composite.jags")
cat("
      model {

      alpha ~ dnorm(1,0.1)
      beta[1] = 1
      beta[2] ~ dnorm(0,1)
      beta[3] ~ dnorm(0,1)
      beta[4] ~ dnorm(0,1)
      sigma ~ dgamma(1,1)
      tau = pow(sigma, -2)
      
      for (i in 1:n){
        # body condition is a composite covarite with no variance      
        b[i] ~ dnorm(beta[1]*m[i] + beta[2]*f[i] + beta[3]*c[i], tau)
        
        # Poisson regression for number of eggs
        e[i] ~ dpois(exp(alpha + beta[4] * b[i]))
      }

      
      }
      ",fill = TRUE)
sink()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bundle data and run both models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jags.data <- list(n = n, m = as.numeric(scale(mass)), 
                  c = cort,
                  f = fat, e = eggs)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list(beta = c(NA, 1,-1, 1))}  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters monitored
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c('alpha','beta','sigma','b')

nc <- 4
nt <- 10
ni <- 25000
nb <- 10000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
Sys.time()
m1 <- jags(jags.data, inits, parameters, "m1_latent.jags", 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
          parallel = T)

m2 <- jags(jags.data, inits, parameters, "m2_composite.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()


print(m1, digits = 3)
print(m2, digits = 3)








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
# 1) Here are estimates of 'Body condition' from each model
par(mfrow = c(1,1))
plot(m1$q50$b ~ m2$q50$b, ylab = "latent 'b'", xlab = "composite 'b'")
#    1a) Why do they vary in scale?
#    1b) How does that affect beta[4] (i.e., the relationship between body condition and eggs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.iter <- length(m1$sims.list$alpha[,1])
res <- 100
lv.b <- seq(min(m1$q50$b),max(m1$q50$b),length.out = res)
cv.b <- seq(min(m2$q50$b),max(m2$q50$b),length.out = res)
lv.e <- matrix(NA, n.iter, res)
cv.e <- matrix(NA, n.iter, res)
qlv.e <- matrix(NA, res, 3)
qcv.e <- matrix(NA, res, 3)
for (j in 1:res){
  lv.e[,j] = exp(m1$sims.list$alpha[,4] + m1$sims.list$beta[,4] * lv.b[j])
  cv.e[,j] = exp(m2$sims.list$alpha + m2$sims.list$beta[,4] * cv.b[j])  
  qlv.e[j,] <- quantile(lv.e[,j], c(0.025,0.5,0.975))
  qcv.e[j,] <- quantile(cv.e[,j], c(0.025,0.5,0.975))  
}
par(mfrow = c(1,2))
plot(qlv.e[,2] ~ lv.b, type = 'l', ylim = c(1,7), lwd = 3,
     xlab = 'Latent condition', ylab = 'Expected clucth size')
lines(qlv.e[,1]~ lv.b, lty = 2, lwd = 2)
lines(qlv.e[,3]~ lv.b, lty = 2, lwd = 2)
plot(qcv.e[,2] ~ cv.b, type = 'l', ylim = c(1,7), lwd = 3,
     xlab = 'Composite condition', ylab = 'Expected clucth size')
lines(qcv.e[,1]~ cv.b, lty = 2, lwd = 2)
lines(qcv.e[,3]~ cv.b, lty = 2, lwd = 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
# 2) Make a similar plot for the relationship between clutch size 
#    and mass. Your first step for creating this relationship with the 
#    latent variable model should be creating a range of conditions (lv.b above)
#    Your first step for the composite covariate model should be creating a range
#    of masses (see jags.data$m). Why?
#
# 2b) Do these relationships differ? Why?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
# 3) Are intercepts necessary for the latent variable model (e.g., 
#    alpha[1], alpha[2], and alpha[3])? Why or why not?
#    
#    Run the model without them and see what happens?!
#    to do this, just delete alpha[1:3] from the relevant dnorm() models.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some rough analogous examples in blavaan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here are two models in blavaan, note there is no longer a Poisson distribution
# for the number of eggs, we just model each log-transformed clutch using
# an identity linke (normal distribution)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bsem1 <- bsem('eggs ~ 1 + condition
              condition =~ mass + fat + cort',
              data = dat)

bsem2 <- bsem('eggs ~ 1 + condition
              condition <~ 1*mass + fat + cort',
              data = dat)
summary(bsem1); m1$q50$beta; m1$q50$sigma
summary(bsem2)

