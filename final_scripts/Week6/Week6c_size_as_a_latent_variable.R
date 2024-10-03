# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we'll use data based on some real brant data from Alaska
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(degreenet)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# s: Size, latent size
# c: Culmen, a bill length measurement in mm
# t: Tarsus, a leg length measurement in cm
# e: clutch size, or number of Eggs in the clutch
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 300
s <- rnorm(n, 0, 1)
par(mfrow = c(2,2), mar = c(5.1,5.1,2.1,2.1))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate number of eggs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha1 <- 4.15; beta1 <- 0.05
t <- rlnorm(n, alpha1 + beta1 * s, 0.025)
plot(t ~ s, ylab = 'Tarsus (mm)', xlab = 'Latent size', cex.lab = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate number of eggs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha2 <- 3.33; beta2 <- 0.05
c <- rlnorm(n, alpha2 + beta2 * s, 0.025)
plot(c ~ s, ylab = 'Culmen (mm)', xlab = 'Latent size', cex.lab = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate number of eggs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha3 <- 1.4; beta3 <- 0.15
# this is a very hacky way to use something akin to 'ordinal regression'
# to generate clutch sizes. The Conway-Maxwell-Poisson would probably
# be more appropriate, but this works...
e <- floor(exp(alpha3 + beta3 * s)) 
table(e)
plot(jitter(e) ~ s, ylab = 'Clutch size (eggs)', xlab = 'Latent size', cex.lab = 2)
plot(c ~ t, xlab = 'Tarsus length (mm)', ylab = 'Culmen length (mm)', cex.lab = 2)






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# note some pieces are missing!
# I didn't supply means for the intercepts... think about what appropriate
# priors might look like with a classmate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("build_a_SEM.jags")
cat("
    model {

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # intercepts need means
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alpha[1] ~ dnorm(, 1)
    alpha[2] ~ dnorm(, 1)
    alpha[3] ~ dnorm(, 1)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # regression parameters are ok
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    beta[1] = 1
    beta[2] ~ dnorm(0, 0.1)
    beta[3] ~ dnorm(0, 0.1)

    for (j in 1:3){
      sigma[j] ~ dgamma()
      tau[j] = 1/(sigma[j] * sigma[j])
    }


    for (i in 1:n){
    
      # modeling latent size as an unobserved, centered random variable
      s[i] ~ dnorm(0, tau[3])
      
      # modeling measured tarsus length and culment
      c[i] ~ dlnorm(alpha[1] + beta[1] * s[i], tau[1])
      t[i] ~ dlnorm(alpha[2] + beta[2] * s[i], tau[2])
      
      # modeling clutch size
      e[i] ~ dpois(exp(alpha[3] + beta[3] * s[i])) 

    }


    }
    ",fill = TRUE)
sink()

jags.data <- list(e = e, c = c, t = t, n = n)

inits <- function(){list()}  

parameters <- c('alpha','beta','sigma','s')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000


# Call JAGS from R 
# 50k iterations takes 3 seconds on an i9
library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "build_a_SEM.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

print(m)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework
#
# 1) Make a plot or two (or three) until you generally feel good about it.
#    If you don't, ask a friend, raise your hand or swing by office hours!
#    Suggestions:
#    plot expected clutch size against expected tarsus length
#    plot expected clutch size against size
#    plot expected culmen length against size
#    overlay the actual data points on top? Does the model fit 'well'?
#
# 2) change which beta in the model above is fixed. Observe how the parameter
#    estimates change, and compare them to your understanding from what we
#    just discussed in class.
#
# 3) Use blavaan below to generate SEM code and try to change the parameter
#    values and distributions to make it work more effectively 
#    (following the data generating process above)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





library(blavaan)
data = data.frame(e,t,c)
names(data) <- c('eggs','tarsus','culmen')
mb <- bsem(model = '
          size =~ tarsus + culmen
          eggs ~ 1 + size',
           data = dat,
           mcmcfile = T,
           target = 'jags')
summary(mb)
fitmeasures(mb)