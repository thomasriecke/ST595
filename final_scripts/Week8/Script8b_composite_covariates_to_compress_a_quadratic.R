# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quadratic driver of multiple processes (composite covariate...)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# flow rate
# invertebrates (gamma) ~ flow + flow^2
# nest success (phi) ~ flow + flow^2
# chick survival ~ invertebrates
# recruitment (rho) ~ n * nest success * clutch size * chick survival
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(12345)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate flow rate across 40 years, and z-standardize it (f)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 40
flow <- rlnorm(n, 4, 0.25)
par(mfrow = c(1,1), mar = c(5.1,7.1,2.1,2.1))
plot(flow, type = 'b', xlab = 'Breeding season', cex = 1.5, cex.lab = 1.5, las = 1,
     ylab = expression(Flow~rate~(frac(m^3,s))))
f <- as.numeric(scale(flow))
# plot z-scored flow against actual measurements
plot(f ~ flow, cex.lab = 1.5, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parameters to simulate nest success, invertebrate density, and chick survival
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(4, 1.5, -2)
beta <- c(0, -0.35, 0, -0.35, 0.05)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate invertebrate abundance using log-normal regression
# and make some figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
invert <- rlnorm(n, alpha[1] + beta[1] * f + beta[2] * f^2, 0.1)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(invert ~ flow, ylab = 'Invertebrate abundance', cex = 1.5,
     xlab = 'Flow', cex.lab = 2, las = 1)
plot(invert, type = 'b', xlab = 'Breeding season', cex = 1.5, cex.lab = 1.5, las = 1,
     ylab = 'Invertebrate abundance')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate nest survival (psi) and make some figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psi <- plogis(alpha[2] + beta[3]*f + beta[4]*f^2 + rnorm(n,0,0.25))
plot(psi ~ flow, ylab = expression(Nest~survival~(psi)), ylim = c(0,1),
     cex = 1.5, las = 1, cex.lab = 1.5, las = 1, xlab = 'Flow')

plot(psi, ylab = expression(Nest~survival~(psi)), type = 'b', ylim = c(0,1),
     cex = 1.5, las = 1, cex.lab = 1.5, las = 1, xlab = 'Breeding season')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate a composite covariate and make some figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
delta <- NULL
delta[1] <- -0.1
delta[2] <- -1
composite <- gamma[1]*f + gamma[2]*f*f

plot(composite ~ flow, xlab = 'Flow', ylab = 'Composite', cex.lab = 2, las = 1)
plot(invert ~ composite, 
     ylab = 'Invertebrate abundance', xlab = 'Composite', cex.lab = 2, las = 1)
plot(invert ~ flow, 
     ylab = 'Invertebrate abundance', xlab = 'Flow', cex.lab = 2, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# chick survival (gamma)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gamma <- plogis(alpha[3] + beta[5]*invert + rnorm(n,0,0.25))
plot(gamma ~ invert, ylab = expression(Chick~survival~(gamma)), ylim = c(0,1),
     cex = 1.5, las = 1, cex.lab = 1.5, las = 1, xlab = 'Invertebrate abundance')
plot(gamma, ylab = expression(Chick~survival~(gamma)), type = 'b', ylim = c(0,1),
     cex = 1.5, las = 1, cex.lab = 1.5, las = 1, xlab = 'Breeding season')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate recruitment process
# r: potential Recruits
# h: Hatched nests
# zeta: clutch size (always 4)
# N: abundance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N <- NULL
h <- NULL
r <- NULL
N[1] <- 100

# juvenile and adult survival (I randomly chose numbers to make the population
# more or less stable in response to environmental fluctuation)
phi <- c(0.3125,0.5)

# mean clutch size
zeta <- 4

# first year
h <- rbinom(1, N[1], psi)
r <- rpois(1, h[1] * zeta * gamma[1])

for (t in 2:n){
  N[t] <- rpois(1, r[t-1] * phi[1] + N[t-1] * phi[2])
  h[t] <- rbinom(1, N[t], psi[t])
  r[t] <- rpois(1, h[t] * zeta * gamma[t])
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make figures of abundance, population growth rate (lambda)
# the number of recruits (r), and breeding pairs and hatched nests
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(N, type = 'b', cex = 1.5, ylab = 'Breeding pairs', 
     cex.lab = 1.5, las = 1, xlab = 'Breeding season')

lambda <- N[2:n]/N[1:(n-1)]
plot(lambda ~ flow[1:(n-1)], cex.lab = 1.5, las = 1, xlab = 'Flow',
     ylab = expression(Pop~growth~rate~(lambda)), cex = 1.5)

plot(r, type = 'b', cex = 1.5, ylab = '3-mo. old offspring', 
     cex.lab = 1.5, las = 1, xlab = 'Breeding season')


plot(N, ylab = 'Breeding pairs and hatched nests', ylim = c(0,500), type = 'b')
points(h, pch = 19, type = 'b')

plot(N, ylab = 'Breeding pairs and recruits', ylim = c(0,700), type = 'b')
points(r, pch = 19, type = 'b')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data we'll supply to JAGS are:
# 1) number of nests/breeding pairs per year (y)
# 2) number of hatched nests per year (h)
# 3) flow rate (f)
# 4) invertebrate abundance (i)
# 5) total number of fledged potential recruits (r)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("composite_flow.jags")
cat("
      model {

      delta[1] ~ dnorm(0, 1)
      delta[2] = -1

      alpha[1] ~ dnorm(4, 1)
      alpha[2] ~ dnorm(2, 1)
      alpha[3] ~ dnorm(-1, 1)
      
      beta[1] ~ dnorm(0.5, 1)
      beta[2] ~ dnorm(0.5, 1)
      beta[3] ~ dnorm(0.05, 1)
      
      sigma ~ dgamma(1,1)
      tau = 1/(sigma * sigma)


      for (t in 1:n){
      
        flow.effect[t] = delta[1]*f[t] + delta[2]*f[t]*f[t]
        i[t] ~ dlnorm(alpha[1] + beta[1]*flow.effect[t], tau[1])
        logit(psi[t]) = alpha[2] + beta[2]*flow.effect[t]
        
        h[t] ~ dbin(psi[t], y[t])
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # note we simulated gamma as 'chick survival'
        # here we model it as 'the number of female chicks produced per breeding season'
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        gamma[t] = exp(alpha[3] + beta[3]*i[t])
        r[t] ~ dpois(gamma[t] * h[t])
        
      }

      
      }
      ",fill = TRUE)
sink()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bundle data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jags.data <- list(y = N, n = n, i = invert, r = r, f = f, h = h)

# Initial values
inits <- function(){list(alpha = c(4,2,1))}  

# Parameters monitored
parameters <- c('alpha','beta','delta','sigma','gamma','flow.effect',
                'psi')

nc <- 4
nt <- 10
ni <- 25000
nb <- 10000

library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "composite_flow.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()
print(m, digits = 3)



plot(m$q50$flow.effect ~ flow)
plot(m$q50$psi ~ psi)
plot(m$q50$psi ~ flow)

plot(m$sims.list$delta[,1])


plot(m$sims.list$alpha[,1])
plot(m$sims.list$alpha[,2])
plot(m$sims.list$alpha[,3])

plot(m4$sims.list$beta[,1])
plot(m4$sims.list$beta[,2])
plot(m4$sims.list$beta[,3])
plot(m4$sims.list$deviance)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 1) Start by simulating a range of reasonable values for 'flow' 
#    (remember that we z-standardized this for use in the model) and picking
#    a number of initial breeding adults (perhaps 100)
#    
#    From there, simulate invertebrate abundance (i) and
#    the predicted number of hatched nests (h) given flow.
#    Then simulate the expected number of recruits given the number of hatched
#    nests and invertebrate abundance
#
# 2) Try to write the model without the composite (i.e., just put quadratic
#    effects of flow on invertebrate abundance and nest survival). What changes
#    from both a parameter estimate perspective and a prediction perspective?
#    Should both change, or just one of the two?
#
# 3) Explore the blavaan_issue.R script, examining the 'buggy' JAGS model
#    created using composite variable syntax (c <~ x1 + 1*x2)
#    what is the problem? (reference against the custom, correct model at the bottom)
# 
#    Also explore the revised model suggested by Ed Merkle (it works almost exactly
#    the same!). Why?
#
# 4) Go back to the top of this script. Put a trend on flow (make it decline over
#    time). Simulate the population. What happens? Try this a few times... heck
#    make the river dry up?
#    
#    Now go back and add additional variance over time, have the river either flood
#    or be very low? What does that do?! What do you predict it will do...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



