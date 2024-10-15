# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quadratic driver of multiple processes
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
# simulate flow rate, and z-standardize it (f)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 40
flow <- rlnorm(n, 4, 0.25)
par(mfrow = c(1,1), mar = c(5.1,6.1,2.1,2.1))
plot(flow, type = 'b',
     ylab = expression(Flow~rate~(frac(m^3,s))))
f <- as.numeric(scale(flow))
plot(flow ~ f)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parameters to simulate nest success and invertebrate density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(4, 1.5, -2)
beta <- c(0, -0.35, 0, -0.35, 0.05)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate invertebrate abundance using log-normal regression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
invert <- rlnorm(n, alpha[1] + beta[1] * f + beta[2] * f^2, 0.1)
plot(invert ~ flow)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate nest survival (psi) using beta regression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psi <- plogis(alpha[2] + beta[3]*f + beta[4]*f^2 + rnorm(n,0,0.25))
plot(psi ~ flow)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# chick survival (gamma)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gamma <- plogis(alpha[3] + beta[5]*invert + rnorm(n,0,0.25))
plot(gamma ~ invert)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate recruitment process
# r: potential Recruits
# h: Hatched nests
# zeta: clutch size
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N <- NULL
h <- NULL
r <- NULL
N[1] <- 100

# juvenile and adult survival
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

plot(N, type = 'b')
plot(r ~ flow)
plot(c(r/N) ~ flow)
plot(c(r/N) ~ invert)
plot(c(r/N) ~ psi)

cor(invert,psi)

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

      gamma[1] ~ dnorm(0, 1)
      gamma[2] = -1

      alpha[1] ~ dnorm(4, 1)
      alpha[2] ~ dnorm(2, 1)
      alpha[3] ~ dnorm(-1, 1)
      
      beta[1] ~ dnorm(0.5, 1)
      beta[2] ~ dnorm(0.5, 1)
      beta[3] ~ dnorm(0.05, 1)
      
      sigma[1] ~ dgamma(1,1)
      tau[1] = 1/(sigma[1] * sigma[1])


      for (t in 1:n){
      
        flow[t] = gamma[1]*f[t] + gamma[2]*f[t]*f[t]
        i[t] ~ dlnorm(alpha[1] + beta[1]*flow[t], tau[1])
        logit(psi[t]) = alpha[2] + beta[2]*flow[t]
        
        h[t] ~ dbin(psi[t], y[t])
        r[t] ~ dpois(exp(alpha[3] + beta[3]*i[t]) * h[t])
        
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
parameters <- c('alpha','beta','sigma','gamma','flow',
                'psi')

nc <- 4
nt <- 10
ni <- 25000
nb <- 10000

library(jagsUI)
Sys.time()
m4 <- jags(jags.data, inits, parameters, "composite_flow.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()
print(m4, digits = 3)



plot(m4$q50$flow ~ flow)
plot(m4$q50$psi ~ psi)
plot(m4$q50$psi ~ flow)


plot(m4$sims.list$gamma[,1])


plot(m4$sims.list$alpha[,1])
plot(m4$sims.list$alpha[,2])
plot(m4$sims.list$alpha[,3])

plot(m4$sims.list$beta[,1])
plot(m4$sims.list$beta[,2])
plot(m4$sims.list$beta[,3])
plot(m4$sims.list$deviance)






