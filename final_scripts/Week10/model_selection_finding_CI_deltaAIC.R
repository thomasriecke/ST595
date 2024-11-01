library(latex2exp)
library(MuMIn)

n = c(25,50,100)
sims <- 5000

aic1 <- matrix(NA, sims, length(n))
aic2 <- matrix(NA, sims, length(n))
delta <- matrix(NA, sims, length(n))
beta1x <- matrix(NA, sims, length(n))
beta2x <- matrix(NA, sims, length(n))
beta2m <- matrix(NA, sims, length(n))
p1x <- matrix(NA, sims, length(n))
p2x <- matrix(NA, sims, length(n))
p2m <- matrix(NA, sims, length(n))

# this will run 15000 simulations (5000 each for different sample sizes [n])
# each simulation will generate a covariate (x)
# a response (y), and a meaningless covariate
for (ii in 1:sims){
  for (jj in 1:3){ # use ii and jj to help avoid accidentally overwriting inside loop

    # simulate covariate (x), response (y), and meaningless covariate (m)
    x <- rnorm(n[jj], 0, 1)
    y <- rnorm(n[jj], x, 1)
    m <- rnorm(n[jj], 0, 1)
    
    m1 <- glm(y ~ x)
    m2 <- glm(y ~ x + m)
    
    aic1[ii,jj] <- AICc(m1)
    aic2[ii,jj] <- AICc(m2)
    
    delta[ii,jj] <- aic1[ii,jj] - aic2[ii,jj]
    
    beta1x[ii,jj] <- summary(m1)$coefficients[2,1]
    beta2x[ii,jj] <- summary(m2)$coefficients[2,1]
    beta2m[ii,jj] <- summary(m2)$coefficients[3,1]    

    p1x[ii,jj] <- summary(m1)$coefficients[2,4]
    p2x[ii,jj] <- summary(m2)$coefficients[2,4]
    p2m[ii,jj] <- summary(m2)$coefficients[3,4]     
    
  }
}

plot(delta[,1] ~ p2m[,1], col = 'red',
     ylab = TeX("$\\delta$AICc"), 
     xlab = TeX("p-value for $\\beta_2$"),
     las = 1, cex.lab = 2, xlim = c(0,0.5), cex = 0.25)
points(delta[,2] ~ p2m[,2], col = 'purple', cex = 0.25)
points(delta[,3] ~ p2m[,3], col = 'blue', cex = 0.25)
abline(h = 0)
# these are 'best guesses', not exact points
abline(v = 0.15, col = 'blue')
abline(v = 0.13, col = 'purple')
abline(v = 0.115, col = 'red')
# note that they'll change with different model likelihoods/error distributions
# key point here is that the alpha value for significance for an additional 
# covariate that leads to an improvement in AIC will change with changing 
# sample sizes

# proportion of simulations where 'm' leads to an 
# improved model even though it's random/meaningless
length(which(delta[,1] > 0))/sims
length(which(delta[,2] > 0))/sims
length(which(delta[,3] > 0))/sims
