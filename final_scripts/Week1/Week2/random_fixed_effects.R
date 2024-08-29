# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# future needs: 
#
# use a red fox (or other organism) home range shape file to 
# generate realistic points and demonstrate how to use spatially-structured
# random effects, Bergmann's and Allen's rules
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In our first example, we'll simulate a certain number of populations or
# locations (l) to sample individuals from. 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
l <- 10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we'll imagine that there is some 'perfectly average' value for the body size
# of an adult (female?) fox, let's say 4 kg
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mu.star <- 4

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's imagine that there is some variation in size among each population!
# in the future, this will be based on ecological concepts, but for now it's
# UTTERLY RANDOM... we will call the average adult size for each population
# mu, and we will call the variance among populations sigma.star^2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigma.star <- 0.25
mu <- rnorm(l, mu.star, sigma.star)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ok, to rehash
# we just simulated some number (l = 20) of imaginary populations
# we simulated a 'grand mean' of female fox body size (mu* = 4)
# we simulated some variation among populations (sigma.star = 0.25 kg)
# and then randomly sampled 20 means for each population (mu)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hist(mu, breaks = 15, main = NULL, xlab = expression(Mean~mass~(mu*';'~kg)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's create some 'real' foxes to catch and measure! 
# There are two (really lots more) ways we
# could do this:
#
# 1) simulate the entirety of each population and then subsample each population
# 2) determine a number of samples and then simulate data for those samples
# 
# it shouldn't matter? and option 2 is often easier, so we'll do that
# 'n' will be the number of foxes we're going to catch at each site,
# so it will be a vector of 'l' length
# lambda will be the average number of foxes we catch at each site,
# for now let's set that to 20
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lambda <- 20
n <- rpois(l, lambda)

# total number of foxes we caught
sum(n)

# now let's create a vector (s) that tells us the site where each fox was caught
s <- rep(1:l, times = n)
s

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in my simulation, there were 21 foxes at site 1 (i.e., n[1] = 21)
# the first 21 values in s are 1
# similarly, there were 21 foxes sampled at site 2, the next 21 values in s
# are 2, i.e., s[22:42] = 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's simulate the mass of each fox, since these are our data
# we'll call them 'y'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first create an empty vector
y <- matrix(NA, sum(n))

# define a 'within population' variance (sigma^2)
sigma <- 0.15

# we could do this 1) in a for-loop, 2) in a vector
# for both, we'll used 'nested indexing'
for (i in 1:sum(n)){
  y[i] = rnorm(1, mu[s[i]], sigma)
}

# this is equivalent
y = rnorm(sum(n), mu[s], sigma)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a data .csv
# note, you can do the same thing for model output! we'll revisit this idea soon,
# for now just remember, please don't ever copy values out of R by hand!
# it introduces errors and wastes your valuable time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path <- 'C:\\Users\\thomas.riecke\\Desktop\\SEM_TVR\\output\\Week2\\fox.data.csv'
df <- data.frame(cbind(y,s))
names(df) <- c('mass','site')
write.csv(df, path)
rm(df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make some figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1))

# distribution of site-level means
vioplot::vioplot(mu, names = '', ylab = 'Mean mass of each population', 
                 col = 'coral2', side = 'right', ylim = c(min(y),max(y)))
points(jitter(rep(0.9,l)), mu, col = 'black', cex = 1.5)

# all masses
vioplot::vioplot(y, names = '',ylab = 'Fox mass (kg)', 
                 col = 'coral2', side = 'right', ylim = c(min(y),max(y)))
points(y ~ jitter(rep(0.9, sum(n))), col = 'black', cex = 1)

# distributions of foxes in each population
vioplot::vioplot(y ~ s, xlab = 'Population', ylab = 'Fox mass (kg)', 
                 col = 'coral2', side = 'right', ylim = c(min(y),max(y)))
points(y ~ jitter(c(s-0.2)), col = 'black', cex = 0.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ok, it's time to write a JAGS model (or two)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









