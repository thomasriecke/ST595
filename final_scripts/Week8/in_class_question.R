library(lavaan)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define sample size and generate collinear covariates (x[,1:3])
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 150
x <- matrix(NA, n, 3)
x[,1] <- rnorm(n, 0, 1)
x[,2] <- x[,1] + rnorm(n, 0, 1)
x[,3] <- x[,1] + rnorm(n, 0, 1)
cor(x)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define slopes (intercept = 0)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta <- c(1,0.5,2)

y <- rnorm(n, beta[1] * x[,1] + beta[2] * x[,2] + beta[3] * x[,3], 0.25)
plot(y ~ x[,1])

composite = x %*% beta
plot(y ~ composite)

data <- data.frame(y,x)
summary(m1 <- glm(y ~ x))

library(lavaan)
summary(m2 <- sem('y ~ c
             c <~ X1 + 1 * X2 + X3',
            data = data))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note the relationship between the composite loadings (under composites below)
# and the beta values? 
# Are the relative effects the same? 
# Absolute? (effect of c on y under regressions)
# What happens if we fix the loading of X1 on the composite to 1 on line 28?
# Why?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
