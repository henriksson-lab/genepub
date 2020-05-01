library(greta)

# data
x <- as_data(iris$Petal.Length)
y <- as_data(iris$Sepal.Length)

# variables and priors
int <- normal(0, 1)
coef <- normal(0, 3)
sd <- student(3, 0, 1, truncation = c(0, Inf))

# operations
mean <- int + coef * x

# likelihood
distribution(y) <- normal(mean, sd)

# defining the model
m <- model(int, coef, sd)

# plotting
plot(m)

# sampling
draws <- mcmc(m, n_samples = 1000)


#################
#################
#################
#################


data(attitude)
design <- as.matrix(attitude[, 2:7])

int <- normal(0, 10)
coefs <- normal(0, 10, dim = ncol(design))
sd <- cauchy(0, 3, truncation = c(0, Inf))

# matrix multiplication is more efficient than multiplying the coefficients
# separately
mu <- int + design %*% coefs

distribution(attitude$rating) <- normal(mu, sd)





