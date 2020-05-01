library(reshape)
library(greta)

allfeat <- read.csv("for_gabija.csv", stringsAsFactors = FALSE)

###################################################################
#### Construct polynomial design matrix
get_poly_mat <- function(px, num_basepoints = 5){
  #px <- allfeat$rank_exp
  bpr <- range(px)
  bpr[1] <- bpr[1] - 0.001
  bpr[2] <- bpr[2] + 0.001
  
  ## Decide base points ... return?
  bp_x <- bpr[1] + (bpr[2]-bpr[1])*((1:num_basepoints)-1)/(num_basepoints-1)
  bp_dx <- bp_x[2]-bp_x[1]
  
  ## Which two base points are the closest?
  ind <- floor((px-bp_x[1])/bp_dx)+1
  ## Calculate weights
  w1 <- 1-(px-bp_x[ind])/bp_dx
  w2 <- 1- w1
  
  ## Generate values in sparse form
  npoint <- nrow(allfeat)
  pdata <- rbind(
    data.frame(
      row=1:npoint,
      col=ind,
      val=w1
    ),
    data.frame(
      row=1:npoint,
      col=ind+1,
      val=w2
    )
  )
  
  pmat <- cast(pdata, row~col, fill = 0, value="val")[,-1]
  list(bx = bp_x,
       design = pmat)
}


pmat <- get_poly_mat(allfeat$rank_exp)

design <- pmat$design

### All polynomials in one go (faster one at a time? - but harder to construct?)


#### Set up polynomial coefficients such that mean(coefs) = mean(function) = 0    ... actually, this is not totally the same(?)
coefsf <- normal(0, 10, dim = ncol(design)-1)
coefs <- c(coefsf, -sum(coefsf))

### Assembling it

sd <- cauchy(0, 3, truncation = c(0, Inf))
intercept <- normal(0, 10)
mu <- intercept + design %*% coefs
#mu <- intercept + design[,1:4] %*% coefsf -   design[,5] * (-sum(coefsf))

distribution(allfeat$rank_pmid) <- normal(mu, sd)


###################### Solving it all 
m <- model(intercept, coefs, sd)

#plot(m)

# sampling
draws <- mcmc(m, n_samples = 3000, chains = 8)
summary(draws)
s <- summary(draws)


plot_line_stat <- function(dist_coef){
  plot(pmat$bx, dist_coef[,"Mean"],type="l")
  lines(pmat$bx, dist_coef[,"Mean"] + dist_coef[,"SD"],type="l", col="red")
  lines(pmat$bx, dist_coef[,"Mean"] - dist_coef[,"SD"],type="l", col="red")
}
plot_line_stat(s$statistics[grep("coefs",rownames(s$quantiles)),])

#mean_coef <- s$quantiles[grep("coefs",rownames(s$quantiles)),"50%"]


##############################
##############################
##############################
##############################
##############################

## intersect has been removed already. functions should have a mean=0  ---  how do we do this here?



library (bayesplot)
mcmc_trace(draws, facet_args = list(nrow = 6, ncol = 2))
mcmc_intervals(draws)



