library(reshape)
library(greta)
library(fda)


###############################################
## Function to calculate RMSE
calc_rmse <- function(residuals){
  RSS <- c(crossprod(residuals))
  MSE <- RSS / length(residuals)
  RMSE <- sqrt(MSE)
  RMSE
}



###############################################
## Replace NA in each column with a neutral value
replace_feat_na <- function(allfeat){
  for(i in 1:ncol(allfeat)){
    nalines <- is.na(allfeat[,i])
    allfeat[nalines,i] <- median(allfeat[!nalines,i])
  }
  allfeat
}



###############################################
## Read all data
totfeature <- read.csv("totfeature.csv", stringsAsFactors = FALSE)

cell_type <- "T cell"


###### Choose one cell type to work with
allfeat <- totfeature[totfeature$ct==cell_type,]
allfeat_meta <- allfeat[,c("gene","ct")]
allfeat <- allfeat[,!(colnames(allfeat) %in% c("gene","ct"))]
#colMeans(is.na(allfeat_red))  #low is good


###### Fill in NAs for now
allfeat <- replace_feat_na(allfeat)
allfeat <- as.data.frame(scale(allfeat))

###### TODO special treatment for the family index column. Cap it




#allfeat <- read.csv("for_gabija.csv", stringsAsFactors = FALSE)


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

###############################################
## Assemble all polynomials
## now mean(coefs) = mean(function) = 0    ... actually, this is not totally the same! 
use_features <- colnames(allfeat)[colnames(allfeat)!="rank_pmid"]
num_basepoints <- 5
design <- matrix(nrow=nrow(allfeat), ncol=0)
coefs <- normal(0, 10, dim = 0)
poly_bp <- list()
for(i in 1:length(use_features)){
  #Get basepoints and design
  pmat <- get_poly_mat(allfeat[,use_features[i]], num_basepoints = num_basepoints)
  design <- cbind(design, pmat$design)
  poly_bp[[use_features[i]]] <- pmat$bx
  
  #Collect coefficients
  coefsf <- normal(0, 10, dim = num_basepoints-1)
  coefs <- rbind(coefs, c(coefsf, -sum(coefsf)))
}





#########################################
### Add other variables
sd <- cauchy(0, 3, truncation = c(0, Inf))
intercept <- normal(0, 10)
mu <- intercept + design %*% coefs

distribution(allfeat$rank_pmid) <- normal(mu, sd)


#########################################
## Solve for the posterior. Store solution
m <- model(intercept, coefs, sd)
draws <- mcmc(m, n_samples = 3000, chains = 8)
s <- summary(draws)
s

write.csv(s$statistics, sprintf("greta/design/%s",design))
write.csv(s$statistics, sprintf("greta/fitted/%s",cell_type))
## store use_features?

## TODO Write an R object instead

#########################################
## Split up polynomial coefficients
poly_stat <- list()
poly_columns <- list()
colnames(design) <- rownames(s$statistics)
for(i in 1:length(use_features)){
  ## Get fitted statistics
  pstat <- s$statistics[grep("coefs",rownames(s$statistics)),]
  ## Get x-coordinates
  thei <- (1+(i-1)*num_basepoints):(i*num_basepoints)
  poly_stat[[use_features[i]]] <- pstat[thei,] 
  poly_columns[[use_features[i]]] <- rownames(pstat)[thei]
}


#########################################
## Calculate contributions of features to each gene
contribs <- matrix(nrow=nrow(allfeat), ncol=length(use_features))
for(i in 1:length(use_features)) {
  the_feature <- use_features[i]
  the_columns <- poly_columns[[the_feature]]
  contribs[,i] <- design[the_columns,poly] %*% s$statistics[the_columns,"Mean"]
}




#########################################
## Plot the function for one feature
plot_line_stat <- function(f, title=""){
  
  bp_x <- poly_bp[[f]]
  dist_coef <- poly_stat[[f]]
  
  ## Plot values to fit
  plot(allfeat[,f], allfeat[,"rank_pmid"], 
       pch=19, col="gray",cex=0.5, 
       xlab = f, ylab = "rank_pmid")
  
  ## Plot the polynomial
  lines(bp_x, dist_coef[,"Mean"],type="l")#, ylim = c(-1,1))#, title=title)
  lines(bp_x, dist_coef[,"Mean"] + dist_coef[,"SD"],type="l", col="red")
  lines(bp_x, dist_coef[,"Mean"] - dist_coef[,"SD"],type="l", col="red")
  
}


#plot_line_stat("rank_exp")
# plot_line_stat("perc_dependent_cells")

#########################################
## Plot all features functions
par(mfrow=c(3,3))
for(i in 1:length(use_features)){
  plot_line_stat(use_features[i])
}



##############################
#### MCMC quality control ####
##############################
# library (bayesplot)
# mcmc_trace(draws, facet_args = list(nrow = 6, ncol = 2))
# mcmc_intervals(draws)



