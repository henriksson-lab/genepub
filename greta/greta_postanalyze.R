library(reshape)
library(greta)
library(splines)



#########################################
## Function: Calculate RMSE
calc_rmse_on_fit <- function(the_fit){
  pstat <- the_fit$statistics[grep("coefs",rownames(the_fit$statistics)),]
  res <- the_fit$allfeat$rank_pmid - the_fit$statistics["intercept","Mean"] - the_fit$design %*% pstat[,"Mean"]
  calc_rmse(res)
}



#########################################
## Function: Plot the function for one feature
plot_line_stat <- function(the_fit, f, xlab=f){
  
  all_poly <- the_fit$all_poly
  allfeat <- the_fit$allfeat
  
  dist_coef <- all_poly[[f]]$stat
  
  ## Plot values to fit
  plot(allfeat[,f], allfeat[,"rank_pmid"], 
       pch=19, col="gray",cex=0.5, 
       xlab = xlab, ylab = "rank_pmid")
  
  ## Plot the polynomial
  lines(
    all_poly[[f]]$plot_px,
    all_poly[[f]]$plot_pmat %*% dist_coef[,"Mean"],
    col="red")
  ## Approx error bounds using sd. Need otherwise bootstrap over all polynomials, messy
  lines(
    all_poly[[f]]$plot_px,
    all_poly[[f]]$plot_pmat %*% (dist_coef[,"Mean"] + dist_coef[,"SD"]),
    col="red")  
  lines(
    all_poly[[f]]$plot_px,
    all_poly[[f]]$plot_pmat %*% (dist_coef[,"Mean"] - dist_coef[,"SD"]),
    col="red")  
  
}





#########################################################
## Process all cell types
all_cell_type <- list.files("greta/out")
for(i in all_cell_type) {
  cell_type <- all_cell_type[i]

  #########################################
  ## Read data  
  the_fit <- readRDS(sprintf("greta/out/%s",cell_type))
  use_features <- names(the_fit$all_poly)
  calc_rmse_on_fit(the_fit)
  

  #########################################
  ## Plot all feature-functions
  par(mfrow=c(3,3))
  for(i in 1:length(use_features)){
    plot_line_stat(the_fit, use_features[i])
  }
  
  
  ##TODO show all the functions across cell types in one diagram.
  ## this might take some normalization on the y-axis ... or maybe that is already sorted actually
  
}







