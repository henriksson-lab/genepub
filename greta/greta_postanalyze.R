library(reshape)
library(greta)
library(splines)
library(ggplot2)
library(plotly)
library(directlabels)



###############################################
## Function to calculate RMSE
calc_rmse <- function(residuals){
  RSS <- c(crossprod(residuals))
  MSE <- RSS / length(residuals)
  RMSE <- sqrt(MSE)
  RMSE
}


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
all_poly_fits <- NULL
all_rmse <- NULL
for(i in 1:length(all_cell_type)) {
  cell_type <- all_cell_type[i]
  #print(cell_type)
  
  #########################################
  ## Read data  
  the_fit <- readRDS(sprintf("greta/out/%s",cell_type))
  use_features <- names(the_fit$all_poly)
  
  all_rmse <- rbind(
    all_rmse,
    data.frame(ct=cell_type, rmse=calc_rmse_on_fit(the_fit))
  )
  
  
  # 
  #   #########################################
  #   ## Plot all feature-functions
  #   par(mfrow=c(3,3))
  #   for(i in 1:length(use_features)){
  #     plot_line_stat(the_fit, use_features[i])
  #   }
  # 
  #   
  features <- names(the_fit$all_poly)
  for(the_feature in features){
    one_poly <- the_fit$all_poly[[the_feature]]
    
    x <- one_poly$plot_px
    x <- x-min(x)
    x <- x/max(x)
    
    y <- one_poly$plot_pmat %*% one_poly$stat[,"Mean"]

    all_poly_fits <- rbind(
      all_poly_fits,
      data.frame(
        ct=cell_type,
        feature=the_feature,
        x=x,
        y=y
      )
    )
  }
  
}

# 
# #############################################################
# ## Function: Run a function for data for each ct
# apply_per_groupcol <- function(dat, col, groupcol, func){
#   #for(col in colnames(dat)){
#   for(ct in unique(dat$groupcol)){
#     dat[dat$groupcol==ct,col] <- func(dat[dat$groupcol==ct,col])
#   }
#   #}
#   dat
# }
# 
# normalize01 <- function(dat){
#   dat <- dat-min(dat)
#   dat <- dat/max(dat)
#   dat
# }

plot_heatmap_for_feature <- function(feature) {
  dat <- all_poly_fits[all_poly_fits$feature==feature,]
  all_cell_type <- list.files("greta/out")
  xout <- (0:29)/29
  hmap <- matrix(0, nrow=length(all_cell_type), ncol=length(xout))
  for(i in 1:length(all_cell_type)) {
    cell_type <- all_cell_type[i]
    hmap[i, ] <- approx(dat[dat$ct==cell_type,"x"],dat[dat$ct==cell_type,"y"], xout=xout)$y
  }
  rownames(hmap) <- all_cell_type
  #hmap <- hmap[,2:30]
  heatmap(hmap,Colv = NA,xlab = feature)  
}

unique(all_poly_fits$feature)
plot_heatmap_for_feature("perc_dependent_cells")
plot_heatmap_for_feature("nearby_pmid")
plot_heatmap_for_feature("ppi")
plot_heatmap_for_feature("coexp10")
plot_heatmap_for_feature("family_indexdiff")
plot_heatmap_for_feature("rank_exp")




plot_lines_for_feature <- function(feature) {
  dat <- all_poly_fits[all_poly_fits$feature==feature,]
  ggplot(data=dat, aes(x=x, y=y, group=ct)) +
    geom_line() +
    geom_dl(aes(label = ct), method = list(dl.combine("first.points", "last.points"), cex = 0.8)) 
}
plot_lines_for_feature("rank_exp")

#plot_ly(data = dat[sample(1:nrow(dat), 1000),], x = ~x, y = ~y, line_group=~ct,  text = ~paste("g: ",ct))


# # Change the line type
# ggplot(data=df, aes(x=dose, y=len, group=1)) +
#   geom_line(linetype = "dashed")+
#   geom_point()
# 
# 
# 
# 
# 
# 
# 
# 
# 
