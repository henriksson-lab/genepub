library(reshape2)



lm_rmse <- function(res){
  #res<- thelm$residuals
  RSS <- c(crossprod(res$residuals))
  MSE <- RSS / length(res$residuals)
  RMSE <- sqrt(MSE)
  RMSE
}






allfeat <- read.csv("totfeature.csv", stringsAsFactors = FALSE)
head(allfeat)

colMeans(is.na(allfeat[,-(1:2)]))

###############################################
## Replace NA in each column with a neutral value
replace_feat_na <- function(allfeat){
  for(i in 1:ncol(allfeat)){
    nalines <- is.na(allfeat[,i])
    allfeat[nalines,i] <- median(allfeat[!nalines,i])
  }
  allfeat
}


###### Choose one cell type to work with
allfeat_red <- allfeat[allfeat$ct=="B cell",]
allfeat_red_meta <- allfeat_red[,c("gene","ct")]
allfeat_red <- allfeat_red[,!(colnames(allfeat_red) %in% c("gene","ct"))]
colMeans(is.na(allfeat_red))  #low is good

###### Fill in NAs for now
allfeat_red <- replace_feat_na(allfeat_red)
allfeat_red <- as.data.frame(scale(allfeat_red))

write.csv(allfeat_red, "for_gabija.csv", row.names = FALSE)

## Try to explain with a linear model
thecor <- cor(allfeat_red)
round(thecor,digits = 2)
#thecor <- cor(allfeat_red)[1:6,1:6]
# for(i in 1:nrow(thecor)) {
#   thecor[i,i] <- 0
# }
# qplot(x=Var1, y=Var2, data=melt(thecor), fill=value, geom="tile") + 
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours = c("red","blue","white"))

plot(allfeat_red$perc_dependent_cells, allfeat_red$rank_exp)




thelm <- lm(
  rank_pmid ~ rank_exp + coexp10 + perc_dependent_cells + nearby_pmid + ppi + family_index + family_founder_rank_pmid, 
  allfeat_red)
round(thelm$coefficients, digits = 5)

lm_rmse(thelm)

#compcoef <- c("rank_exp","percentage_of_dep_cells","nearby_pmid")






######colnames(feature_genefam) <- c("gene","family_index","rank_founder_pmid")
#"gene","family_index","family_founder_rank_pmid"
feature_genefam$family_founder_rank_pmid <- rank(feature_genefam$family_founder_rank_pmid)
feature_genefam$new_founder_rank_pmid <- 0.8^(feature_genefam$family_index-1)*feature_genefam$family_founder_rank_pmid
feature_genefam <- feature_genefam[,c("gene","new_founder_rank_pmid")]
print(length(unique(feature_genefam$gene)))



################################################
############# Test linear model ################
################################################

## Try to explain with a linear model
dat <- dat[dat$ct=="T cell",]
thelm <- lm(rank_pmid~rank_exp + percentage_of_dep_cells + nearby_pmid,dat)
compcoef <- c("rank_exp","percentage_of_dep_cells","nearby_pmid")

### full model
# thelm <- lm(rank_pmid~rank_exp + percentage_of_dep_cells + nearby_pmid + family_founder_rank_pmid + family_index,dat)
# compcoef <- c("rank_exp","percentage_of_dep_cells","nearby_pmid", "family_founder_rank_pmid","family_index")
thelm <- lm(rank_pmid~rank_exp + percentage_of_dep_cells + nearby_pmid + new_founder_rank_pmid,dat)
compcoef <- c("rank_exp","percentage_of_dep_cells","nearby_pmid", "new_founder_rank_pmid")

## Check the contributions
datred <- dat[,compcoef]
for(i in compcoef){
  datred[,i] <- datred[,i]*thelm$coefficients[i]/dat$rank_pmid
}
apply(datred,2,mean)

# 1] "gene"                     "ct"                      
# [3] "rank_exp"                 "percentage_of_dep_cells" 
# [5] "rank_pmid"                "family_index"            
# [7] "family_founder_rank_pmid" "nearby_pmid" 



