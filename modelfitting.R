library(reshape2)

#############################################################
############## Prepare data #################################
#############################################################

lm_rmse <- function(res){
  RSS <- c(crossprod(res$residuals))
  MSE <- RSS / length(res$residuals)
  RMSE <- sqrt(MSE)
  RMSE
}


allfeat <- read.csv("totfeature.csv", stringsAsFactors = FALSE)

##QC
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
allfeat_red <- allfeat[allfeat$ct=="T cell",]
allfeat_red_meta <- allfeat_red[,c("gene","ct")]
allfeat_red <- allfeat_red[,!(colnames(allfeat_red) %in% c("gene","ct"))]
colMeans(is.na(allfeat_red))  #low is good

###### Fill in NAs for now
allfeat_red <- replace_feat_na(allfeat_red)


#apply(allfeat_red,2,min)
#allfeat_red <- log10(allfeat_red)   #Implicitly on log-scale already

### Rescale; keep the year before scaling
allfeat_red_year <- allfeat_red$first_year
allfeat_red <- as.data.frame(scale(allfeat_red))

#############################################################
############## Simple global linear model ###################
#############################################################


## Try to explain with a linear model
thecor <- cor(allfeat_red)
round(thecor,digits = 2)



################ Simplest linear model attempt

thelm <- lm(
  rank_pmid ~ rank_exp + coexp10 + perc_dependent_cells + nearby_pmid + ppi + first_year + family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas, 
  allfeat_red)
round(thelm$coefficients, digits = 5)
lm_rmse(thelm)

# 
# (Intercept)             rank_exp              coexp10 perc_dependent_cells          nearby_pmid                  ppi 
# 0.00000              0.08114              0.04286              0.06232              0.11542              0.02327 
# first_year     family_indexdiff  family_founder_rank          rank_cosmic            rank_gwas 
# -0.59176             -0.08730              0.16013              0.02263             -0.02824 

#############################################################
############## Multiple linear models #######################
#############################################################


# 
# plot_staple_coef <- function(thelm){
#   thelm_coef <- data.frame(
#     sgn=sign(thelm$coefficients[-1]),
#     coef=thelm$coefficients[-1],
#     condition=names(thelm$coefficients)[-1]
#   )
#   thelm_coef$coef <- thelm_coef$coef/sum(thelm_coef$coef[thelm_coef$coef>0])
#   ggplot(thelm_coef, aes(fill=condition, y=coef,x=1)) + 
#     geom_bar(position="stack", stat="identity")
# }

###############################################
## Function: Plot linear model coefficients for multiple models, as staple diagram
plot_staple_coef_multi <- function(list_thelm){
  thelm_coef <- NULL
  for(n in names(list_thelm)){
    thelm <- list_thelm[[n]]
    print(n)
    one_coef <- data.frame(
      n=n,
      sgn=sign(thelm$coefficients[-1]),
      coef=thelm$coefficients[-1],
      condition=names(thelm$coefficients)[-1]
    )
    one_coef$coef <- one_coef$coef / sum(one_coef$coef[one_coef$coef>0])
    
    
    thelm_coef <- rbind(
      thelm_coef,
      one_coef
    )
  }
  
  ggplot(thelm_coef, aes(fill=condition, y=coef,x=n)) + 
    geom_bar(position="stack", stat="identity")
}


plot_staple_coef(thelm)

#compcoef <- c("rank_exp","percentage_of_dep_cells","nearby_pmid")




############# Early genes vs late genes

fit_lm <- function(allfeat_red){
  #Should definitely not rescale - the scale should be representative. any shifts in points goes into the intercept
  thelm <- lm(
    rank_pmid ~ rank_exp + coexp10 + perc_dependent_cells + nearby_pmid + ppi + first_year + family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas, 
    allfeat_red)
  thelm
}

# thelm_early <- fit_lm(allfeat_red[allfeat_red_year>=1980 & allfeat_red_year<1995,])
# thelm_late  <- fit_lm(allfeat_red[allfeat_red_year>=1995 & allfeat_red_year<2010,])
thelm_early <- fit_lm(allfeat_red[allfeat_red_year>=1980 & allfeat_red_year<2000,])
thelm_late  <- fit_lm(allfeat_red[allfeat_red_year>=2000 & allfeat_red_year<2020,])
plot_staple_coef_multi(list(
  early=thelm_early,
  late=thelm_late
))


fit_year <- NULL
coef_year <- NULL
for(i in 1970:2000){
  f <- fit_lm(allfeat_red[allfeat_red_year>=i & allfeat_red_year<i+20,])
  fit_year <- rbind(
    fit_year,
    data.frame(year=i, rmse=lm_rmse(f))
    )
  coef_year <- rbind(coef_year, f$coefficients)
}
plot(fit_year$year, fit_year$rmse, type="l")

coef_year <- as.data.frame(coef_year)
for(i in 1:nrow(coef_year)){
#  coef_year[i,] <- coef_year[i,]/sum(abs(coef_year[i,]))
#    coef_year[i,] <- coef_year[i,]/sum(coef_year[i,])
#  coef_year[i,] <- coef_year[i,]/sd(coef_year[i,])
  coef_year[i,] <- coef_year[i,]/sd(coef_year[i,-7])
}

#coef_year$pe
plot(fit_year$year, coef_year$first_year, type="l")
plot(fit_year$year, coef_year$perc_dependent_cells, type="l")
plot(fit_year$year, coef_year$rank_exp, type="l")
plot(fit_year$year, coef_year$ppi, type="l")
plot(fit_year$year, coef_year$coexp10, type="l")
plot(fit_year$year, coef_year$nearby_pmid, type="l")
plot(fit_year$year, coef_year$family_founder_rank, type="l")
#plot(fit_year$year, coef_year$, type="l")


hist(allfeat_red_year)





























######colnames(feature_genefam) <- c("gene","family_index","rank_founder_pmid")
#"gene","family_index","family_founder_rank_pmid"
feature_genefam$family_founder_rank_pmid <- rank(feature_genefam$family_founder_rank_pmid)
feature_genefam$new_founder_rank_pmid <- 0.8^(feature_genefam$family_index-1)*feature_genefam$family_founder_rank_pmid
feature_genefam <- feature_genefam[,c("gene","new_founder_rank_pmid")]
print(length(unique(feature_genefam$gene)))


