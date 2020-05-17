library(reshape2)
library(gridExtra)


lm_rmse <- function(res){
  RSS <- c(crossprod(res$residuals))
  MSE <- RSS / length(res$residuals)
  RMSE <- sqrt(MSE)
  RMSE
}





get_lm_weights <- function(thelm){
  return(thelm$coefficients[-1]) ####################
  
  if(FALSE){
    b <- thelm$coefficients[-1]
    outv <- allfeat_red[,names(b)]
    for(i in 1:ncol(outv)){
      outv[,i] <- outv[,i]*b[i]/allfeat_red$rank_pmid       #dividing here makes this model rather unstable. low rank genes cause issues
    }
    cm <- colMeans(outv)
    cm
  } else {
    b <- thelm$coefficients[-1]
    outv <- allfeat_red[,names(b)]
    for(i in 1:ncol(outv)){
      outv[,i] <- outv[,i]*b[i]
    }
    cm <- colMeans(outv)
    cm
  }
}



get_lm_relweights <- function(thelm){
  cm <- get_lm_weights(thelm)
#  cm <- cm/sd(cm)
  cm <- cm/sum(pmax(cm,0))
  #cm <- cm/-cm["first_year"]
  #cm <- cm[names(cm)!="first_year"]
  #/-min(colMeans(outv))
  cm
}

###############################################
## Function: Model fitting
fit_lm <- function(allfeat_red){
  
  print(sprintf("Rows to fit on: %s",nrow(allfeat_red)))
  # thelm <- lm(
  #   rank_pmid ~ rank_exp + coexp10 +   ppi + 
  #     nearby_pmid +
  #     essentiality_global +# essentiality_ct + 
  #     first_year + #first_year2 + first_year3 + # first_year4+
  #     #rank_gwas_1 +
  #     family_indexdiff + 
  #     family_founder_rank + #confounded with index diff
  #     rank_cosmic + 
  #     rank_gwas + homology_pmid, 
  #   allfeat_red)
  
  
  thelm <- lm(
    rank_pmid ~ rank_exp + coexp10 +   ppi + 
      nearby_pmid +
      essentiality_global +
      first_year + 
      family_indexdiff + 
      family_founder_rank + 
      rank_cosmic + 
      rank_gwas + homology_pmid, 
    allfeat_red)
  
  
  if("founder_fitted_pmid" %in% colnames(allfeat_red)){
    thelm <- lm(
      rank_pmid ~ rank_exp + coexp10 +   ppi + 
        nearby_pmid +
        essentiality_global +
        first_year + 
        founder_fitted_pmid +
        family_indexdiff +
        rank_cosmic + 
        rank_gwas + homology_pmid, 
      allfeat_red)
    
  }
  
  
  #Should definitely not rescale - the scale should be representative. any shifts in points goes into the intercept
  # thelm <- lm(
  #   rank_pmid ~ rank_exp + coexp10 + essentiality_global + nearby_pmid + ppi + first_year + 
  #     family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas + homology_pmid, 
  #   allfeat_red)
  thelm
  
  #essentiality_ct  or   essentiality_global
}




###############################################
## Function: Model fitting
fit_lm_nofounder <- function(allfeat_red){
  
  thelm <- lm(
    rank_pmid ~ rank_exp + coexp10 +   ppi + 
      nearby_pmid +
      essentiality_global +
      first_year + 
      #family_indexdiff + 
      #family_founder_rank + 
      rank_cosmic + 
      rank_gwas + homology_pmid, 
    allfeat_red)
  
  thelm
}


###############################################
## Nicely plot a correlation matrix
ggplot_cor <- function(thecor){
  thecor <- round(thecor,digits = 2)
  for(i in 1:nrow(thecor)){
    thecor[i,i] <- 0
  }
  tm <- melt(thecor)
  colnames(tm) <- c("Var1","Var2","value")
  ggp <- ggplot(tm, aes(Var1, Var2, fill=value)) + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="red", high="blue") +
    geom_tile() +
    #  geom_text(aes(label = round(value, digits = 2))) 
    geom_text(aes(label = value)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  plot(ggp)
  ggp
}

###############################################
## Function: Plot linear model coefficients for multiple models, as staple diagram
plot_staple_coef_multi <- function(list_thelm){
  thelm_coef <- NULL
  for(n in names(list_thelm)){
    thelm <- list_thelm[[n]]
    print(n)
    #print(thelm)
    w <- get_lm_relweights(thelm)
    
    one_coef <- data.frame(
      n=n,
      sgn=sign(w),
      coef=w,
      condition=names(w)
    )
    #one_coef$coef <- one_coef$coef / sum(one_coef$coef[one_coef$coef>0])
    
    thelm_coef <- rbind(
      thelm_coef,
      one_coef
    )
  }
  
  ggplot(thelm_coef, aes(fill=condition, y=coef,x=n)) + 
    geom_bar(position="stack", stat="identity") +
    geom_hline(yintercept = 0) +
    geom_text(aes(label = condition), position = position_stack(vjust = 0.5)) +
    theme(legend.position = "none")
}






#############################################################
############## Simple global linear model ###################
#############################################################


list.files("greta/feature/")
cell_type <- "T cell"
#cell_type <- "B cell"
#cell_type <- "fibroblast"
#cell_type <- "epithelial cell"
allfeat <- read.csv(sprintf("greta/feature/%s.csv",cell_type))

### before 1975, not so much. definitely avoid <1950
allfeat <- allfeat[allfeat$orig_year>=1970,]  ################# note, affects before after!


meta_names <- c("gene","ct","founder_name",   "orig_year", "has_founder")
allfeat_meta <- allfeat[,meta_names]
allfeat_red <- allfeat[,!(colnames(allfeat) %in% meta_names)]
allfeat_meta$has_founder[is.na(allfeat_meta$has_founder)] <- FALSE

#hist(allfeat$rank_pmid)

if(FALSE){
  plot(
    rank(allfeat$rank_gwas),
    allfeat$rank_pmid,
    cex=0.5)

  plot(
    allfeat$rank_exp,
    allfeat$rank_pmid,
    cex=0.5)
  
  plot(
    allfeat$essentiality_global,   #Clear trend
    allfeat$rank_pmid,
    cex=0.5)
  
  plot(
    allfeat$essentiality_ct,   #this seems broken?  or just because few cell types?
    allfeat$rank_pmid,
    cex=0.5)  
  
  plot(
    allfeat$nearby_pmid,   #looks fine - just weak trend
    allfeat$rank_pmid,
    cex=0.5)  

  plot(
    allfeat$rank_cosmic,   #clear trend
    allfeat$rank_pmid,
    cex=0.5)  
  
  plot(
    (allfeat$homology_pmid),   #trend is more clear if using rank. but there in either case
    allfeat$rank_pmid,
    cex=0.5)  
  

  plot(
    (allfeat$family_founder_rank),  #super clear ... but could be an age effect
    allfeat$rank_pmid,
    cex=0.5)  

  plot(
    (allfeat$family_indexdiff)+runif(nrow(allfeat),0,0.2),   #can't say that much
    allfeat$rank_pmid,
    cex=0.5)  
  

  plot(
    rank(allfeat$coexp10),   #yup, a trend. better if using rank .. maybe
    allfeat$rank_pmid,
    cex=0.5)  
  
  plot(
    rank(allfeat$ppi),   #meah.
    allfeat$rank_pmid,
    cex=0.5)  
  
}


if(FALSE){
  #thecor <- cor(allfeat_red)
  thecor <- cor(allfeat_red[,!(colnames(allfeat_red) %in% c("pmid_count","rank_pmid"))], method = "spearman")
  ggplot_cor(thecor)
  #ggsave("plots/lm_corr.pdf", ggplot_cor(thecor))
  
  
  thecor <- cor(allfeat_red[,!(colnames(allfeat_red) %in% c("pmid_count"))], method = "spearman")
  ggplot_cor(thecor)
  
}


################ Simplest linear model attempt

thelm  <- fit_lm(allfeat_red)
round(thelm$coefficients, digits = 5)
anova(thelm)


lm_rmse(thelm)
plot_staple_coef_multi(list(foo=thelm))
#ggsave("plots/out_lm_all.pdf",plot = plot_staple_coef_multi(list(foo=thelm)))

### Since data not centered, can now rescale


get_lm_relweights(thelm)
get_lm_weights(thelm)


### Plot time vs weight used
# 
# plot(
#   allfeat_meta$orig_year,
#   allfeat_red$first_year*thelm$coefficients["first_year"]
# )
# 
# plot(
#   allfeat_meta$orig_year,
#   allfeat_red$rank_exp*thelm$coefficients["rank_exp"]
# )
# 
# plot(
#   allfeat_meta$orig_year,
#   allfeat_red$essentiality_global*thelm$coefficients["essentiality_global"]
# )
# 
# plot(
#   allfeat_meta$orig_year,
#   allfeat_red$rank_cosmic*thelm$coefficients["rank_cosmic"]
# )
# 
# plot(
#   allfeat_meta$orig_year,
#   allfeat_red$family_founder_rank*thelm$coefficients["family_founder_rank"]
# )
# 
# plot(
#   allfeat_meta$orig_year,
#   allfeat_red$founder_fitted_pmid#*thelm$coefficients["founder_fitted_pmid"]
# )

#density(allfeat_meta$orig_year, allfeat_red$family_founder_rank*thelm$coefficients["family_founder_rank"])

allfeat_redmeta <- allfeat_red #$orig_year
allfeat_redmeta$orig_year <- allfeat_meta$orig_year
ggplot(allfeat_redmeta, aes(orig_year,founder_fitted_pmid)) + geom_point() + geom_smooth()
ggplot(allfeat_redmeta, aes(orig_year,rank_cosmic)) + geom_point() + geom_smooth()
ggplot(allfeat_redmeta, aes(orig_year,essentiality_global)) + geom_point() + geom_smooth()


if(FALSE){
  
  
  b <- thelm$coefficients[-1]
  outv <- allfeat_red[,names(b)]
  for(i in 1:ncol(outv)){
    outv[,i] <- outv[,i]*b[i]#/allfeat_red$rank_pmid       #dividing here makes this model rather unstable. low rank genes cause issues
  }
  
  i<-1
  colnames(outv)[i]
  
  va <- data.frame(
    val=outv[,i],
    year=allfeat_meta$orig_year)
  va <- sqldf("select avg(val) as val, year from va group by year order by year")
  plot(allfeat_meta$orig_year + runif(nrow(allfeat_meta)), outv[,i], ylab=colnames(outv)[i], xlab="year")
  lines(va$year, va$val, col="red")
  
  plot(
    rank_nogap(allfeat_red$rank_exp),
    rank_nogap(allfeat_red$rank_pmid), cex=0.5)
  #rank_pmid always more than cosmic - cosmic starts later
  #not the case with gwas
  #why with exp?  whyyyyyyyy?
  
  plot(
    rank_nogap(allfeat_red$rank_exp),
    rank_nogap(allfeat_red$rank_pmid - outv$first_year), cex=0.5)
  
  #d <- density(x, bw = 5, 
               
               # "nrd0", adjust = 1,
               # kernel = c("gaussian", "epanechnikov", "rectangular",
               #            "triangular", "biweight",
               #            "cosine", "optcosine"),
               # weights = NULL, window = kernel, width,
               # give.Rkern = FALSE,
               # n = 512, from, to, cut = 3, na.rm = FALSE, ...))
  
  #cm <- colMeans(outv)
  #cm
  
}

# for(i in 1:nrow(outv)){
#   outv[i,] <- outv[i,]/allfeat_red$rank_pmid[i]
# }
#as.matrix(allfeat_red[,names(b)]) %*% t(t(b))

# 
# (Intercept)             rank_exp              coexp10 perc_dependent_cells          nearby_pmid                  ppi 
# 0.00000              0.08114              0.04286              0.06232              0.11542              0.02327 
# first_year     family_indexdiff  family_founder_rank          rank_cosmic            rank_gwas 
# -0.59176             -0.08730              0.16013              0.02263             -0.02824 

#############################################################
############## Fitted founder ###############################
#############################################################


##If we do this then it need to be done on all features to be fair. better not get into this!
if(FALSE){

  
  ## First fit everything
  thelm_nof  <- fit_lm_nofounder(allfeat_red)
  
  
  
  ## Take predicted pmid, plug into those that have it as founder
  # founder_fit <- data.frame(
  #   founder_name=allfeat_meta$gene,
  #   fitted_pmid=thelm_nof$fitted.values
  # )
  # founder_fit <- sqldf("select founder_name, avg(fitted_pmid) as fitted_pmid from founder_fit group by founder_name")
  # rownames(founder_fit) <- founder_fit$founder_name
  
  ## Predicted pmid too aggresive. should just compensate by removing all other factors from it
  founder_fit <- data.frame(
    founder_name=allfeat_meta$gene,
    fitted_pmid=allfeat_red$rank_pmid - thelm_nof$fitted.values
  )
  founder_fit <- sqldf("select founder_name, avg(fitted_pmid) as fitted_pmid from founder_fit group by founder_name")
  rownames(founder_fit) <- founder_fit$founder_name
  
  
  ## Create fitted founder
  allfeat_red$founder_fitted_pmid <- founder_fit[allfeat$gene,]$fitted_pmid
  #mean(is.na(allfeat$founder_fitted_pmid))
  allfeat_meta$is_founder <- allfeat$family_indexdiff==min(allfeat$family_indexdiff)   #hack
  allfeat_red$founder_fitted_pmid[allfeat_red$is_founder] <- NA
  allfeat_red$founder_fitted_pmid[!allfeat_meta$has_founder] <- NA
  
  
  ##Fill in the NAs for those without a founder
  allfeat_red$founder_fitted_pmid[is.na(allfeat_red$founder_fitted_pmid)] <- mean(allfeat_red$founder_fitted_pmid, na.rm=TRUE)
  
  
  thelm_fittedf <- lm(
    rank_pmid ~ rank_exp + coexp10 +   ppi + 
      nearby_pmid +
      essentiality_global +
      first_year + 
      founder_fitted_pmid +
      family_indexdiff +
      rank_cosmic + 
      rank_gwas + homology_pmid, 
    allfeat_red)
  
  plot_staple_coef_multi(list(
    nof=thelm_nof,
    ff=thelm_fittedf))
  #plot_staple_coef_multi(list(foo=thelm_nof))
  

    
}



#############################################################
############## Linear model over time #######################
#############################################################


############# Early genes vs late genes
thelm_early <- fit_lm(allfeat_red[allfeat_meta$orig_year>=1970 & allfeat_meta$orig_year<1996,])
thelm_late  <- fit_lm(allfeat_red[allfeat_meta$orig_year>=1996 & allfeat_meta$orig_year<2020,])
plot_staple_coef_multi(list(
  early=thelm_early,
  late=thelm_late
))


#### Windowed fit
fit_year <- NULL
coef_year <- NULL
n_paper <- NULL
for(i in 1970:2000){
  allfeat_red_w <- allfeat_red[allfeat_meta$orig_year>=i & allfeat_meta$orig_year<i+20& allfeat_meta$orig_year<last(allfeat_meta$orig_year),]
  f <- fit_lm(allfeat_red_w)
  fit_year <- rbind(
    fit_year,
    data.frame(year=i, rmse=lm_rmse(f))
    )
  n_paper <- c(n_paper, nrow(allfeat_red_w))
  coef_year <- rbind(coef_year, f$coefficients[-1])
}

#### Rescale coefficients
coef_year <- as.data.frame(coef_year)
for(i in 1:nrow(coef_year)){
#  coef_year[i,] <- coef_year[i,]/sum(abs(coef_year[i,]))
#    coef_year[i,] <- coef_year[i,]/sum(coef_year[i,])
#  coef_year[i,] <- coef_year[i,]/sd(coef_year[i,])
#  coef_year[i,] <- coef_year[i,]/sd(coef_year[i,-7])   #what is 7?
}

#### RMSE over time
plot(fit_year$year, fit_year$rmse, type="l")
plot(fit_year$year, n_paper, type="l")   ### beep. 1970-1980, not many!

#### Plot every feature over time
list_plots <- list()
for(the_feature in colnames(coef_year)){
  
  dat <- data.frame(
    coef=coef_year[,the_feature],
    year=fit_year$year)
  
  list_plots[[the_feature]] <- ggplot(dat, aes(y=coef,x=year)) +  ################## just not possible to do it this way! #######################
    geom_line(size=1, colour="blue") +
    ylim(min(c(dat$coef,0)), max(c(dat$coef,0))) +
    xlab("") +
    ylab(the_feature) +
    scale_x_continuous(breaks=c(1970,2000)) +
    theme(
      #panel.background = element_rect(fill = "white"),
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      #axis.text.y=element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
}
ggp <- grid.arrange(grobs=list_plots, nrow=2)
ggsave("plots/lm_feature_over_time.pdf",ggp)


############### read.csv("website/data/feature_long_name.csv")   extend this one! and map names above




#############################################################
###### Simple global linear model, all cell types ###########
#############################################################

all_cell_types <- str_split_fixed(list.files("greta/feature/"),".csv",2)[,1]
all_ct_lm_coef <- NULL
num_cor <- 0
sum_cor <- 0
for(cell_type in all_cell_types){
  #cell_type <- "T cell"
  # cell_type <- "B cell"
  #cell_type <- "fibroblast"
  #cell_type <- "epithelial cell"
  print(cell_type)
  
  allfeat <- read.csv(sprintf("greta/feature/%s.csv",cell_type))
  allfeat_meta <- allfeat[,c("gene","ct","orig_year")]
  allfeat_red <- allfeat[,!(colnames(allfeat) %in% c("gene","ct","orig_year"))]
  
  if(!any(is.na(allfeat_red)) & nrow(allfeat_red>3000)) {

    ## Try to explain with a linear model
    thecor <- cor(allfeat_red[,-(1)])
    #round(thecor,digits = 2)
    sum_cor <- sum_cor + thecor
    num_cor <- num_cor + 1

    ### Fit model
    thelm <- fit_lm(allfeat_red)
    # thelm <- lm(
    #   rank_pmid ~ rank_exp + coexp10 + perc_dependent_cells + nearby_pmid + ppi + first_year + family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas, 
    #   allfeat_red)
    #round(thelm$coefficients, digits = 5)
    #lm_rmse(thelm)
    #plot_staple_coef_multi(list(foo=thelm))
    
    all_ct_lm_coef <- rbind(
      all_ct_lm_coef,
      data.frame(
        ct=cell_type,
        t(get_lm_relweights(thelm)))
    )
    
  } else {
    print("   skipping, NA data")
  }
}
#rownames(all_ct_lm_coef) <- NULL

#########################
## Look at average correlation
ggplot_cor(sum_cor/num_cor)
#ggsave("plots/lm_corr_avgct.pdf", ggplot_cor(sum_cor/num_cor))


#########################
## Plot it all on one page
list_plots <- list()
for(the_feature in setdiff(colnames(all_ct_lm_coef),c("ct"))){
  
  dat <- data.frame(
    coef=all_ct_lm_coef[,the_feature],
    ct=all_ct_lm_coef$ct)
  
  list_plots[[the_feature]] <- ggplot(dat, aes(y=coef,x=ct)) +
    geom_bar(position="stack", stat="identity") + 
    coord_flip() +
    ylab(the_feature)
}
ggp <- grid.arrange(grobs=list_plots, nrow=2)
plot(ggp)
ggsave("plots/lm_feature_vs_ct_gg.pdf",ggp,width = 25, height = 10)
ggp


# 
# 
# 
# melt(all_ct_lm_coef[,c("ct","rank_exp")])
# #ggplot(melt(all_ct_lm_coef[,c("ct","rank_exp")]), aes(fill=ct, y=value)) 
# 
# ggplot(melt(all_ct_lm_coef[,c("ct","rank_exp")]), aes(fill=ct, y=value,x=ct)) +
#   geom_bar(position="stack", stat="identity")
  
# geom_bar(position="stack", stat="identity") +
  # geom_hline(yintercept = 0) +
  # geom_text(aes(label = condition), position = position_stack(vjust = 0.5)) +
  # theme(legend.position = "none")


#ggsave("plots/out_lm_all.pdf",plot = plot_staple_coef_multi(list(foo=thelm)))






















######colnames(feature_genefam) <- c("gene","family_index","rank_founder_pmid")
#"gene","family_index","family_founder_rank_pmid"
feature_genefam$family_founder_rank_pmid <- rank(feature_genefam$family_founder_rank_pmid)
feature_genefam$new_founder_rank_pmid <- 0.8^(feature_genefam$family_index-1)*feature_genefam$family_founder_rank_pmid
feature_genefam <- feature_genefam[,c("gene","new_founder_rank_pmid")]
print(length(unique(feature_genefam$gene)))


