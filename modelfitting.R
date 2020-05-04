library(reshape2)
library(gridExtra)


lm_rmse <- function(res){
  RSS <- c(crossprod(res$residuals))
  MSE <- RSS / length(res$residuals)
  RMSE <- sqrt(MSE)
  RMSE
}


###############################################
## Nicely plot a correlation matrix
ggplot_cor <- function(thecor){
  thecor <- round(thecor,digits = 2)
  for(i in 1:nrow(thecor)){
    thecor[i,i] <- 0
  }
  ggp <- ggplot(melt(thecor), aes(Var1, Var2, fill=value)) + 
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
    geom_bar(position="stack", stat="identity") +
    geom_hline(yintercept = 0) +
    geom_text(aes(label = condition), position = position_stack(vjust = 0.5)) +
    theme(legend.position = "none")
}




#############################################################
############## Prepare data #################################
#############################################################

list.files("greta/feature/")
cell_type <- "T cell"
# cell_type <- "B cell"
#cell_type <- "fibroblast"
#cell_type <- "epithelial cell"

allfeat <- read.csv(sprintf("greta/feature/%s.csv",cell_type))
allfeat_meta <- allfeat[,c("gene","ct","orig_year")]
allfeat_red <- allfeat[,!(colnames(allfeat) %in% c("gene","ct","orig_year"))]

#############################################################
############## Simple global linear model ###################
#############################################################


## Try to explain with a linear model



#thecor <- cor(allfeat_red)
thecor <- cor(allfeat_red[,!(colnames(allfeat_red) %in% c("pmid_count","rank_pmid"))])
ggsave("plots/lm_corr.pdf", ggplot_cor(thecor))

################ Simplest linear model attempt

thelm <- lm(
  rank_pmid ~ rank_exp + coexp10 + perc_dependent_cells + nearby_pmid + ppi + first_year + family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas + homology_pmid, 
  allfeat_red)
round(thelm$coefficients, digits = 5)
lm_rmse(thelm)
plot_staple_coef_multi(list(foo=thelm))
ggsave("plots/out_lm_all.pdf",plot = plot_staple_coef_multi(list(foo=thelm)))


# 
# (Intercept)             rank_exp              coexp10 perc_dependent_cells          nearby_pmid                  ppi 
# 0.00000              0.08114              0.04286              0.06232              0.11542              0.02327 
# first_year     family_indexdiff  family_founder_rank          rank_cosmic            rank_gwas 
# -0.59176             -0.08730              0.16013              0.02263             -0.02824 





#############################################################
############## Linear model over time #######################
#############################################################

fit_lm <- function(allfeat_red){
  #Should definitely not rescale - the scale should be representative. any shifts in points goes into the intercept
  thelm <- lm(
    rank_pmid ~ rank_exp + coexp10 + perc_dependent_cells + nearby_pmid + ppi + first_year + family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas + homology_pmid, 
    allfeat_red)
  thelm
}


############# Early genes vs late genes
thelm_early <- fit_lm(allfeat_red[allfeat_meta$orig_year>=1970 & allfeat_meta$orig_year<2000,])
thelm_late  <- fit_lm(allfeat_red[allfeat_meta$orig_year>=2000 & allfeat_meta$orig_year<2020,])
plot_staple_coef_multi(list(
  early=thelm_early,
  late=thelm_late
))


#### Windowed fit
fit_year <- NULL
coef_year <- NULL
for(i in 1970:2000){
  f <- fit_lm(allfeat_red[allfeat_red_year>=i & allfeat_red_year<i+20,])
  fit_year <- rbind(
    fit_year,
    data.frame(year=i, rmse=lm_rmse(f))
    )
  coef_year <- rbind(coef_year, f$coefficients[-1])
}

#### Rescale coefficients
coef_year <- as.data.frame(coef_year)
for(i in 1:nrow(coef_year)){
#  coef_year[i,] <- coef_year[i,]/sum(abs(coef_year[i,]))
#    coef_year[i,] <- coef_year[i,]/sum(coef_year[i,])
#  coef_year[i,] <- coef_year[i,]/sd(coef_year[i,])
  coef_year[i,] <- coef_year[i,]/sd(coef_year[i,-7])
}

#### RMSE over time
plot(fit_year$year, fit_year$rmse, type="l")

#### Plot every feature over time
list_plots <- list()
for(the_feature in colnames(coef_year)){
  
  dat <- data.frame(
    coef=coef_year[,the_feature],
    year=fit_year$year)
  
  list_plots[[the_feature]] <- ggplot(dat, aes(y=coef,x=year)) +
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
    thecor <- cor(allfeat_red[,-(1:2)])
    #round(thecor,digits = 2)
    sum_cor <- sum_cor + thecor
    num_cor <- num_cor + 1

    ### Fit model
    thelm <- lm(
      rank_pmid ~ rank_exp + coexp10 + perc_dependent_cells + nearby_pmid + ppi + first_year + family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas, 
      allfeat_red)
    #round(thelm$coefficients, digits = 5)
    #lm_rmse(thelm)
    #plot_staple_coef_multi(list(foo=thelm))
    
    all_ct_lm_coef <- rbind(
      all_ct_lm_coef,
      data.frame(
        ct=cell_type,
        t(as.data.frame(thelm$coefficients[-1])))
    )
    
  } else {
    print("   skipping, NA data")
  }
}
rownames(all_ct_lm_coef) <- NULL

#########################
## Look at average correlation
ggsave("plots/lm_corr_avgct.pdf", ggplot_cor(sum_cor/num_cor))
ggplot_cor(sum_cor/num_cor)


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


