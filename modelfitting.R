library(reshape2)
library(gridExtra)


lm_rmse <- function(res){
  RSS <- c(crossprod(res$residuals))
  MSE <- RSS / length(res$residuals)
  RMSE <- sqrt(MSE)
  RMSE
}





get_lm_weights <- function(thelm){
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
  #cm <- cm/sd(cm)
  cm <- cm/-cm["first_year"]
  cm <- cm[names(cm)!="first_year"]
  #/-min(colMeans(outv))
  cm
}

###############################################
## Function: Model fitting
fit_lm <- function(allfeat_red){
  
  allfeat_red$first_year2 <- allfeat_red$first_year**2
  allfeat_red$first_year3 <- allfeat_red$first_year**3
  allfeat_red$first_year4 <- allfeat_red$first_year**4
  
  thelm <- lm(
    rank_pmid ~ rank_exp + coexp10 + essentiality_global + nearby_pmid + ppi + 
      first_year + first_year2 + first_year3 + first_year4+
      family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas + homology_pmid, 
    allfeat_red)
  
  #Should definitely not rescale - the scale should be representative. any shifts in points goes into the intercept
  # thelm <- lm(
  #   rank_pmid ~ rank_exp + coexp10 + essentiality_global + nearby_pmid + ppi + first_year + 
  #     family_indexdiff + family_founder_rank + rank_cosmic + rank_gwas + homology_pmid, 
  #   allfeat_red)
  thelm
  
  #essentiality_ct  or   essentiality_global
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
#cell_type <- "T cell"
# cell_type <- "B cell"
#cell_type <- "fibroblast"
cell_type <- "epithelial cell"
allfeat <- read.csv(sprintf("greta/feature/%s.csv",cell_type))
allfeat_meta <- allfeat[,c("gene","ct","orig_year")]
allfeat_red <- allfeat[,!(colnames(allfeat) %in% c("gene","ct","orig_year"))]


if(FALSE){
  #thecor <- cor(allfeat_red)
  thecor <- cor(allfeat_red[,!(colnames(allfeat_red) %in% c("pmid_count","rank_pmid"))])
  ggplot_cor(thecor)
  #ggsave("plots/lm_corr.pdf", ggplot_cor(thecor))
  
  
  thecor <- cor(allfeat_red[,!(colnames(allfeat_red) %in% c("pmid_count"))])
  ggplot_cor(thecor)
  
}


################ Simplest linear model attempt

thelm  <- fit_lm(allfeat_red)
round(thelm$coefficients, digits = 5)


lm_rmse(thelm)
plot_staple_coef_multi(list(foo=thelm))
#ggsave("plots/out_lm_all.pdf",plot = plot_staple_coef_multi(list(foo=thelm)))

### Since data not centered, can now rescale


get_lm_relweights(thelm)
get_lm_weights(thelm)

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
############## Linear model over time #######################
#############################################################


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
n_paper <- NULL
for(i in 1970:2000){
  allfeat_red_w <- allfeat_red[allfeat_meta$orig_year>=i & allfeat_meta$orig_year<i+20,]
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
  coef_year[i,] <- coef_year[i,]/sd(coef_year[i,-7])   #what is 7?
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


