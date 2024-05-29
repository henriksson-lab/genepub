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
 
  thelm <- lm(
    rank_pmid ~ rank_exp + coexp10 +   ppi + 
      nearby_pmid +
      essentiality_global +
      first_year + 
      family_indexdiff + 
      family_founder_rank + 
      rank_cosmic + 
      #rank_tad + 
      rank_gwas + homology_pmid, 
    allfeat_red)

  thelm
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
#allfeat <- allfeat[allfeat$orig_year>=1970,]  ################# note, affects before after!

allfeat <- allfeat[allfeat$orig_year>=1950,]  ################# note, affects before after!


meta_names <- c("gene","ct","founder_name",   "orig_year", "has_founder")
allfeat_meta <- allfeat[,meta_names]
allfeat_red <- allfeat[,!(colnames(allfeat) %in% meta_names)]
allfeat_meta$has_founder[is.na(allfeat_meta$has_founder)] <- FALSE


#allfeat_red <- allfeat_red[allfeat_red$rank_exp > quantile(allfeat_red$rank_exp,probs = 0.95),]

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
  
  # plot(
  #   rank(allfeat$ppi),   #meah.
  #   allfeat$rank_pmid,
  #   cex=0.5)  
  
}


if(FALSE){
  #thecor <- cor(allfeat_red)
  thecor <- cor(allfeat_red[,!(colnames(allfeat_red) %in% c("pmid_count","rank_pmid"))], method = "spearman")
  ggplot_cor(thecor)
  #ggsave("plots/lm_corr.pdf", ggplot_cor(thecor))

  thecor <- cor(allfeat_red[,!(colnames(allfeat_red) %in% c("pmid_count","rank_gwas_1","rank_tad","essentiality_ct","rank_pmid"))], method = "spearman")
  ggplot_cor(thecor)
  ggsave("plots/out_lm_all_corr.pdf", width = 8, height = 7)
  
}


################ Simplest linear model attempt

thelm  <- fit_lm(allfeat_red)
round(thelm$coefficients, digits = 5)
anova(thelm)

lm_rmse(thelm)
plot_staple_coef_multi(list(foo=thelm))
#ggsave("plots/out_lm_all.pdf",plot = plot_staple_coef_multi(list(foo=thelm)))

(thelm$coefficients["rank_gwas"]+thelm$coefficients["essentiality_global"]+thelm$coefficients["rank_cosmic"])/sum(abs(thelm$coefficients))
(thelm$coefficients["rank_exp"] + thelm$coefficients["coexp10"]+thelm$coefficients["ppi"])/sum(abs(thelm$coefficients))
(abs(thelm$coefficients["family_indexdiff"])+thelm$coefficients["homology_pmid"]+thelm$coefficients["family_founder_rank"]+thelm$coefficients["nearby_pmid"])/sum(abs(thelm$coefficients))
(thelm$coefficients["first_year"])/sum(abs(thelm$coefficients))
45+10+20+25

#(thelm$coefficients["coexp10"]+thelm$coefficients["ppi"])/sum(abs(thelm$coefficients))



#(thelm$coefficients["rank_expcoexp10"]+thelm$coefficients["ppi"])/sum(abs(thelm$coefficients))

### Since data not centered, can now rescale
# get_lm_relweights(thelm)
# get_lm_weights(thelm)




#############################################################
############## Linear model over time #######################
#############################################################


# thelm_early <- fit_lm(allfeat_red[allfeat_meta$orig_year>=1970 & allfeat_meta$orig_year<1996,])
# thelm_late  <- fit_lm(allfeat_red[allfeat_meta$orig_year>=1996 & allfeat_meta$orig_year<=2018,])

############# Early genes vs late genes
thelm_early <- fit_lm(as.data.frame(scale(allfeat_red[allfeat_meta$orig_year>=1970 & allfeat_meta$orig_year<1990,])))
thelm_late  <- fit_lm(as.data.frame(scale(allfeat_red[allfeat_meta$orig_year>=1991 & allfeat_meta$orig_year<=2010,])))
# thelm_early <- fit_lm(allfeat_red[allfeat_meta$orig_year>=1970 & allfeat_meta$orig_year<1990,])
# thelm_late  <- fit_lm(allfeat_red[allfeat_meta$orig_year>=1991 & allfeat_meta$orig_year<=2010,])
plot_staple_coef_multi(list(
  early=thelm_early,
  late=thelm_late
))
ggsave("plots/linmod_early_vs_late.pdf")


w <- cbind(
  allfeat_meta[allfeat_meta$orig_year>=1991,],
  allfeat_red[allfeat_meta$orig_year>=1991,])
w <- w[order(w$rank_pmid, decreasing = TRUE),]
w


#### Windowed fit
fit_year <- NULL
coef_year <- NULL
n_paper <- NULL
for(i in 1970:1995){
  allfeat_red_w <- allfeat_red[allfeat_meta$orig_year>=i & allfeat_meta$orig_year<i+25,]
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
  
  coef_year[i,] <- coef_year[i,]/sum(pmax(coef_year[i,],0))
  
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
######### Focus on genes expressed in a cell type ###########
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

### Only consider more expressed genes; need to rescale the variable after this!
allfeat_red <- allfeat_red[allfeat_red$rank_exp > quantile(allfeat_red$rank_exp,probs = 0.75),]
#allfeat_red <- allfeat_red[allfeat_red$rank_exp != min(allfeat_red$rank_exp),]
allfeat_red$rank_exp <- scale(allfeat_red$rank_exp)


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





#############################################################
######### Store residuals for website #######################
#############################################################

all_cell_types <- unique(totfeature$ct)
all_res_coord <- NULL
for(cell_type in all_cell_types) {
  print(cell_type)
  ##### Fit the model
  allfeat <- read.csv(sprintf("greta/feature/%s.csv",cell_type))
  allfeat <- allfeat[allfeat$orig_year>=1950,]  ################# note, affects before after!
  meta_names <- c("gene","ct","founder_name",   "orig_year", "has_founder")
  allfeat_meta <- allfeat[,meta_names]
  allfeat_red <- allfeat[,!(colnames(allfeat) %in% meta_names)]
  allfeat_meta$has_founder[is.na(allfeat_meta$has_founder)] <- FALSE
  
  thelm  <- fit_lm(allfeat_red)
  
  ##### Extract residuals
  dat <- data.frame(gene=allfeat_meta$gene)
  dat[,sprintf("y_%s", celltype)] <- thelm$residuals
  dat <- dat[order(dat[,sprintf("y_%s", celltype)]),]
  dat[,sprintf("x_%s", celltype)] <- 1:nrow(dat)
  
  if(is.null(all_res_coord)){
    all_res_coord <- dat
  } else {
    all_res_coord <- merge(all_res_coord, dat, all=TRUE)
  }
}
store_website_coordinates("residual",all_res_coord)







###############################################################################################################
########## Produce plot: NN fit ###############################################################################
###############################################################################################################

dat <- read.csv("plots/SimpleNN_out_weight.csv")[,-1]

thelm_coef <- data.frame(
  n="NN",
  coef=dat$val,
  condition=dat$feature
)
#one_coef$coef <- one_coef$coef / sum(one_coef$coef[one_coef$coef>0])


ggplot(thelm_coef, aes(fill=condition, y=coef,x=n)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = condition), position = position_stack(vjust = 0.5)) +
  theme(legend.position = "none")
ggsave("plots/SimpleNN_out_weight.pdf")





###############################################################################################################
########## Produce plot: XG fit ###############################################################################
###############################################################################################################

dat <- read.csv("plots/SimpleXG_out_weight.csv")[,-1]

thelm_coef <- data.frame(
  n="NN",
  coef=dat$val,
  condition=dat$feature
)
#one_coef$coef <- one_coef$coef / sum(one_coef$coef[one_coef$coef>0])


ggplot(thelm_coef, aes(fill=condition, y=coef,x=n)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = condition), position = position_stack(vjust = 0.5)) +
  theme(legend.position = "none")
ggsave("plots/SimpleXG_out_weight.pdf")





