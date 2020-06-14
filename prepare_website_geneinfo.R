dat <- read.csv("input/nih/Mus_musculus.gene_info", sep="\t", stringsAsFactors = FALSE)
dat <- dat[,c("GeneID","Symbol","description")]

dat2 <- read.csv("input/nih/gene2ensembl.mouse", sep="\t", header=FALSE, stringsAsFactors = FALSE)
dat2 <- dat2[,c(2,3)]
colnames(dat2) <- c("GeneID","Ensembl_gene_identifier")

dat3 <- read.csv("features/feature_minyear_gene_ct.csv", stringsAsFactors = FALSE)
colnames(dat3) <- c("Symbol","firstyear")
#feature_minyear_gene_ct

dat4 <- read.csv("features/feature_ranked_pmid.csv", stringsAsFactors = FALSE)
dat4$rank_pmid <- 10^(dat4$rank_pmid)-1
colnames(dat4) <- c("Symbol","numcitations")

  
############################################################
######################### XGboost ##########################
############################################################

all_dat5 <- NULL
all_xgboost_residual <- NULL
for(fname in list.files("plots/xgboost_prediction/")){
  cell_type <- str_split_fixed(fname,pattern = "\\.",2)[1]
  print(cell_type)

  dat5 <- read.csv(sprintf("plots/xgboost_prediction/%s",fname),stringsAsFactors = FALSE)
  colnames(dat5) <- c("Symbol","xgscore")
  
  ### Unscale predicted #citations
  ### Note: final_score = (rank_pmid - mean(allfeat$rank_pmid)) / sd(allfeat$rank_pmid)
  allfeat <- totfeature[totfeature$ct=="T cell",]
  the_sd <- sd(allfeat$rank_pmid)
  the_mean <- mean(allfeat$rank_pmid)
  dat5$xgscore <- dat5$xgscore * the_sd + the_mean
  dat5$xgscore <- 10^dat5$xgscore - 1

  ### Generate coordinate system  
  temp <- merge(dat5, dat4)
  xgboost_residual <- data.frame(
    gene=temp$Symbol
  )
  xgboost_residual[,sprintf("y_%s",cell_type)] <- log10(temp$xgscore+1) - log10(temp$numcitations+1)
  xgboost_residual <- xgboost_residual[order(xgboost_residual[,sprintf("y_%s",cell_type)]),]
  xgboost_residual[,sprintf("x_%s",cell_type)] <- 1:nrow(xgboost_residual)
  
  if(is.null(all_xgboost_residual)){
    all_xgboost_residual <- xgboost_residual
  } else {
    all_xgboost_residual <- merge(all_xgboost_residual, xgboost_residual)
  }
  
  #store_website_coordinates("residual_xg",xgboost_residual)
  dat5$ct <- cell_type
  all_dat5 <- rbind(
    all_dat5, dat5
  )
}
store_website_coordinates("residual_xg",all_xgboost_residual)

# dat5 <- read.csv("plots/xgboost.csv",stringsAsFactors = FALSE)
# colnames(dat5) <- c("Symbol","xgscore")
# 
# ### Scale back predicted #citations
# #final_score = (rank_pmid - mean(allfeat$rank_pmid)) / sd(allfeat$rank_pmid)
# #final_score * sd(allfeat$rank_pmid) + mean(allfeat$rank_pmid) = rank_pmid 
# allfeat <- totfeature[totfeature$ct=="T cell",]
# the_sd <- sd(allfeat$rank_pmid)
# the_mean <- mean(allfeat$rank_pmid)
# 
# dat5$xgscore <- dat5$xgscore * the_sd + the_mean
# dat5$xgscore <- 10^dat5$xgscore - 1
# 
# temp <- merge(dat5, dat4)
# 
# xgboost_residual <- data.frame(
#   gene=temp$Symbol,
#   y=log10(temp$xgscore+1) - log10(temp$numcitations+1)
# )
# xgboost_residual <- xgboost_residual[order(xgboost_residual$y),]
# xgboost_residual$x <- 1:nrow(xgboost_residual)




############################################################
###### merge it all
############################################################


mdat <- unique(merge(merge(merge(merge(dat2,dat),dat3,all.x = TRUE),dat4,all.x = TRUE),all_dat5))
mdat$firstyear[is.na(mdat$firstyear)] <- -10
mdat$numcitations[is.na(mdat$numcitations)] <- 0
mdat$numcitations <- as.integer(mdat$numcitations)

colnames(mdat) <- c("symbol","nih_geneid","ensembl","description","firstyear","numcitations","xgscore","ct")

library(RSQLite)
con <- dbConnect(SQLite(), dbname = "website/data/geneinfo.sqlite")
dbWriteTable(con, "geneinfo", mdat, overwrite=TRUE)
dbDisconnect(con)


#need to undo the scaling!

############################################################
######################### Genelist #########################
############################################################

write.table(
  unique(data.frame(
    Ensembl.Gene.ID=mdat$ensembl,
    Associated.Gene.Name=mdat$symbol
  )), "website/data/genelist.csv", row.names = FALSE, quote = TRUE, sep=",")




############################################################
## Histogram of papers for each gene
############################################################


df <- feature_pmidcount
#df <- feature_pmidcount_ct[feature_pmidcount_ct$ct=="T cell",]

pmid_term <- read.csv("features/keywords_pmids.csv",stringsAsFactors = FALSE)[,-1]
colnames(pmid_term) <- c("ct","pmid")
head(pmid_term)

g2phm_genespresent<- read.csv(file = "input/g2phm_genespresent.csv",header = T, stringsAsFactors = FALSE)[,-1]
colnames(g2phm_genespresent) <- c("mouseid","gene","pmid")

pub_year <- read.csv("input/pubyear.csv", header=F, stringsAsFactors = F, sep="\t")[,1:2] 
colnames(pub_year)<- c("pmid", "year")

df_geneyear <- sqldf("select count(pmid) as citationcount, mouseid as ensembl, year from pmid_term natural join g2phm_genespresent natural join pub_year group by ensembl,year")
mat_geneyear <- cast(df_geneyear,ensembl~year, value = "citationcount",fill = 0)  
rownames(mat_geneyear) <- mat_geneyear$ensembl
mat_geneyear <- mat_geneyear[,-1]

#Cumulative sum
matcum_geneyear <- mat_geneyear*0
matcum_geneyear[,1] <- mat_geneyear[,1]
for(i in 2:ncol(matcum_geneyear)){
  matcum_geneyear[,i] <- mat_geneyear[,i] + matcum_geneyear[,i-1]
}
matcum_geneyear$ensembl <- rownames(matcum_geneyear)

matcum_geneyear <- melt(as.matrix(matcum_geneyear))
colnames(matcum_geneyear) <- c("ensembl","year","citations")
con <- dbConnect(SQLite(), dbname = "website/data/citations_per_year.sqlite")
dbWriteTable(con, "citationsperyear", matcum_geneyear, overwrite=TRUE)
dbDisconnect(con)


