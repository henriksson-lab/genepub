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

  
mdat <- unique(merge(merge(merge(dat2,dat),dat3,all.x = TRUE),dat4,all.x = TRUE))
mdat$firstyear[is.na(mdat$firstyear)] <- -10
mdat$numcitations[is.na(mdat$numcitations)] <- 0
mdat$numcitations <- as.integer(mdat$numcitations)

colnames(mdat) <- c("symbol","nih_geneid","ensembl","description","firstyear","numcitations")

library(RSQLite)
con <- dbConnect(SQLite(), dbname = "website/data/geneinfo.sqlite")
dbWriteTable(con, "geneinfo", mdat, overwrite=TRUE)
dbDisconnect(con)


