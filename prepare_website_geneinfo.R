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





##############################
## Histogram of papers for each gene


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

# import plotly.express as px
# df = px.data.tips()
# fig = px.histogram(df, x="total_bill")
# fig.show()


# 
# 
# dcc.Graph(
#   id='papertime-histogram',
#   figure={
#     'data': [
#       {
#         'x': df['year'],
#         'text': df['STRCITY'],
#         'customdata': df['storenum'],
#         'name': 'Open Date',
#         'type': 'histogram'
#       }
#       ],
#     'layout': {}
#   }
# ),
