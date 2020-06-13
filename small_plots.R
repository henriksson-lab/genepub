## feature_qc
library(ggplot2)
library(sqldf)
library(gridExtra)
library(reshape)


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


###############################################################################################################
########## Produce plot: Histograms of paper distributions ####################################################
###############################################################################################################


feature_pmidcount <- read.csv("features/feature_ranked_pmid.csv", stringsAsFactors=FALSE)

### First the two histograms
df <- feature_pmidcount
p1 <- ggplot(df, aes(x=rank_pmid)) + 
  geom_histogram(bins = 50,color="blue", fill="blue")+
  xlab("#Citations (all)")+
  ylab("Frequency")

df <- feature_pmidcount_ct[feature_pmidcount_ct$ct=="T cell",]
p2 <- ggplot(df, aes(x=rank_pmid)) + 
  geom_histogram(bins = 15,color="darkgreen", fill="darkgreen")+
  xlab("#Citations (T cell)")+
  ylab("Frequency")

ggp <- grid.arrange(p1,p2, nrow=1)
ggsave("plots/paper_histograms.pdf", ggp, width = 6, height = 3)


### Now the pareto fit, and other fits

library(ParetoPosStable)
df <- feature_pmidcount_ct[feature_pmidcount_ct$ct=="T cell",]
pdf("plots/fit_pareto.pdf")
plot(pareto.fit(10^df$rank_pmid-1))
dev.off()

require(MASS)
ex <- 10^df$rank_pmid-1
fit1 <- fitdistr(ex, "exponential") 
df("plots/fit_exp.pdf")
hist(ex, freq = FALSE, breaks = 500, xlim = c(0, quantile(ex, 0.99)))
curve(dexp(x, rate = fit1$estimate), from = 0, col = "red", add = TRUE)
dev.off()

###############################################################################################################
########## Produce plot: Expression vs papers #################################################################
###############################################################################################################

#### Comparison expression vs paper count

cell_type <- "T cell"

dat <- merge(feature_pmidcount, feature_exp)
feat_orig <- dat[dat$ct==cell_type & dat$rank_pmid>1, ]
# plot(
#   feat_orig$rank_exp,
#   feat_orig$rank_pmid,cex=0.5)
feat_orig$density <- get_density(feat_orig$rank_exp, feat_orig$rank_pmid, n = 100)
gp <- ggplot(feat_orig, aes(rank_exp, rank_pmid, color=density)) +
  geom_point() +
  xlab("Normalized RNA Expression level") +
  ylab("Normalized Citations") +
  scale_color_viridis_c()
gp
ggsave("plots/final_exp_vs_pmid.pdf", gp)




###############################################################################################################
########## Produce plot: GWAS vs papers #######################################################################
###############################################################################################################



######################### Just the trend line
dat <- merge(feature_pmidcount, feature_gwas)
feat_orig <- dat[dat$rank_pmid>1, ]
feat_orig$density <- get_density(feat_orig$rank_gwas, feat_orig$rank_pmid, n = 20)
# plot(
#   feat_orig$rank_gwas,
#   feat_orig$rank_pmid,cex=0.5)
gp <- ggplot(feat_orig, aes(rank_gwas, rank_pmid)) +  #, color=density
  geom_point() +
  xlab("Normalized GWAS significance") +
  ylab("Normalized Citations") +
  #  scale_color_viridis_c() +
  #  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x)
gp
ggsave("plots/final_gwas_vs_pmid.pdf", gp)






#########################################
#########################################
#########################################

# dat <- merge(merge(feature_founder_rank, feature_minyear_gene_ct), feature_pmidcount)
# #dim(dat)
# 
# is.na(dat$family_founder_rank)
# 
# sqldf("select family_founder_rank")
# dat$family_indexdiff
# dat$gene
# dat$family_founder_rank
# 
# ggplot
# dat$first_year
# 
# 
# 
# 
# 


###############################################################################################################
########## Produce plot: #paper/#gene vs years ################################################################
###############################################################################################################

pub_year <- read.csv("input/pubyear.csv", header=F, stringsAsFactors = F, sep="\t")[,1:2] 
colnames(pub_year)<- c("pmid", "year")

papers_per_year <- sqldf("select count(pmid) as count_paper, year from pub_year group by year")
genes_per_year <- sqldf("select distinct first_year as year, count(gene) as count_gene  from feature_minyear_gene_ct group by year")

dat <- merge(papers_per_year, genes_per_year, all=TRUE)
dat$count_paper[is.na(dat$count_paper)] <- 0
dat$count_gene[is.na(dat$count_gene)] <- 0


dat$count_paper <- log10(1+dat$count_paper)
dat$count_gene <- log10(1+dat$count_gene)

dat <- dat[dat$year>1900 & dat$year<=2018,]
dat <- dat[order(dat$year),]

plot(dat$year, dat$count_paper, type="l")
plot(dat$year, dat$count_gene, type="l")


g1 <- ggplot(dat, aes(x=year, y=count_paper)) + 
  ylab("Log10 1+Number of publications") +
  geom_area(color="darkblue",
            fill="lightblue")
g2 <- ggplot(dat, aes(x=year, y=count_gene)) + 
  ylab("Log10 1+Number of new genes") +
  geom_area(color="#C4961A",
            fill="#FFDB6D")

ggp <- grid.arrange(g1,g2, nrow=2)
ggsave("plots/final_year_vs_papers.pdf", ggp)


###############################################################################################################
########## Produce plot: #paper vs chromosome #################################################################    # for playing around
###############################################################################################################

feature_exp <- read.csv("features/feature_geneexp.csv", stringsAsFactors=FALSE)
feature_chromloc <- read.csv("features/feature_chromloc.csv", stringsAsFactors=FALSE) ####################### remove old feature!
feature_minyear_gene_ct <- read.csv("features/feature_minyear_gene_ct.csv", stringsAsFactors = FALSE)
chromdat <- read.csv("plots/data_chromloc.csv", stringsAsFactors = FALSE)
#feature_founder_rank <- read.csv("features/feature_founder_fam_rankpmid.csv", stringsAsFactors = FALSE)

feature_founder_rank

### Merge data
dat <- merge(chromdat, feature_minyear_gene_ct)
#dat <- merge(dat, feature_founder_rank, all.x = TRUE)
cur_chrom <- "7"
dat <- dat[dat$chrom==cur_chrom,] #10 cool
dat$genetype <- "ZZZ"

#chrom 4 has a dip around 400

### Compensate for first_year
onelm <- lm(rank_pmid~first_year, dat)
dat$adj_pmid <- dat$rank_pmid - onelm$fitted.values
dat$adj_pos <- rank(dat$pos)
dat <- dat[order(dat$adj_pos),]


####################
#for 9, olfr
# dat[dat$adj_pos>280 & dat$adj_pos<320,]
# dat[dat$adj_pos>70 & dat$adj_pos<80,]
# dat[dat$adj_pos>800 & dat$adj_pos<900,]  #10
# dat[dat$adj_pos>380 & dat$adj_pos<400,]  #10
# dat[dat$rank_pmid<1.2,]  #10
dat[dat$rank_pmid<1,]  #8
# dat[dat$rank_pmid>2.5,]



dat$genetype[grep("Olfr",dat$gene)] <- "Olfr"
dat$genetype[grep("Vmn",dat$gene)] <- "Vmn"
dat$genetype[grep("Taar",dat$gene)] <- "Taar"
dat$genetype[grep("Zfp",dat$gene)] <- "Zfp"
dat$genetype[grep("Defa",dat$gene)] <- "Defa"   #more studied than the rest
dat$genetype[grep("Defb",dat$gene)] <- "Defb"   #https://en.wikipedia.org/wiki/Beta_defensin   chrom8
dat$genetype[grep("Psg",dat$gene)] <- "Psg"     #chrom 7, many
dat$genetype[grep("Obox",dat$gene)] <- "Obox"  #4
dat$genetype[grep("Zscan",dat$gene)] <- "Zscan"
#dat$genetype <- factor(dat$genetype)


#chrom1
dat$genetype[grep("Ugt",dat$gene)] <- "Ugt" #decent size
dat$genetype[grep("Fcgr",dat$gene)] <- "Fcgr"
dat$genetype[grep("Sel",dat$gene)] <- "Sel"
dat$genetype[grep("Il",dat$gene)] <- "Il"   #decent, localized

#chrom2
dat$genetype[grep("Lcn",dat$gene)] <- "Lcn"    #chrom 2, low
dat$genetype[grep("Wfdc",dat$gene)] <- "Wfdc" 

#chrom3
dat$genetype[grep("Lce",dat$gene)] <- "Lce" #skin
dat$genetype[grep("Tdpoz",dat$gene)] <- "Tdpoz"  #in embryos hard to study 
dat$genetype[grep("Sprr",dat$gene)] <- "Sprr" #skin

#chrom4
dat$genetype[grep("Mup",dat$gene)] <- "Mup"


#Tdpoz
#Sprr




dat[dat$adj_pmid < 1.2,]  #10
tail(dat[dat$adj_pmid < 1.2,],n=100)

#View(dat[dat$adj_pmid < 1,])
#10
#tail(dat[dat$adj_pmid < 1,],n=100)


#g1 <- 
ggplot(dat, aes(x=adj_pos, y=rank_pmid, color=genetype)) + 
  ylab("Normalized Citations")+
  xlab(sprintf("Rank[ Chromosome %s position ]",cur_chrom))+
  geom_point()




ggplot(dat, aes(x=adj_pos, y=adj_pmid, color=genetype)) + 
  #ggplot(dat, aes(x=pos, y=rank_pmid)) + 
  #ylab("Log10 1+Number of publications") +
  ylab("Normalized Citations") +
  xlab("Rank[ Chromosome 10 position ]") +
  geom_point()



############# Merge gene families. Do we still see patterns?  

dat_f <- dat
dat_f$common_name <- str_extract(dat_f$gene, "[a-zA-Z]*") #letters until first number
dat_f <- sqldf("select distinct common_name, avg(pos) as pos, avg(rank_pmid) as rank_pmid, genetype from dat_f group by common_name order by pos")
dat_f$adj_pos <- rank(dat_f$pos)
dat_f$rand_pos <- sample(dat_f$adj_pos, nrow(dat_f))

ggplot(dat_f, aes(x=adj_pos, y=rank_pmid)) + 
  ylab("Normalized Citations")+
  xlab("Rank[ Chromosome position, avg family ]")+
  geom_point()+
  geom_smooth()

ggplot(dat_f, aes(x=rand_pos, y=rank_pmid)) +
  ylab("Normalized Citations")+
  xlab("Random position")+
  geom_point()+
  geom_smooth()

### some spatial statistic test

# alli <- 1:20
# acorr <- NULL
# for(i in alli){
#   v1 <- dat_f$rank_pmid
#   v2 <- lag(dat_f$rank_pmid,i)
#   v1 <- v1[20:(length(v1)-20)]
#   v2 <- v2[20:(length(v2)-20)]
#   acorr <- c(acorr, cor(v1,v2))
# }
# plot(acorr,type="l")
#   
# plot(
#   dat_f$rank_pmid,
#   lag(dat_f$rank_pmid,3))
# 

dat_f[dat_f$rank_pmid>3,]
#chrom 11, immune genes toward the end


#more citations toward the ends of chromosomes?




#genes_per_year <- sqldf("select distinct first_year as year, count(gene) as count_gene  from feature_minyear_gene_ct group by year")

dat <- merge(fpapers_per_year, genes_per_year, all=TRUE)

chromdat <- read.csv("plots/data_chromloc.csv", stringsAsFactors = FALSE)


dat_newfam <- dat
dat_newfam$common_name <- str_extract(dat_newfam$gene, "[a-zA-Z]*") #letters until first number
dat_f <- sqldf("select distinct common_name, avg(pos) as pos, sum(rank_pmid) as sum_rank_pmid from dat_f group by common_name order by pos")
dat_f$adj_pos <- rank(dat_f$pos)
dat_f$rand_pos <- sample(dat_f$adj_pos, nrow(dat_f))




###############################################################################################################
########## Produce plot: #paper vs chromosome #################################################################    # final plot
###############################################################################################################

feature_exp <- read.csv("features/feature_geneexp.csv", stringsAsFactors=FALSE)
feature_chromloc <- read.csv("features/feature_chromloc.csv", stringsAsFactors=FALSE) ####################### remove old feature!
feature_minyear_gene_ct <- read.csv("features/feature_minyear_gene_ct.csv", stringsAsFactors = FALSE)
chromdat <- read.csv("plots/data_chromloc.csv", stringsAsFactors = FALSE)


### Merge data
dat <- merge(chromdat, feature_minyear_gene_ct)
cur_chrom <- "7"
dat <- dat[dat$chrom==cur_chrom,] #10 cool
dat$genetype <- "Other"


### Compensate for first_year
onelm <- lm(rank_pmid~first_year, dat)
dat$adj_pmid <- dat$rank_pmid - onelm$fitted.values
dat$adj_pos <- rank(dat$pos)
dat <- dat[order(dat$adj_pos),]



dat$genetype[grep("Olfr",dat$gene)] <- "Olfr"
dat$genetype[grep("Vmn",dat$gene)] <- "Vmn"
dat$genetype[grep("Taar",dat$gene)] <- "Taar"
dat$genetype[grep("Zfp",dat$gene)] <- "Zfp"
#dat$genetype[grep("Defa",dat$gene)] <- "Defa"   #more studied than the rest
#dat$genetype[grep("Defb",dat$gene)] <- "Defb"   #https://en.wikipedia.org/wiki/Beta_defensin   chrom8
dat$genetype[grep("Psg",dat$gene)] <- "Psg"     #chrom 7, many
#dat$genetype[grep("Obox",dat$gene)] <- "Obox"  #4
#dat$genetype[grep("Zscan",dat$gene)] <- "Zscan"
#dat$genetype <- factor(dat$genetype)


ggplot(dat, aes(x=adj_pos, y=rank_pmid, color=genetype)) + 
  ylab("Normalized Citations")+
  xlab(sprintf("Rank[ Chromosome %s position ]",cur_chrom))+
  ylab("Log10 1+Number of publications") +
  geom_point()
ggsave("plots/chrom7.pdf", width = 7, height = 3)

dat[dat$rank_pmid>3.1,]


###############################################################################################################
########## Produce plot: #paper vs essentiality ###############################################################
###############################################################################################################


feature_essentiality <- read.csv("features/feature_essentiality_global.csv", stringsAsFactors=FALSE)


dat <- merge(feature_pmidcount, feature_essentiality)
#feat_orig <- dat[dat$ct==cell_type & dat$rank_pmid>1, ]
feat_orig <- dat[dat$rank_pmid>1, ]


feat_orig$density <- get_density(feat_orig$essentiality_global, feat_orig$rank_pmid, n = 100)
gp <- ggplot(feat_orig, aes(essentiality_global, rank_pmid, color=density)) +
  geom_point() +
  xlab("Essentiality") +
  ylab("Normalized Citations") +
  #geom_smooth(method='lm', formula= y~x) +
  scale_color_viridis_c()
gp
ggsave("plots/final_essentiality_vs_pmid.pdf", gp)





###############################################################################################################
########## Produce plot: #paper vs gene family ################################################################
###############################################################################################################


feature_founder_rank <- read.csv("features/feature_founder_fam_rankpmid.csv", stringsAsFactors = FALSE)
feature_founder_rank <- feature_founder_rank[feature_founder_rank$family_indexdiff<=10,]
feature_minyear_gene_ct <- read.csv("features/feature_minyear_gene_ct.csv", stringsAsFactors = FALSE)

dat <- merge(merge(feature_pmidcount, feature_founder_rank),feature_minyear_gene_ct)

#Get year of the founder
dat <- merge(
  dat, 
  data.frame(
    founder_name=feature_minyear_gene_ct$gene,
    founder_year=feature_minyear_gene_ct$first_year)
)

#Get relative year
dat$relative_year <- dat$first_year-dat$founder_year

#Get the relative rank
dat$rel_rank <- dat$rank_pmid - dat$family_founder_rank

#Regress out the year dependence
ylm <- lm(rel_rank~relative_year,dat)
dat$rel_rank_noyear <- dat$rel_rank - dat$relative_year*ylm$coefficients[2]

# plot(dat$family_indexdiff, dat$rel_rank)
# plot(dat$family_indexdiff, dat$first_year)  #strong dependence
# plot(dat$family_indexdiff, dat$relative_year)  #some dependence

dat2 <- rbind(
  sqldf("select family_indexdiff, avg(rel_rank) as rel_rank, 'Uncorrected' as grp from dat group by family_indexdiff"),
  sqldf("select family_indexdiff, avg(rel_rank_noyear) as rel_rank, 'Year-corrected' as grp from dat group by family_indexdiff"))

gp <- ggplot(dat2, aes(family_indexdiff, rel_rank, fill=grp)) +
  geom_bar(stat="identity", position=position_dodge())
#  geom_line()
gp


dat3 <- sqldf("select family_indexdiff, avg(rel_rank) as rel_rank, 'Uncorrected' as grp from dat group by family_indexdiff")
dat4 <- sqldf("select family_indexdiff, avg(rel_rank_noyear) as rel_rank, 'Year-corrected' as grp from dat group by family_indexdiff")
g1 <- ggplot(dat3, aes(family_indexdiff, rel_rank)) +
  geom_bar(stat="identity", position=position_dodge(), color="black",fill="lightblue") +
  ylab("Relative Log10 1+Number of publications") +
  ylim(-0.5, 0)
#  geom_area(color="darkblue",
#            fill="lightblue")
g2 <- ggplot(dat4, aes(family_indexdiff, rel_rank)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",fill="#FFDB6D") +
  ylab("Relative Log10 1+Number of publications") +
  ylim(-0.5, 0)
ggp <- grid.arrange(g1,g2, nrow=2)
ggsave("plots/family_rel_papers.pdf", ggp)






###############################################################################################################
########## Produce plot: umap of genes ########################################################################
###############################################################################################################

feature_exp <- read.csv("features/feature_geneexp.csv", stringsAsFactors=FALSE)
dat_coexp <- read.csv("plots/data_umap_coexp_k10_T cell.csv")
colnames(dat_coexp) <- c("gene","x","y")

dat <- merge(merge(feature_pmidcount, feature_exp[feature_exp$ct=="T cell",]), dat_coexp)
dat <- dat[order(dat$rank_pmid, decreasing = FALSE),]

g1 <- ggplot(dat, aes(x, y, color=rank_pmid)) +
  geom_point(size=2) +
  xlab("UMAP x") +
  ylab("UMAP y") +
  scale_color_gradientn(colours = c("#000000","#000000","#FFFFFF"))
#  scale_color_viridis_c()
g2 <- ggplot(dat, aes(x, y, color=rank_exp)) +
  geom_point(size=2) +
  xlab("UMAP x") +
  ylab("UMAP y") +
  scale_color_viridis_c()
ggp <- grid.arrange(g1,g2, ncol=2)
ggsave("plots/one_coexp.pdf", ggp, width = 10, height = 4)


w <- dat[dat$y>16,]
w[order(w$y),]


###############################################################################################################
########## Produce plot: histogram of gene coverage ###########################################################
###############################################################################################################


dat <- feature_pmidcount
dat$rank_pmid <- 10^(dat$rank_pmid)-1

hist(dat$rank_pmid,breaks=100)


p<-ggplot(dat, aes(x=rank_pmid)) + 
  geom_histogram(color="black", fill="white") +
  scale_x_log10() +
  xlab("Genes") +
  ylab("Number of publications") 
p
ggsave("plots/histogram_all_genes.pdf", p, width = 3, height = 3)





###############################################################################################################
########## Produce plot: XXX feature vs time ##################################################################
###############################################################################################################




# ############# Expression
# dat <- merge(feature_exp[feature_exp$ct=="T cell",], feature_minyear_gene_ct)
# plot_tscore_time(dat,"rank_exp")   #quite unaffected
# 
# plot(dat$first_year, dat$rank_exp)
# dat$density <- get_density(dat$first_year, dat$rank_exp, n = 100)
# gp <- ggplot(dat, aes(dat$first_year, dat$rank_exp, color=density)) +
#   geom_point() +
# #  xlab("Normalized RNA Expression level") +
# #  ylab("Normalized Citations") +
#   scale_color_viridis_c()
# gp
# 
# 
# ############# Essentiality
# dat <- merge(feature_essentiality, feature_minyear_gene_ct)
# plot_tscore_time(dat,"essentiality_global")   #essentiality direction depends on the weight. up if not weighted. 
# 
# plot(dat$first_year, dat$essentiality_global)
# dat$density <- get_density(dat$first_year, dat$essentiality_global, n = 100)
# gp <- ggplot(dat, aes(dat$first_year, dat$essentiality_global, color=density)) +
#   geom_point() +
#   #  xlab("Normalized RNA Expression level") +
#   #  ylab("Normalized Citations") +
#   scale_color_viridis_c()
# gp
# 
# 
# 
# 
# 
# ############# Cosmic
# dat <- merge(feature_cosmic, feature_minyear_gene_ct)
# plot_tscore_time(dat,"rank_cosmic")  #weird stuff around 2000 ... going down
# 
# plot(dat$first_year, dat$rank_cosmic)
# dat$density <- get_density(dat$first_year, dat$rank_cosmic, n = 100)
# gp <- ggplot(dat, aes(first_year, rank_cosmic, color=density)) +
#   geom_point() +
#   #  xlab("Normalized RNA Expression level") +
#   #  ylab("Normalized Citations") +
#   scale_color_viridis_c()
# gp
# 
# 
# 
# 
# ############# gwas
# dat <- merge(feature_gwas, feature_minyear_gene_ct)
# plot_tscore_time(dat,"rank_gwas")  #trend downwards   ...
# 
# plot(dat$first_year, dat$rank_gwas)    ##### why 2 gwas ???? problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# dat$density <- get_density(dat$first_year, dat$rank_gwas, n = 100)
# gp <- ggplot(dat, aes(first_year, rank_gwas, color=density)) +
#   geom_point() +
#   #  xlab("Normalized RNA Expression level") +
#   #  ylab("Normalized Citations") +
#   scale_color_viridis_c()
# gp
# 



feature_founder_rank <- read.csv("features/feature_founder_fam_rankpmid.csv", stringsAsFactors = FALSE)
feature_exp <- read.csv("features/feature_geneexp.csv", stringsAsFactors=FALSE)
feature_essentiality <- read.csv("features/feature_essentiality_global.csv", stringsAsFactors=FALSE)
feature_cosmic <- read.csv("features/feature_cosmic.csv", stringsAsFactors = FALSE)
feature_gwas <- read.csv("features/feature_gwas.csv", stringsAsFactors = FALSE)



#############
plot_tscore_time <- function(dat, fname,yname=fname){
  wdat <- data.frame(
    year=dat$first_year,
    f=dat[,fname]#$rank_gwas
  )
  wdat$f <- (wdat$f - mean(wdat$f))/sd(wdat$f)
  
  wdat <- sqldf("select avg(f) as avgf, count(f) as cf, year from wdat group by year")
  wdat$norm <- wdat$avgf#*sqrt(wdat$cf)
  
  gp <- ggplot(wdat, aes(year, norm)) +
    geom_point(aes(size=cf)) + 
    geom_line() +
    geom_smooth(method='loess', formula= y~x, mapping = aes(weight = cf)) +
    scale_x_continuous(breaks=seq(1900,2020,20),limits = c(1900,2010)) +
    ylab(yname)
  gp
}
# 
# 
# ####################################### Only using known values
# 
# ############# Summary
# dat <- merge(feature_exp[feature_exp$ct=="T cell",], feature_minyear_gene_ct)
# p1 <- plot_tscore_time(dat,"rank_exp","Expression")   #quite unaffected
# dat <- merge(feature_essentiality, feature_minyear_gene_ct)
# p2 <- plot_tscore_time(dat,"essentiality_global","Essentiality")   #essentiality direction depends on the weight. up if not weighted. 
# dat <- merge(feature_cosmic, feature_minyear_gene_ct)
# p3 <- plot_tscore_time(dat,"rank_cosmic","COSMIC")  #weird stuff around 2000 ... going down
# dat <- merge(feature_gwas, feature_minyear_gene_ct)
# p4 <- plot_tscore_time(dat,"rank_gwas","GWAS")  #trend downwards   ...
# dat <- merge(feature_founder_rank, feature_minyear_gene_ct)
# p5 <- plot_tscore_time(dat,"family_indexdiff","Index")
# # p6 <- plot_tscore_time(dat,"family_founder_rank","FR") #goes down trivially
# # p6
# ggp <- grid.arrange(p1,p2,p3,p4,p5, ncol=1)
# ggp
# 
# ggsave("plots/feature_vs_time.pdf", ggp,width = 5,height=10)

####################################### Imputed

dat <- totfeature[totfeature$ct=="T cell",]
dat <- dat[,c("rank_exp","essentiality_global","rank_cosmic","rank_gwas","family_indexdiff", "first_year")]
dat$family_indexdiff[is.na(dat$family_indexdiff)] <- 0
dat$essentiality_global[is.na(dat$essentiality_global)] <- median(dat$essentiality_global,na.rm = TRUE)
dat$rank_cosmic[is.na(dat$rank_cosmic)] <- min(dat$rank_cosmic,na.rm = TRUE)
dat$rank_gwas[is.na(dat$rank_gwas)] <- min(dat$rank_gwas,na.rm = TRUE)

############# Summary
p1 <- plot_tscore_time(dat,"rank_exp","Expression")   #quite unaffected
p2 <- plot_tscore_time(dat,"essentiality_global","Essentiality")   #essentiality direction depends on the weight. up if not weighted. 
p3 <- plot_tscore_time(dat,"rank_cosmic","COSMIC")  #weird stuff around 2000 ... going down
p4 <- plot_tscore_time(dat,"rank_gwas","GWAS")  #trend downwards   ...
p5 <- plot_tscore_time(dat,"family_indexdiff","Index")
# p6 <- plot_tscore_time(dat,"family_founder_rank","FR") #goes down trivially
# p6
ggp <- grid.arrange(p1,p2,p3,p4,p5, ncol=2)
ggp

ggsave("plots/feature_vs_time.pdf", ggp,width = 10,height=10)





###############################################################################################################
########## Produce plot: Expression level over tissues correlation with papers ################################
###############################################################################################################

library(lineup)


feature_pmidcount_ct <- read.csv("features/feature_ranked_pmid_ct.csv", stringsAsFactors = FALSE)
unique(feature_pmidcount_ct$ct)
dat <- merge(feature_exp,feature_pmidcount_ct)
#dat$pmid_count <- log10(1+dat$pmid_count)

dat1 <- cast(dat[,c("ct","gene","rank_exp")], ct~gene, value="rank_exp")
dat2 <- cast(dat[,c("ct","gene","pmid_count")], ct~gene, value="pmid_count") 
datcor <- corbetw2mat(
  dat1[,-1],
  log10(1+dat2[,-1]))

dat2[,-1]

df <- data.frame(
  cor=datcor,
  gene=colnames(dat1)[-1],
  sumpaper=log10(1+apply(dat2[,-1],2,sum))
)

df$text <- as.character(df$gene)
df$text[df$density>0.007] <- ""

df$density <- get_density(df$cor, df$sumpaper, n = 100)
gp <- ggplot(df, aes(cor, sumpaper, color=density, label=text)) +
  geom_point() +
  xlab("Correlation Expression vs Citations") +
  ylab("Log10 1+Number of citations") +
  geom_text(color="red") +
  scale_color_viridis_c()
gp
ggsave("plots/correlation_exp_citation_ct.pdf", gp, width = 4, height = 4)

median(df$cor)
median(df$sumpaper)

#More than 10 papers ideally
#10^1-1     =9
#10^1.2-1   =15
#10^1.4-1   =24  *****
#10^1.5-1   =30



df[df$cor>0.5 & df$sumpaper<1,]


# plot(sort(df$cor))
# plot(sort(df$cor/df$sumpaper))
# plot(df$cor,df$sumpaper)


# df$ncor <- df$cor / apply(datcor,1,sd)
# plot(sort(df$ncor))

df <- df[order(df$cor, decreasing = TRUE),]
df
#acor <- diag(datcor) #/ apply(datcor,1,sd)

plot_gene_vs_ct <- function(gname){
  sdat <- dat[dat$gene==gname,]
  plot(sdat$rank_exp, sdat$pmid_count,cex=0)
  text(sdat$rank_exp, sdat$pmid_count, label=sdat$ct)
  
  ggplot(sdat,aes(rank_exp, pmid_count,label=ct)) +
    xlab("Expression") +
    ylab("Log10 1+Number of citations") +
    geom_point(color="red") +
    geom_text(angle=90, hjust="left")
}

# sdat <- dat[dat$gene=="Foxk1",]
# plot(sdat$rank_exp, sdat$pmid_count,cex=0)
# text(sdat$rank_exp, sdat$pmid_count, label=sdat$ct)
# 
# ggplot(sdat,aes(rank_exp, pmid_count,label=ct)) +
#   xlab("Expression") +
#   ylab("Log10 1+Number of citations") +
#   geom_point(color="red") +
#   geom_text(angle=90, hjust="left")
plot_gene_vs_ct("Foxk1")
ggsave("plots/ct_exp_corr_Foxk1.pdf", width = 5, height = 3)

plot_gene_vs_ct("Pou5f1")
ggsave("plots/ct_exp_corr_Pou5f1.pdf", width = 5, height = 3)









