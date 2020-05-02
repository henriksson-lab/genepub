library(reshape2)
library(rentrez)
library(irlba)
library(dplyr)
library(sqldf)
library(stringr)
library(umap)
library(smooth)
library(plotly)
library(Matrix)
library(stringr)
library(Rtsne)
library(FNN)
library(RSQLite)



#############################################################
## Function: Rank a column within each cell type (column ct)
rank_by_ct <- function(dat, col){
  for(ct in unique(dat$ct)){
    dat[dat$ct==ct,col] <- rank(dat[dat$ct==ct,col])
  }
  dat[,col]
}


#############################################################
## Function: Store coordinates for the website
store_website_coordinates <- function(name,coord) {

  ## Integrate mouse ENSMUSG* IDs  
  map_ensmus_sym <- read.csv("input/map_ensmus_symbol.csv", stringsAsFactors = FALSE)
  colnames(map_ensmus_sym) <- c("ensmus","gene")

  ##Write to file
  con <- dbConnect(SQLite(), dbname = sprintf("website/data/coord_%s.sqlite", name))
  dbWriteTable(con, "coord", merge(coord, map_ensmus_sym), overwrite=TRUE)
  dbDisconnect(con)
}


#############################################################
########## How many papers per gene and cell type? ##########
#############################################################


### Read map: PMID - mouse ensmusg - mouse symbol
g2phm_genespresent<- read.csv(file = "input/g2phm_genespresent.csv",header = T, stringsAsFactors = FALSE)[,-1]
colnames(g2phm_genespresent) <- c("mouseid","gene","pmid")

if(file.exists("features/keywords_pmids.csv")){
  
  ### Read map: keyword - PMID
  pmid_term <- read.csv("features/keywords_pmids.csv",stringsAsFactors = FALSE)[,-1]
  
} else {

  #Get list of keywords we use to search for papers about cell types: c("celltype","collapsed")
  keywords_master_dict <- read.csv("input/celltypes_tissue_to_collapse.csv",stringsAsFactors = F)[,3:4]
  keywords_list_for_query <- unique(keywords_master_dict$collapsed)
  
  #Query pubmed for each keyword
  pmid_term <- NULL
  for (i in 1:length(keywords_list_for_query)){
    print(i)
    
    #Get how many entries
    tst_count<- entrez_search(db="pubmed", term=keywords_list_for_query[i], email="iose0001@student.umu.se")$count
    print(tst_count)
    
    #Get the entries
    tst2<- entrez_search(db="pubmed", term=keywords_list_for_query[i], email="iose0001@student.umu.se", retmax= tst_count)
    #pmids[[i]] <- tst2$ids
    
    pmid_term <- rbind(
      pmid_term,
      data.frame(
        term=keywords_list_for_query[i],   ############ beep. is this correct? should use collapsed type
        pmid=tst2$ids, 
        stringsAsFactors = FALSE)
    )
    
  }
  
  ### Write map: keyword - PMID
  write.csv(pmid_term,"features/keywords_pmids.csv")
}

############# Create #paper for each gene and cell type


if(file.exists("features/feature_ranked_pmid.csv")) {
  
  ### Load from disk
  feature_pmidcount <- read.csv("features/feature_ranked_pmid.csv", stringsAsFactors=FALSE)
  
} else {
  
  ### Calculate #paper / ct,gene
  ranked_pmid_ct_gene <- merge(pmid_term, g2phm_genespresent)
  feature_pmidcount <- sqldf("select gene, term as ct, count(pmid) as pmid_count from ranked_pmid_ct_gene group by gene, term")
  
  ### Store  
  write.csv(feature_pmidcount, "features/feature_ranked_pmid.csv", row.names = F)

}


############# Calculate rank #paper

feature_pmidcount$rank_pmid <- rank_by_ct(feature_pmidcount, "pmid_count")
feature_pmidcount <- within(feature_pmidcount, rm(pmid_count))






#############################################################
######### Common code for RNAseq data #######################
#############################################################



### Read what cell types we reduce to
update_for_FACS_onth <- unique(read.csv(file = "input/celltypes_tissue_to_collapse.csv", stringsAsFactors = FALSE)[,-(1:2)])
colnames(update_for_FACS_onth)[1]<- "cell_ontology_class"

### Read cell type annotation
FACS_meta<- read.csv("input/tabula_muris/annotations_facs.csv", header = T, stringsAsFactors = FALSE, na.strings=c("","NA"))
FACS_onthology<- FACS_meta[,3:4]
FACS_onthology <- FACS_onthology[!is.na(FACS_onthology$cell_ontology_class),]
new_FACS_ontology <- merge(FACS_onthology, update_for_FACS_onth)

dim(new_FACS_ontology)


### Read one count table into a sparse matrix; normalize counts
count_norm_colmeans<- function(x){
  input<-read.csv(file=paste("input/tabula_muris/counts/",x, sep = ""), header=T, stringsAsFactors = F, row.names = 1)
  input_mtx <- as.matrix(input)
  input_sparse_mtx = Matrix(input_mtx,sparse = T)
  cm <- colSums(input_sparse_mtx)
  normcounts <- t(t(input_sparse_mtx)/cm)
  return(normcounts)
}

### Read all count tables into memory
FACS_norm_arrays<- list()
for(tissue in list.files("input/tabula_muris/counts")){
  renamed_tissue <- str_split_fixed(tissue,"-counts",2)[1]
  print(renamed_tissue)
  FACS_norm_arrays[[renamed_tissue]]<- count_norm_colmeans(tissue)
}


#############################################################
### Figure out which tissue a cell type is most common in ###
#############################################################

### Count #cells for each cell type across tissues
tissue_celltype_count <- NULL
for(n in names(FACS_norm_arrays)){
  x <- as.data.frame(table(merge(
    data.frame(cell=colnames(FACS_norm_arrays[[n]]), stringsAsFactors=FALSE),
    new_FACS_ontology)[,c("collapsed")]))
  x$Var1 <- as.character(x$Var1)
  
  tissue_celltype_count <- rbind(
    tissue_celltype_count,
    data.frame(
      ct=x$Var1,
      count=x$Freq,
      tissue=n
    )
  )
}

### Extract which tissue has max # cells for each cell type
tissue_celltype_count_max <- merge(
  sqldf("select ct, max(count) as count, tissue from tissue_celltype_count group by ct"),
  tissue_celltype_count, stringsAsFactors=FALSE)
tissue_celltype_count_max$ct <- as.character(tissue_celltype_count_max$ct)
tissue_celltype_count_max <- tissue_celltype_count_max[order(tissue_celltype_count_max$tissue),]

### Store which tissue has the max #cells for each cell type
write.csv(tissue_celltype_count_max, "features/tabula_tissue_mostabundant_ct.csv", row.names = FALSE)




#############################################################
###### Features from average expression, each cell type #####
#############################################################


if(file.exists("features/feature_geneexp.csv")) {

  ### Load from disk
  feature_exp <- read.csv("features/feature_geneexp.csv", stringsAsFactors=FALSE)
  
} else {

  ### For each cell type, find rank average expression level
  feature_exp <- NULL
  for (i in 1:nrow(tissue_celltype_count_max)){
    the_tissue <- tissue_celltype_count_max$tissue[i]
    the_ct <- tissue_celltype_count_max$ct[i]
    print(sprintf("%s -- %s" , the_tissue, the_ct))
    
    #Extract count table for this cell type  
    #tissue_count <- count_norm_colmeans(tissue)   #if need to load on the fly
    tissue_count <- FACS_norm_arrays[[the_tissue]]
    input <- tissue_count[,colnames(tissue_count) %in% new_FACS_ontology$cell[new_FACS_ontology$collapsed==the_ct]]
    input <- log10(1+input)
    
    ### TODO size factor normalization?
    
    feature_one_exp <- data.frame(
      gene=rownames(input),
      ct=the_ct,
      rank_exp=rank(rowMeans(input))
    )
    
    feature_exp <- rbind(
      feature_exp,
      feature_one_exp
    )
  }

  ### Store feature, rank expression
  write.csv(
    feature_exp[,c("gene","ct","rank_exp")], 
    "features/feature_geneexp.csv", row.names = FALSE)
  
}



#############################################################
######### Features from co-expression #######################
#############################################################



############################################################
## Calculate neighbour avg(PMID), given a graph. 
## graph currently does not have ct, code can be made to use it optionally
calc_neigh_pmid <- function(graph, feature_pmidcount, featurename){
  #symmetrize graph
  graph <- unique(rbind(
    graph,
    data.frame(
      from=graph$to,
      to=graph$from
    )
  ))
  
  ### Removing self interactions
  graph<-graph[graph$from!=graph$to,]
  
  ## Assemble from-to neighbour PMID table. Can do all cell types in one go
  neigh_pmid <- merge(
    data.frame(
      neigh=feature_pmidcount$gene,
      ct=feature_pmidcount$ct,
      neigh_rank_pmid=feature_pmidcount$rank_pmid
    ),
    data.frame(
      gene=graph$from,
      neigh=graph$to
    )
  )
  
  
  #Calculate neighbour average PMID
  neigh_pmid <- sqldf("select gene, ct,  avg(neigh_rank_pmid) as neigh_pmid from neigh_pmid group by gene,ct")
  
  #Rename new feature
  colnames(neigh_pmid)[colnames(neigh_pmid)=="neigh_pmid"] <- featurename
  
  neigh_pmid
}



if(file.exists("features/feature_coexp.csv")) {
  
  ### Read from disk
  feature_coexp <- read.csv("features/feature_coexp.csv", stringsAsFactors = FALSE)
  
} else {
  
  ### Clear count tables to save memory. Load on demand later
  rm(FACS_norm_arrays)
  
  ### For each cell type, find regulatory network (co-expression) neighbours
  for (i in 1:nrow(tissue_celltype_count_max)){   #quick hack
    #i<-17
    the_tissue <- tissue_celltype_count_max$tissue[i]
    the_ct <- tissue_celltype_count_max$ct[i]
    print(sprintf("%s -- %s" , the_tissue, the_ct))
    
    #Extract count table for this cell type  
    #tissue_count <- FACS_norm_arrays[[the_tissue]]   
    tissue_count <- count_norm_colmeans(sprintf("%s-counts.csv",the_tissue))
    input <- tissue_count[,colnames(tissue_count) %in% new_FACS_ontology$cell[new_FACS_ontology$collapsed==the_ct]]
    input <- log10(1+ input)
    
    #input <- input[1:100,]
    
    ### Reduce by PCA - later calculations become faster and much variation is removed. Should really rather use Facebook stochastic PCA
    print("  PCA")
    red_count <- na.omit(t(scale(t(input), scale = T, center = T)))
    #pca.norm = prcomp(red_count)
    pca.norm = prcomp_irlba(red_count, n = 6, retx = TRUE, center = TRUE)  #Lanczos-based PCA 
    
    
    
    calc_exp_knn <- function(k=10) {
      print("  umap")
      custom.settings = umap.defaults
      custom.settings$n_neighbors = k
      custom.settings$random_state = 666
      umap.norm <- umap(pca.norm$x, config = custom.settings)
      
      rownames(umap.norm$layout) <- rownames(red_count)
      write.csv(umap.norm$layout, file = sprintf("plots/data_umap_coexp_k%s_%s.csv", k, the_ct))
      
      ## Extract neighbour graph
      gene_names <- rownames(red_count)
      knn_named <- data.frame(
        from=rep(gene_names,ncol(umap.norm$knn$indexes)),
        to=gene_names[as.integer(umap.norm$knn$indexes)]
      )
      knn_named
    }  
    
    # Generate k-NN neighbour graph 
    calc_exp_knn_tsne <- function(k=20) {
      ### Reduce by tSNE 
      print("  tSNE")
      set.seed(666) #for reproducibility
      tsne.norm = Rtsne(pca.norm$x[,1:6], pca = FALSE, check_duplicates = F, perplexity = 100)
      nonlin_red <- tsne.norm$Y
      
      #separate kNN
      knn.norm = get.knn(as.matrix(nonlin_red), k = k)
      knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
                                       k), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
      knn_named <- data.frame(
        from=rownames(red_count)[knn.norm$from],
        to=rownames(red_count)[knn.norm$to]
      )
      knn_named
    }
    
    ### Generate feature: Co-expression neighbour #PMID  (for this cell type)
    k <- 10
    the_graph <- calc_exp_knn(k)
    feature_one_coexp <- calc_neigh_pmid(the_graph, feature_pmidcount[feature_pmidcount$ct==the_ct,], sprintf("coexp%d",k))
    
    ### Write one result at a time to handle potential crashes
    write.csv(
      feature_one_coexp[,c("gene","ct","coexp10")], 
      sprintf("features/feature_coexp_part/feature_coexp_%s",the_ct), row.names = FALSE)
  }
  
  
  ### Store feature, Co-expression neighbour #PMID
  feature_coexp <- NULL
  for(f in list.files("features","feature_coexp_part/feature_coexp_*")){
    feature_coexp <- rbind(
      feature_coexp,
      read.csv(sprintf("features/feature_coexp_part/%s",f)))
  }
  write.csv(
    feature_coexp[,c("gene","ct","coexp10")],
    "features/feature_coexp.csv", row.names = FALSE)
  nrow(feature_coexp)

  
  
  
  ## Store umap coordinates for website
  all_cell_types <- unique(feature_pmidcount$ct)
  all_umap_coord <- NULL
  for(celltype in all_cell_types){  #list.files("plots","data_umap_coexp_k10_*"
    the_file <- sprintf("plots/data_umap_coexp_k10_%s.csv",celltype)
    
    dat <- read.csv(the_file, stringsAsFactors = FALSE)#,sep=",")[,-1]#, col.names = FALSE)  
    colnames(dat) <- c("gene",sprintf("x_%s", celltype), sprintf("y_%s", celltype))
    if(is.null(all_umap_coord)){
      all_umap_coord <- dat
    } else {
      all_umap_coord <- merge(all_umap_coord, dat, all=TRUE)
    }
  }
  store_website_coordinates("coexp", all_umap_coord)

}




#############################################################
############## Features from chromatin structure ############
#############################################################

if(file.exists("features/feature_chromloc.csv")) {
  
  feature_chromloc <- read.csv("features/feature_chromloc.csv", stringsAsFactors = FALSE)
  
} else {

  
  
  ### Read chromosome locations of genes
  gtf <- read.table("input/Mus_musculus.GRCm38.97.gtf",sep="\t", stringsAsFactors=FALSE)
  gtf <- gtf[gtf$V3 =="gene",]
  extract_genename <- function(s) {
    gs <- str_extract(s,"gene_name [a-zA-Z0-9]+")
    if(is.na(gs)){
      return("")
    } else {
      return(str_sub(gs, 11))
    }
  }
  
  ### Compute gene midpoints
  feature_chromloc <- data.frame(
    chrom=gtf$V1,
    pos  = (gtf$V4 + gtf$V5)/2,
    gene = unname(sapply(gtf$V9, extract_genename)),
    stringsAsFactors = FALSE
  )
  feature_chromloc <- feature_chromloc[order(feature_chromloc$pos),]
  length(unique(feature_chromloc$gene))  #how can this be 53 000 ?????
  
  
  ###########################################
  ## Prepare coordinates for website
  coord_chrom <- merge(
    feature_chromloc,
    data.frame(
      chrom=c(1:19, "X","Y"),  #Only care about the "proper" chromosomes
      y=1:21
    ))
  coord_chrom <- data.frame(
    gene=coord_chrom$gene,
    x=coord_chrom$pos,
    y=coord_chrom$y)
  store_website_coordinates("chrom",coord_chrom)

  ###########################################
  ## Calculate avg(#PMID) of the k nearest neigbours along the chromosome.
  ## Assumes only one cell type
  calc_nearby_pmid <- function(dat, nearby_range=2){
    # For each gene, find the closest two genes. Since we only care about average positions,
    # this means we can simply consider a shorted list
    datsub <- dat[order(dat$pos),]
    datchrom_all <- NULL
    #Run for each chromosome
    for(curchrom in unique(datsub$chrom)){
      datchrom <- datsub[datsub$chrom==curchrom,]
      datchrom <- datchrom[order(datchrom$pos),]
      #For each gene
      datchrom$nearby_pmid <- NA
      for(i in 1:nrow(datchrom)){
        this_pmid  <- datchrom$rank_pmid[i]
        other_index <- c(max(1,i-nearby_range):min(i+nearby_range, nrow(datchrom)))
        other_index <- other_index[other_index!=i]
        #print(c(i,other_index))
        other_pmid <- datchrom$rank_pmid[other_index]
        datchrom$nearby_pmid[i] <- mean(other_pmid)
      }
      datchrom_all <- rbind(
        datchrom_all, 
        datchrom)
    }
    datchrom_all <- datchrom_all[!is.na(datchrom_all$nearby_pmid),]
    datchrom_all
  }
  
  #### For testing
  # datchrom <- calc_nearby_pmid(dat)  #for all
  # unique(dat$ct)
  # datchrom <- calc_nearby_pmid(dat[dat$ct=="T cell",])
  # plot(datchrom$rank_pmid, datchrom$nearby_pmid, 
  #      xlab="Rank #PMID", ylab="Nearby gene Rank #PMID")
  # cor(datchrom$rank_pmid, datchrom$nearby_pmid)#, method = "spearman")
  
  #### Link location - #PMID
  dat <- merge(feature_chromloc, feature_pmidcount)

  #### Calculate neighbour #PMID, for all cell types (need to do one ct at a time)
  datchrom_all <- NULL
  for(ct in unique(dat$ct)){
    print(ct)
    dat_for_ct <- dat[dat$ct==ct,]
    neigh_for_ct <- calc_nearby_pmid(dat_for_ct)
    datchrom_all <- rbind(
      datchrom_all,
      neigh_for_ct)  
  }

  
  #### Generate feature: Chromatin structure #PMID
  feature_chromloc <- datchrom_all[,c("gene","ct","nearby_pmid")]
  length(unique(feature_chromloc$gene))
  
  #### Store feature
  write.csv(
    feature_chromloc, 
    "features/feature_chromloc.csv", row.names = FALSE)

}


#### Quick plotting of correlation
# rs <- round(runif(10000, min=1, max=nrow(datchrom_all))) 
# rs <- which(datchrom_all$ct=="B cell")
# plot(datchrom_all$rank_pmid[rs], datchrom_all$nearby_pmid[rs],
#      xlab="RPC", ylab="Neigbour RPC")
# cor(datchrom_all$rank_pmid[rs], datchrom_all$nearby_pmid[rs])


####################### Plotting of expression & PMID along one chromosome


if(FALSE){

  dat <- merge(feature_chromloc, feature_pmidcount)
  dat <- merge(dat, feature_exp)
  length(unique(dat$gene))  #17000
  
  datsub <- dat[dat$ct=="T cell",]
  datsub <- datsub[datsub$chrom==10,]
  datsub <- datsub[order(datsub$pos),]
  
  thespan <- 0.1
  
  loessMod_pmid <- loess(rank_pmid ~ pos, data=datsub, span=thespan)
  smoothed_pmid <- predict(loessMod_pmid) / mean(predict(loessMod_pmid)) * 10000 
  
  loessMod_exp <- loess(rank_exp ~ pos, data=datsub, span=thespan)
  smoothed_exp <- predict(loessMod_exp) / mean(predict(loessMod_exp) ) * 10000
  
  plot(datsub$pos, datsub$rank_pmid,pch=19,cex=1, col="gray", ylim = c(0,15000),
       xlab="Chromosome position", ylab="Rank #PMID")
  #plot(datsub$pos, datsub$rank_exp)
  lines(datsub$pos, smoothed_pmid, col="red")
  lines(datsub$pos, smoothed_exp, col="green")
  
  # fig <- plot_ly(data = datsub, x = ~pos, y = ~rank_pmid,
  #                text = ~paste(gene))
  # fig

}


#############################################################
###### Features from gene homology graph neighbours #########
#############################################################

if(file.exists("features/feature_homology_pmid.csv")) {
  
  feature_homology_pmid <- read.csv("features/feature_homology_pmid.csv", stringsAsFactors = FALSE)
  
} else {

  ## Full graph; don't do this
  #graph_homology <- read.table("homology_graph.csv", sep="\t")[,1:2]#, col.names = FALSE)#[,c("from","to")]
  #colnames(graph_homology) <- c("from","to")
  
  ## Read reduced homology graph
  graph_homology <- read.csv("homology/homology_graph.red.csv")[,c("from","to")]
  
  ## Generate feature: homology neighbour #PMID
  feature_homology_pmid <- calc_neigh_pmid(graph_homology, feature_pmidcount, "homology_pmid")
  feature_homology_pmid <- feature_homology_pmid[,c("gene","ct","homology_pmid")]
  
  write.csv(
    feature_homology_pmid, 
    "features/feature_homology_pmid.csv", row.names = FALSE)
  
}




## Quick check of correlation
if(FALSE){
  foo <- feature_homology_pmid[feature_homology_pmid$ct=="fibroblast",]
  plot(foo$rank_pmid, foo$homology_pmid)
  cor(foo$rank_pmid, foo$homology_pmid)
  lm(homology_pmid~rank_pmid, foo)
  
}




#############################################################
############## Features from PPI neighbours #################
#############################################################

### Read HuRI data (ENSG*)
huri_raw <- data.frame(read.table("input/HuRI.tsv.txt",stringsAsFactors = F, sep= "\t"))
colnames(huri_raw) <- c("from","to")

### Create mapping ENSG* -> mouse symbol
map_human_mouse= read.csv("input/human_mouse.txt", stringsAsFactors = F)
colnames(map_human_mouse) <- c("ensmus","ensg")
map_ensmus_sym <- read.csv("input/map_ensmus_symbol.csv", stringsAsFactors = FALSE)
map_ensg_mousesym <- merge(map_ensmus_sym, map_human_mouse)

### Transform HuRI to mouse ENSMUS*, then to Mouse gene symbol
huri_symbol <- merge(merge(
  huri_raw,
  data.frame(
    from=map_ensg_mousesym$ensg,
    from_sym=map_ensg_mousesym$mouse_symbol
  )),
  data.frame(
    to=map_ensg_mousesym$ensg,
    to_sym=map_ensg_mousesym$mouse_symbol
  ))[,c("from_sym","to_sym")]
colnames(huri_symbol) <- c("from","to")

### Generate feature: PMID of PPI neighbours
feature_ppi <- calc_neigh_pmid(huri_symbol, feature_pmidcount, "ppi")  ################ TODO - normalize by the number of neighbours? *sqrt(n)?

feature_ppi <- feature_ppi[,c("gene","ct","ppi")]
write.csv(
  feature_ppi, 
  "features/feature_ppi.csv", row.names = FALSE)


### Quick check of correlation
if(FALSE){
  dat <- merge(feature_pmidcount, feature_ppi)
  dat2 <- dat[dat$ct=="T cell",]
  
  plot(dat2$rank_pmid, dat2$ppi)
  cor(dat2$rank_pmid, dat2$ppi)
  
}



#############################################################
########## Features for gene essentiality ###################
#############################################################


#### Read CRISPR screen data
crispr_data<- read.csv("input/crisprscreen.csv", header = T, stringsAsFactors = FALSE)[-1,]
colnames(crispr_data)[1:3]<-c("HGNC.symbol","perc_dependent_cells" ,"ensg")
crispr_data= crispr_data[!crispr_data$ensg=="",]
crispr_data$perc_dependent_cells <- as.double(crispr_data$perc_dependent_cells)

### Create mapping ENSG* -> mouse symbol
map_human_mouse= read.csv("input/human_mouse.txt", stringsAsFactors = F)
colnames(map_human_mouse) <- c("ensmus","ensg")
map_ensmus_sym <- read.csv("input/map_ensmus_symbol.csv", stringsAsFactors = FALSE)
map_ensg_mousesym <- merge(map_ensmus_sym, map_human_mouse)

#### Map to mouse ENSMUSG* and then to mouse symbol
crispr_data_working <- merge(
  crispr_data,
  map_ensg_mousesym)

mouse_human_ensid <- read.csv("input/human_mouse.txt",stringsAsFactors = F)

#### Keep only first occurence of duplicates, by hgnc symbol (dplyr)
crispr_data_working <- crispr_data_working %>% distinct(HGNC.symbol, .keep_all = T)

#### Extract feature
feature_essentiality <- data.frame(
  gene = crispr_data_working$mouse_symbol,
  perc_dependent_cells=crispr_data_working$perc_dependent_cells
)

### Store feature, gene essentiality
write.csv(
  feature_essentiality,   #note - no cell type
  "features/feature_essentiality.csv", row.names = FALSE)




#############################################################
##############  Feature for gene family index    ############
#############################################################


### Get merge first the celltype pmid and the g2phm_genespresent
pmid_term <- read.csv("features/keywords_pmids.csv",stringsAsFactors = FALSE)[,-1]
colnames(pmid_term) <- c("ct","pmid")
genes_ct_pmid <- merge(pmid_term, g2phm_genespresent)

### Get the ranked PMIDs for each cell type
feature_pmidcount <- read.csv("features/feature_ranked_pmid.csv", stringsAsFactors=FALSE)
feature_pmidcount$rank_pmid <- rank_by_ct(feature_pmidcount, "pmid_count")
feature_pmidcount <- within(feature_pmidcount, rm(pmid_count))

### Split gene names
uniques_g2phm_symbols<- unique(genes_ct_pmid$gene)
genefam<- data.frame(
  gene = uniques_g2phm_symbols,
  letters=str_extract(uniques_g2phm_symbols, "[a-zA-Z]*"), #letters until first number
  numbers=as.numeric(str_extract(uniques_g2phm_symbols, "[0-9]+")), #1st number
  stringsAsFactors = FALSE
)

# Keep genes which match the structure Xyz123
genefam <- genefam[genefam$gene==sprintf("%s%s", genefam$letters, genefam$numbers),]

### Get the founder - the gene with the smallest number in each family
genefam_founder <- merge(feature_pmidcount, genefam)
genefam_min_index <- sqldf("select distinct ct, letters, min(`numbers`) as numbers from genefam_founder group by letters, ct")
genefam_founder <- merge(genefam_founder, genefam_min_index)

gene_with_founder <- merge(
  merge(
    genefam,
    feature_pmidcount),
  data.frame(
    ct=genefam_founder$ct,
    letters=genefam_founder$letters,
    founder_rank_pmid=genefam_founder$rank_pmid,
    founder_index=genefam_founder$numbers
  )
)
gene_with_founder$number_diff <- gene_with_founder$numbers - gene_with_founder$founder_index

### Do not keep genes which have some ridiculous numbers
gene_with_founder <- gene_with_founder[gene_with_founder$numbers<2000,]

### Construct the feature
feature_founder_rank <- data.frame(
  gene=gene_with_founder$gene,
  ct=gene_with_founder$ct,
  family_indexdiff=gene_with_founder$number_diff,
  family_founder_rank=gene_with_founder$founder_rank_pmid
)

### Do not keep genes which are founders
feature_founder_rank <- feature_founder_rank[feature_founder_rank$family_indexdiff!=0,]

### Export feature
nrow(feature_founder_rank)
write.csv(feature_founder_rank, "features/feature_founder_fam_rankpmid.csv", row.names = F)


#############################################################
##############        First year feature          ###########
#############################################################

if(file.file.exists("features/feature_minyear_gene_ct.csv")){

  feature_minyear_gene_ct <- read.csv("features/feature_minyear_gene_ct.csv", stringsAsFactors = FALSE)
  
} else {
 
  ### Read when a paper was published and for which cell type
  pmid_term <- read.csv("features/keywords_pmids.csv",stringsAsFactors = F)[,-1] 
  colnames(pmid_term)<- c("ct", "pmid")
  
  pub_year <- read.csv("input/pubyear.csv", header=F, stringsAsFactors = F, sep="\t")[,1:2] 
  colnames(pub_year)<- c("pmid", "year")
  
  ### Merge tables
  pmids_present <- merge(g2phm_genespresent, pub_year)
  pmids_present <- merge(pmids_present, pmid_term)

  ### Figure out the first year for each gene and cell type
  feature_minyear_gene_ct <- sqldf("select gene, ct, min(year) as first_year from pmids_present group by gene, ct")
  
  ### Export feature
  write.csv(feature_minyear_gene_ct, "features/feature_minyear_gene_ct.csv", row.names = FALSE)

}





###############################################################################################################
########## Merge the features for the final model #############################################################
###############################################################################################################

####### Load all the features 

feature_pmidcount <- read.csv("features/feature_ranked_pmid.csv", stringsAsFactors=FALSE)
feature_pmidcount$rank_pmid <- rank_by_ct(feature_pmidcount, "pmid_count")
feature_pmidcount <- within(feature_pmidcount, rm(pmid_count))

feature_exp <- read.csv("features/feature_geneexp.csv", stringsAsFactors=FALSE)
feature_coexp <- read.csv("features/feature_coexp.csv", stringsAsFactors=FALSE)
feature_essentiality <- read.csv("features/feature_essentiality.csv", stringsAsFactors=FALSE)
feature_chromloc <- read.csv("features/feature_chromatin.csv", stringsAsFactors=FALSE)
feature_ppi <- read.csv("features/feature_ppi.csv", stringsAsFactors=FALSE)
#feature_genefam <- read.csv("features/founder_fam_rankpmid.csv", stringsAsFactors=FALSE)   ### Not yet generated here!   TODO
feature_minyear_gene_ct <- read.csv("features/feature_minyear_gene_ct.csv", stringsAsFactors = FALSE)
feature_founder_rank <- read.csv("features/feature_founder_fam_rankpmid.csv", stringsAsFactors = FALSE)


####### Merge all the features 

allfeat <- feature_pmidcount
allfeat <- merge(allfeat, feature_exp, all=TRUE)
allfeat <- merge(allfeat, feature_coexp, all=TRUE)
allfeat <- merge(allfeat, feature_essentiality, all=TRUE)
allfeat <- merge(allfeat, feature_chromloc, all=TRUE)
allfeat <- merge(allfeat, feature_ppi, all=TRUE)
allfeat <- merge(allfeat, feature_minyear_gene_ct, all=TRUE)
allfeat <- merge(allfeat, feature_founder_rank, all=TRUE)

allfeat <- allfeat[!is.na(allfeat$rank_pmid),]    ########## a surprising number of missing rank PMIDs. how can this be?

####### QC

colnames(allfeat)
print(length(unique(dat$gene)))
print(nrow(allfeat))

####### Write data for classification

write.csv(allfeat, "totfeature.csv", row.names = FALSE)

####### Write data for webserver

feature_long_name <- read.csv("input/feature_long_name.csv")

con <- dbConnect(SQLite(), dbname = "website/data/totfeature.sqlite")
dbWriteTable(con, "feature_matrix", allfeat, overwrite=TRUE)
dbWriteTable(con, "feature_desc", feature_long_name, overwrite=TRUE)
dbDisconnect(con)



############################## For debugging
# allfeat <- feature_pmidcount
# nrow(allfeat)
# allfeat <- merge(allfeat, feature_exp, all=FALSE)
# nrow(allfeat)
# allfeat <- merge(allfeat, feature_coexp, all=FALSE)
# nrow(allfeat)
# allfeat <- merge(allfeat, feature_essentiality, all=FALSE)
# nrow(allfeat)
# allfeat <- merge(allfeat, feature_chromloc, all=FALSE)
# nrow(allfeat)
# allfeat <- merge(allfeat, feature_ppi, all=FALSE)
# nrow(allfeat)
# allfeat <- merge(allfeat, feature_genefam, all=FALSE)
# nrow(allfeat)
# 
# 
