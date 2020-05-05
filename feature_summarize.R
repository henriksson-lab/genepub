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
library(igraph)
library(reshape)

#############################################################
## Function: Rank a column within each cell type (column ct)
rank_by_ct <- function(dat, col){
  for(ct in unique(dat$ct)){
    x <- dat[dat$ct==ct,col]

    #fp <- fitdist(x, "pareto")
    #px <- dinvpareto(x, shape = fp$estimate[1], scale = fp$estimate[2], log = FALSE)
    #dat[dat$ct==ct,col] <- 
    
    dat[dat$ct==ct,col] <- rank(dat[dat$ct==ct,col])

    
  }
  dat[,col]
}

if(FALSE){
  feature_pmidcount$rank_pmid <- rank_by_ct(feature_pmidcount, "pmid_count")
  hist(feature_pmidcount$rank_pmid[feature_pmidcount$ct=="T cell"])
  library(ParetoPosStable)

  x <- feature_pmidcount$pmid_count[feature_pmidcount$ct=="T cell"]
  hist((x^0.01)^0.01)
  hist((exp(x))^0.001)
  hist(rank(x, ties.method = "average")^4)
  fit <- PPS.fit(x)
  fit$estimate
  hist(x^fit$estimate$nu)
  plot(fit)
  plot(rank(x, ties.method = "average")^2)

  fp <- fitdist(x, "pareto")#, start=list(shape = 1, scale = 500))
  hist(dinvpareto(x, shape = fp$estimate[1], scale = fp$estimate[2], log = FALSE),breaks=100)

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


###############################################
## Function: Layout a graph (data frame), return coordinates
graph_layout_to_df <- function(graph_df, ct=NULL) {
  
  ## Do the layout
  igg <- graph_from_data_frame(graph_df, directed = FALSE, vertices = NULL)
  #lay <- layout_with_dh(igg, )  # too slow?
  #lay <- layout_with_kk(igg, )  # nice speed. not tested
  lay <- layout_with_fr(igg, )
  
  ## Figure out column names
  if(is.null(ct)){
    xname <- "x"
    yname <- "y"
  } else {
    xname <- sprintf("x_%s",ct)
    yname <- sprintf("y_%s",ct)
  }
  
  ## Construct data frame
  df <- data.frame(
    gene=vertex_attr(igg, "name")
  )
  df[,xname] <- lay[,1]
  df[,yname] <- lay[,2]
  
  df
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
  
  ### Not a full matrix - some papers missed. Fill in with 0s
  if(FALSE){
    tab <- cast(feature_pmidcount,gene~ct, fill = 0, value = "pmid_count")
    feature_pmidcount <- melt(tab)
    colnames(feature_pmidcount) <- c("gene","pmid_count","ct")
    rownames(feature_pmidcount) <- NULL
  }
  
  ### Do the ranking
  feature_pmidcount$rank_pmid <- rank_by_ct(feature_pmidcount, "pmid_count")
  
  ### Store  
  write.csv(feature_pmidcount, "features/feature_ranked_pmid.csv", row.names = F)
  
}


############# Calculate rank #paper

#feature_pmidcount <- within(feature_pmidcount, rm(pmid_count))






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
  for(f in list.files("features/feature_coexp_part","feature_coexp_*")){
    print(f)
    feature_coexp <- rbind(
      feature_coexp,
      read.csv(sprintf("features/feature_coexp_part/%s",f)))
  }
  write.csv(
    feature_coexp[,c("gene","ct","coexp10")],
    "features/feature_coexp.csv", row.names = FALSE)
  nrow(feature_coexp)
  
  
  
  
  ## Store umap coordinates for website
  all_cell_typess <- unique(feature_pmidcount$ct)
  all_umap_coord <- NULL
  for(celltype in all_cell_typess){  #list.files("plots","data_umap_coexp_k10_*"
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
    gs <- str_extract(s,"gene_name [a-zA-Z0-9\\-]+")
    if(is.na(gs)){
      return("")
    } else {
      return(str_sub(gs, 11))
    }
  }
  
  extract_ensid <- function(s) {
    gs <- str_extract(s,"gene_id [a-zA-Z0-9]+")
    if(is.na(gs)){
      return("")
    } else {
      return(str_sub(gs, 9))
    }
  }
  
  
  extract_gene_version <- function(s) {
    gs <- str_extract(s,"gene_version [a-zA-Z0-9]+")
    if(is.na(gs)){
      return("")
    } else {
      return(str_sub(gs, 14))
    }
  }
  
  
  ## Only keep these "proper" chromosomes
  consider_chrom <- c(1:19, "X","Y")
  
  ### Compute gene midpoints   -- mistake, bad name, reused later
  feature_chromloc <- data.frame(
    chrom=gtf$V1,
    pos  = (gtf$V4 + gtf$V5)/2,
    gene = unname(sapply(gtf$V9, extract_genename)),
    ensid = unname(sapply(gtf$V9, extract_ensid)),
    genever = as.double(unname(sapply(gtf$V9, extract_gene_version))),
    stringsAsFactors = FALSE
  )
  
  ### Filter out earlier gene versions
  feature_chromloc <- merge(feature_chromloc,sqldf("select gene, max(genever) as genever from feature_chromloc group by gene"))
  
  ### Filter out weird gene symbols and chromosomes
  feature_chromloc <- unique(feature_chromloc[feature_chromloc$gene %in% feature_pmidcount$gene & feature_chromloc$chrom %in% consider_chrom, ])
  
  #Crucial ordering for later
  feature_chromloc <- feature_chromloc[order(feature_chromloc$pos),]
  
  
  ###########################################
  ## Prepare coordinates for website
  coord_chrom <- merge(
    feature_chromloc,
    data.frame(
      chrom=consider_chrom, 
      y=1:length(consider_chrom)
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

  nrow((feature_chromloc))
  nrow(unique(feature_chromloc))
  
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
  
  
  ## Store graph for website use
  store_website_coordinates("homology",graph_layout_to_df(graph_homology))
  
  
  
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
length(unique(feature_ppi$gene))   ## only 6700

### Store HuRI layout for website
store_website_coordinates("PPI",graph_layout_to_df(huri_symbol))


### Quick check of correlation
if(FALSE){
  dat <- merge(feature_pmidcount, feature_ppi)
  dat2 <- dat[dat$ct=="T cell",]
  
  plot(dat2$rank_pmid, dat2$ppi)
  cor(dat2$rank_pmid, dat2$ppi)
  
}


#############################################################
########## Features from EBI GWAS ###########################
#############################################################

### Create a mapping human <-> mouse gene symbol
map_ensg_sym <- read.csv("input/map_ensg_symbol.csv", stringsAsFactors = FALSE)
map_human_mouse <- read.csv("input/human_mouse.txt", stringsAsFactors = F)
colnames(map_human_mouse) <- c("ensmus","ensg")
map_ensmus_sym <- read.csv("input/map_ensmus_symbol.csv", stringsAsFactors = FALSE)
map_human_mouse_sym <- unique(merge(merge(map_ensmus_sym, map_human_mouse),map_ensg_sym)[,c("mouse_symbol","human_symbol")])
map_human_mouse_sym <- map_human_mouse_sym[map_human_mouse_sym$mouse_symbol %in% feature_pmidcount$gene,]


## Read GWAS file
dat_gwas <- read.csv("input/gwas/gwas_catalog_v1.0-associations_e98_r2020-03-08.tsv",sep="\t", stringsAsFactors = FALSE, quote = "")
#dat_gwas <- read.csv("input/gwas/temp.csv",sep="\t", stringsAsFactors = FALSE)[1:500,]

## Figure out which mouse genes are annotated on each row
v <- str_split_fixed(dat_gwas$REPORTED.GENE.S.,pattern = ",", 10)
rownames(v) <- sprintf("row_%s",1:nrow(dat_gwas))
v <- melt(v)[,c(1,3)]
colnames(v) <- c("row","human_symbol")
map_row_sym <- unique(merge(v, map_human_mouse_sym)[,c("row","mouse_symbol")])

## Associate metadata to gene
dat_gwas_clean <- merge(
  map_row_sym,
  data.frame(
    row=sprintf("row_%s",1:nrow(dat_gwas)),
    pval=dat_gwas$P.VALUE,
    intergenic=dat_gwas$INTERGENIC
  )
)

## Ignore intergenic
dat_gwas_clean <- dat_gwas_clean[dat_gwas_clean$intergenic==0,]

## Aggregate p-values ... take smallest
dat_gwas_agg <- sqldf("select distinct mouse_symbol as gene, min(pval) as gwas_pval from dat_gwas_clean group by gene")

## How do we normalize? 
#hist(qnorm(dat_gwas_agg$gwas_pval),breaks = 100) 
#hist(rank(log10(dat_gwas_agg$gwas_pval)))

## Generate and store feaure
feature_gwas <- data.frame(
  gene=dat_gwas_agg$gene,
  rank_gwas=qnorm(dat_gwas_agg$gwas_pval)    #superior way. most important genes are negative
#  rank_gwas=rank(log10(dat_gwas_agg$gwas_pval))
)
feature_gwas <- feature_gwas[!is.na(feature_gwas$rank_gwas),]
feature_gwas$rank_gwas[is.infinite(feature_gwas$rank_gwas)] <- -40   ## fix, cannot go lower in probability

write.csv(
  feature_gwas, 
  "features/feature_gwas.csv", row.names = FALSE)


#############################################################
########## Features from COSMIC #############################
#############################################################

## Load and create a mapping ENST* to mouse symbol
map_ensg_transc <- read.csv("input/human_gene_transc.csv",sep="\t")[,1:2]
colnames(map_ensg_transc) <- c("ensg","enst")
map_human_mouse= read.csv("input/human_mouse.txt", stringsAsFactors = F)
colnames(map_human_mouse) <- c("ensmus","ensg")
map_ensmus_sym <- read.csv("input/map_ensmus_symbol.csv", stringsAsFactors = FALSE)
map_enst_mousesym <- unique(merge(map_ensg_transc,merge(map_ensmus_sym, map_human_mouse))[,c("enst","mouse_symbol")])
map_enst_mousesym <- map_enst_mousesym[map_enst_mousesym$mouse_symbol %in% feature_pmidcount$gene,]


## Load COSMIC mutation count and map to mouse gene symbols
dat_cosmic <- read.csv("input/cosmic/genecount.csv", stringsAsFactors = FALSE)
colnames(dat_cosmic) <- c("enst","count","length")
dat_cosmic <- merge(dat_cosmic, map_enst_mousesym)

## Reduce to genes
dat_cosmic_by_gene <- sqldf("select sum(count) as count, avg(length) as length, mouse_symbol from dat_cosmic group by mouse_symbol")

length(unique(dat_cosmic$mouse_symbol))
hist(log10(1+dat_cosmic$count/dat_cosmic$length))
hist(log10(1+dat_cosmic_by_gene$count/dat_cosmic_by_gene$length))

feature_cosmic <- data.frame(
  gene=dat_cosmic_by_gene$mouse_symbol,
#  rank_cosmic=log10(1+dat_cosmic_by_gene$count/dat_cosmic_by_gene$length)
  rank_cosmic=rank(dat_cosmic_by_gene$count/dat_cosmic_by_gene$length)   #better than log10
)
#hist(rank(dat_cosmic$count))
#hist(dat_cosmic$count)

write.csv(
  feature_cosmic, 
  "features/feature_cosmic.csv", row.names = FALSE)


#############################################################
########## Features for gene essentiality ###################
#############################################################


#### Read CRISPR screen data
crispr_data<- read.csv("input/crisprscreen.csv", header = T, stringsAsFactors = FALSE)[-1,]
colnames(crispr_data)[1:3]<-c("HGNC.symbol","essentiality_global" ,"ensg")
crispr_data <- crispr_data[!crispr_data$ensg=="",]
crispr_data$essentiality_global <- as.double(crispr_data$essentiality_global)

### Create mapping ENSG* -> mouse symbol
map_human_mouse= read.csv("input/human_mouse.txt", stringsAsFactors = F)
colnames(map_human_mouse) <- c("ensmus","ensg")
map_ensmus_sym <- read.csv("input/map_ensmus_symbol.csv", stringsAsFactors = FALSE)
map_ensg_mousesym <- merge(map_ensmus_sym, map_human_mouse)
map_ensg_mousesym <- map_ensg_mousesym[map_ensg_mousesym$mouse_symbol %in% feature_pmidcount$gene,]

#### Map to mouse ENSMUSG* and then to mouse symbol
crispr_data_working <- merge(
  crispr_data,
  map_ensg_mousesym)

#mouse_human_ensid <- read.csv("input/human_mouse.txt",stringsAsFactors = F)

#### Keep only first occurence of duplicates, by hgnc symbol (dplyr)
#crispr_data_working <- crispr_data_working %>% distinct(HGNC.symbol, .keep_all = T)   #### why do we need this?

#### Extract feature
feature_essentiality <- sqldf("select mouse_symbol as gene, avg(essentiality_global) as essentiality_global from crispr_data_working group by gene")

### Store feature, gene essentiality
write.csv(
  feature_essentiality,   #note - no cell type
  "features/feature_essentiality_global.csv", row.names = FALSE)


### Adding collapsed cell type terminology as proxy for cancer celltype and keeping feature as feature_crispr_ct_dependency
unique_celltype_mapping <- read.csv("input/unique_celltypes.csv", stringsAsFactors = F)
crispr_data_working_unique_celltype <- crispr_data_working[,-c(1,2,3,17)]
crispr_data_working_unique_celltype <- melt(crispr_data_working_unique_celltype, id= "mouse_symbol")
colnames(crispr_data_working_unique_celltype) <- c("gene","cancer_tissue", "essentiality_ct")

feature_essentiality_ct <- unique(merge(crispr_data_working_unique_celltype, unique_celltype_mapping, by= "cancer_tissue")[,c(2,3,4)])
feature_essentiality_ct$essentiality_ct <- as.double(feature_essentiality_ct$essentiality_ct)

feature_essentiality_ct <- sqldf("select gene, max(essentiality_ct) as essentiality_ct, ct from feature_essentiality_ct group by gene,ct")

### Store feature, gene essentiality for ct
write.csv(
  feature_essentiality_ct,
  "features/feature_essentiality_ct.csv", row.names = FALSE
)







#############################################################
##############  Feature for gene family index    ############
#############################################################


### Get merge first the celltype pmid and the g2phm_genespresent
pmid_term <- read.csv("features/keywords_pmids.csv",stringsAsFactors = FALSE)[,-1]
colnames(pmid_term) <- c("ct","pmid")
genes_ct_pmid <- merge(pmid_term, g2phm_genespresent)

### Get the ranked PMIDs for each cell type
feature_pmidcount <- read.csv("features/feature_ranked_pmid.csv", stringsAsFactors=FALSE)

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





###############################################
## Replace NA in each column with a neutral value
replace_feat_na <- function(allfeat, replace.func=median){
  for(i in 1:ncol(allfeat)){
    nalines <- is.na(allfeat[,i])
    allfeat[nalines,i] <- replace.func(allfeat[!nalines,i])#,na.rm = TRUE)
  }
  allfeat
}


###############################################
## Replace NA in one column with a neutral value
replace_feat_na_one <- function(allfeat, colname, replace.func=median){
  nalines <- is.na(allfeat[,colname])
  allfeat[nalines,colname] <- replace.func(allfeat[!nalines,colname])
  allfeat
}



####### Load all the features 

feature_pmidcount <- read.csv("features/feature_ranked_pmid.csv", stringsAsFactors=FALSE)
feature_pmidcount <- within(feature_pmidcount, rm(pmid_count))

feature_exp <- read.csv("features/feature_geneexp.csv", stringsAsFactors=FALSE)
feature_coexp <- read.csv("features/feature_coexp.csv", stringsAsFactors=FALSE)
feature_essentiality <- read.csv("features/feature_essentiality_global.csv", stringsAsFactors=FALSE)
feature_essentiality_ct<- read.csv("features/feature_essentiality_ct.csv", stringsAsFactors = FALSE)
feature_chromloc <- read.csv("features/feature_chromatin.csv", stringsAsFactors=FALSE)
feature_ppi <- read.csv("features/feature_ppi.csv", stringsAsFactors=FALSE)
feature_minyear_gene_ct <- read.csv("features/feature_minyear_gene_ct.csv", stringsAsFactors = FALSE)
feature_founder_rank <- read.csv("features/feature_founder_fam_rankpmid.csv", stringsAsFactors = FALSE)
feature_cosmic <- read.csv("features/feature_cosmic.csv", stringsAsFactors = FALSE)
feature_gwas <- read.csv("features/feature_gwas.csv", stringsAsFactors = FALSE)
feature_homology <- read.csv("features/feature_homology_pmid.csv", stringsAsFactors = FALSE)



####### Merge all the features 

totfeature <- feature_pmidcount
totfeature <- merge(totfeature, feature_exp, all=TRUE)
totfeature <- merge(totfeature, feature_coexp, all=TRUE)
totfeature <- merge(totfeature, feature_essentiality, all=TRUE)
totfeature <- merge(totfeature, feature_essentiality_ct, all=TRUE)
totfeature <- merge(totfeature, feature_chromloc, all=TRUE)
totfeature <- merge(totfeature, feature_ppi, all=TRUE)
totfeature <- merge(totfeature, feature_minyear_gene_ct, all=TRUE)
totfeature <- merge(totfeature, feature_founder_rank, all=TRUE)
totfeature <- merge(totfeature, feature_cosmic, all=TRUE)
totfeature <- merge(totfeature, feature_gwas, all=TRUE)
totfeature <- merge(totfeature, feature_homology, all=TRUE)


######## Deal with #PMID. 
#In many cases there are on papers about a gene, so need to set to lowest rank
#possible better is to set to 0 before ranking, for comparability
mean(is.na(totfeature$rank_pmid))
totfeature <- replace_feat_na_one(totfeature, "rank_pmid", min)
#older way - nasty! skews the first_year coefficient badly over time
#totfeature <- totfeature[!is.na(totfeature$rank_pmid),]   

sqldf("select ct, count(*) from totfeature group by ct")
### Too many. some cell types explode!


####### QC

colnames(totfeature)
print(nrow(totfeature))

#hist(na.omit(totfeature$ppi))

####### Write data for classification

write.csv(totfeature, "totfeature.csv", row.names = FALSE)

####### Write data for webserver

con <- dbConnect(SQLite(), dbname = "website/data/totfeature.sqlite")
dbWriteTable(con, "feature_matrix", totfeature, overwrite=TRUE)
feature_long_name <- read.csv("website/data/feature_long_name.csv")
dbWriteTable(con, "feature_desc", feature_long_name, overwrite=TRUE)
dbDisconnect(con)



###############################################
## Prepare separate dataset for each cell type, do imputation
#totfeature <- read.csv("totfeature.csv", stringsAsFactors = FALSE)
all_cell_types <- unique(totfeature$ct)
for(cell_type in all_cell_types) {
  #cell_type <- "T cell"
  #cell_type <- "pancreatic D cell"
  print(cell_type)
  
  ###### Choose one cell type to work with
  allfeat <- totfeature[totfeature$ct==cell_type,]
  allfeat_meta <- allfeat[,c("gene","ct")]
  allfeat <- allfeat[,!(colnames(allfeat) %in% c("gene","ct"))]
  #colMeans(is.na(allfeat_red))  #low is good
  nrow(allfeat)
  
  ###### Only consider cells with enough genes
  if(nrow(allfeat)>3000){
    
    ###### Fill in NAs for now
    allfeat <- replace_feat_na(allfeat)
    
    ###### Special treatment for the family index column: Cap it
    #allfeat$family_index[allfeat$family_index>10] <- 10
    allfeat$family_indexdiff[allfeat$family_indexdiff>10] <- 10
    
    ###### Special treatment for expression level: set to 0? no need
    #allfeat$family_indexdiff[allfeat$family_indexdiff>10] <- 10
    #sum(is.na(allfeat$rank_exp))
    #min(allfeat$rank_exp)
    
    ### Rescale; keep the year before scaling
    orig_year <- allfeat$first_year
    #allfeat_scaled <- as.data.frame(scale(allfeat,center = FALSE, scale = TRUE))
    allfeat_scaled <- as.data.frame(scale(allfeat))   #no longer center!
    
    allfeat_final <-  data.frame(
      allfeat_meta, 
      orig_year=orig_year, 
      allfeat_scaled)
    
    ###### Write similar object for python ML
    write.csv(
      allfeat_final,
      sprintf("greta/feature/%s.csv",cell_type), row.names = FALSE)    
  } else {
    print("   too few cells")
  }
}




###############################################
## Quality control of features

isfine_feature <- function(dat){
  length(unique(dat$gene))==length(dat$gene)  
}
isfine_feature(feature_essentiality)
isfine_feature(feature_cosmic)
isfine_feature(feature_gwas)

isfine_feature_ct <- function(dat){
  dat_u <- sqldf("select distinct gene, ct from dat")
  nrow(dat_u)==nrow(dat)  
}
isfine_feature_ct(feature_pmidcount)
isfine_feature_ct(feature_exp)
isfine_feature_ct(feature_coexp)  #TODO
isfine_feature_ct(feature_essentiality_ct)
isfine_feature_ct(feature_chromloc) 
isfine_feature_ct(feature_ppi)
isfine_feature_ct(feature_minyear_gene_ct)
isfine_feature_ct(feature_founder_rank)
isfine_feature_ct(feature_homology)




hist(totfeature$first_year)
hist(feature_minyear_gene_ct$first_year)

hist(totfeature$first_year[totfeature$ct=="T cell"])
hist(totfeature$first_year[totfeature$ct=="astrocyte"])


hist(allfeat_final$rank_pmid)

