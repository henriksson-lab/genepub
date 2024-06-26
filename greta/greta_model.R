library(reshape)
library(greta)
library(splines)

cell_type <- "T cell"



###############################################
## Read all data
#totfeature <- read.csv("totfeature.csv", stringsAsFactors = FALSE)



###############################################
## Fit for every cell type
all_cell_type <- unique(totfeature$ct)
for(cell_type in all_cell_type) {
 #cell_type <- "T cell"
 
  allfeat <- read.csv(sprintf("greta/feature/%s.csv",cell_type))
  allfeat_meta <- allfeat[,c("gene","ct","orig_year")]
  allfeat <- allfeat[,!(colnames(allfeat_red) %in% c("gene","ct","orig_year"))]
  
  
  ###################################################################
  #### Construct polynomial design matrix
  get_poly_mat <- function(px, num_basepoints = 5){
    #px <- allfeat$rank_exp
    bpr <- range(px)
    bpr[1] <- bpr[1] - 0.001
    bpr[2] <- bpr[2] + 0.001
    
    ## Decide base points ... return?
    bp_x <- bpr[1] + (bpr[2]-bpr[1])*((1:num_basepoints)-1)/(num_basepoints-1)
    bp_dx <- bp_x[2]-bp_x[1]
    
    ## Which two base points are the closest?
    ind <- floor((px-bp_x[1])/bp_dx)+1
    ## Calculate weights
    w1 <- 1-(px-bp_x[ind])/bp_dx
    w2 <- 1- w1
    
    ## Generate values in sparse form
    npoint <- nrow(allfeat)
    pdata <- rbind(
      data.frame(
        row=1:npoint,
        col=ind,
        val=w1
      ),
      data.frame(
        row=1:npoint,
        col=ind+1,
        val=w2
      )
    )
    
    pmat <- cast(pdata, row~col, fill = 0, value="val")[,-1]
    list(bx = bp_x,
         design = pmat)
  }
  
  ###################################################################
  #### Construct polynomial design matrix
  get_spline_mat <- function(px, num_basepoints = 5){
    bpr <- range(px)
    
    ## Decide base points ... return?
    bp_x <- bpr[1] + (bpr[2]-bpr[1])*((1:num_basepoints)-1)/(num_basepoints-1)
    bp_dx <- bp_x[2]-bp_x[1]
    
    pmat <- bs(px, Boundary.knots = bpr, df = num_basepoints)
    
    #Design for plotting
    plot_px <- (bpr[2]-bpr[1])*(1:49)/50 + bpr[1]
    plot_pmat <- bs(plot_px, Boundary.knots = bpr, df = num_basepoints)
    
    list(bx = bp_x,
         design = pmat,
         plot_px = plot_px,
         plot_pmat = plot_pmat)
  }
  
  
  
  ###############################################
  ## Assemble all polynomials
  ## now mean(coefs) = mean(function) = 0    ... actually, this is not totally the same! 
  use_features <- colnames(allfeat)[colnames(allfeat)!="rank_pmid"]
  num_basepoints <- 5
  design <- matrix(nrow=nrow(allfeat), ncol=0)
  coefs <- normal(0, 10, dim = 0)
  all_poly <- list()
  for(i in 1:length(use_features)){
    the_feature <- use_features[i]
    #Get basepoints and design
    pmat <- get_spline_mat(allfeat[,the_feature], num_basepoints = num_basepoints)
    design <- cbind(design, pmat$design)
    all_poly[[the_feature]] <- pmat

    #Collect coefficients
    coefsf <- normal(0, 10, dim = num_basepoints-1)
    coefs <- rbind(coefs, c(coefsf, -sum(coefsf)))
  }
  
  
  
  
  
  #########################################
  ### Add other variables
  sd <- cauchy(0, 3, truncation = c(0, Inf))
  intercept <- normal(0, 10)
  mu <- intercept + design %*% coefs
  
  distribution(allfeat$rank_pmid) <- normal(mu, sd)
  
  
  #########################################
  ## Solve for the posterior
  m <- model(intercept, coefs, sd)
  draws <- mcmc(m, n_samples = 3000, chains = 8)  ### more later
  s <- summary(draws)
  
  
  #########################################
  ## Split up polynomial coefficients
  #poly_stat <- list()
  poly_columns <- list()
  pstat <- s$statistics[grep("coefs",rownames(s$statistics)),]
  colnames(design) <- rownames(pstat)
  for(i in 1:length(use_features)){
    the_feature <- use_features[i]
    
    ## Get coefficients for this polynomial
    thei <- (1+(i-1)*num_basepoints):(i*num_basepoints)
    all_poly[[the_feature]]$stat    <- pstat[thei,] 
    all_poly[[the_feature]]$statcol <- rownames(pstat)[thei]
  }
  
  
  #########################################
  ## Calculate contributions of features to each gene
  contribs <- matrix(nrow=nrow(allfeat), ncol=length(use_features))
  for(i in 1:length(use_features)) {
    the_feature <- use_features[i]
    the_columns <- all_poly[[the_feature]]$statcol   
    contribs[,i] <- design[,the_columns] %*% s$statistics[the_columns,"Mean"]
  }
  
  
  
  #########################################
  ## Store it all for later analysis
  the_fit <- list(
    statistics = s$statistics,
    all_poly = all_poly,
    allfeat = allfeat,
    allfeat_meta = allfeat_meta,
    intercept = s$statistics["intercept","Mean"],
    contribs = contribs,
    design = design
  )
  
  saveRDS(the_fit, file = sprintf("greta/out/%s",cell_type)) 

   
}


