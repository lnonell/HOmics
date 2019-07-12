#' Integrates methylation and gene data
#'
#' @param meth.data GenomicRatioSet or  ExpressionSet with methylation data
#' @param pheno.cond.col response variable, column of annotation that contains condition (one of the variable names contained in pheno data obtained with function pData())
#' @param annot.gene.col column of annotation that contains gene names (same annotations as gene.list). Default = "UCSC_RefGene_Name"
#' @param annot.z.col prior information. Column of annotation that contains prior variable (for GenomicRatioSet functions getAnnotation() or fData() for ExpressionSet). Default = "UCSC_RefGene_Group"
#' @param annot.mult.sep separator for annot.gene.col and maybe annot.z.col multiple values. Default = ";"
#' @param z.matrix prior information. Alternative parameter to annot.z.col as matrix with rows the featureNames of meth.data and columns the Zvars can be given. Default = NULL
#' @param pheno.covar.col covariates, column of annotation that contains names of covariates to include in the model (one or more of the variable names contained in pheno data). Default = NULL
#' @param gene.list genes to analyze as a vector with symbols (HUGO)
#' @param seed numerical seed for the use of function set.seed in the generation of the model, for reproducibility
#' @param cores cores in case of parallelization. Default = 1 (no parallelization)

#' @import dplyr
#' @import Biobase
#' @import minfi
#' @import parallel
#' @import doParallel
#' @import foreach
#' 
#' @return an object of class HOmics
#'  
#' @examples to be built

#' @export HOmics.meth


HOmics.meth <- function(meth.data = GRSet, pheno.cond.col="tissue.ch1", annot.gene.col = "UCSC_RefGene_Name",
                      annot.z.col = c("UCSC_RefGene_Group"), annot.mult.sep = ";", z.matrix = NULL, pheno.covar.col = NULL, 
                      gene.list = genes.u, seed=NULL, cores = 1,
                      ...)
{
  
  cont = FALSE
  
  #########################
  ##Argument checking #####
  #########################
  call<-match.call() #to be returned at the end
  
  if (is(meth.data, "GenomicRatioSet")){
    pheno.data <- pData(meth.data)
    pheno.data$id <- rownames(pheno.data)
    pheno.data <- as_tibble(pheno.data)

    annot.data <- getAnnotation(meth.data)
    annot.data$cpg <- rownames(annot.data)
    annot.data <- as_tibble(annot.data)
    
    cpg.data <- getBeta(meth.data)
    
  } else if (is(meth.data, "ExpressionSet")) {
    
    pheno.data <- pData(meth.data)
    pheno.data$id <- rownames(pheno.data)
    pheno.data <- as_tibble(pheno.data)
    
    annot.data <- fData(meth.data)
    annot.data$cpg <- rownames(annot.data) 
    annot.data <- as_tibble(annot.data)
    
    cpg.data <- exprs(meth.data)
    
  } else stop("meth.data should be a GenomicRatioSet or an ExpressionSet or a SummarizedExperiment") 
  
  if (is.null(pheno.cond.col) | !(pheno.cond.col %in% colnames(pheno.data))) stop("cond variable should be a variable in phenoData slot")
 
  if (is.null(annot.gene.col) | !(annot.gene.col %in% colnames(annot.data))) stop("annot.gene.col should be a column in the annotations slot")
  
  #########################
  ## Completeness cpgs ####
  #########################
  if (anyNA(cpg.data)) {
    cpg.data <-  cpg.data[complete.cases(cpg.data),]
    annot.data <- annot.data %>% filter(cpg %in% rownames(cpg.data))
    cat("some missing values were detected, only complete features will be selected\n" )
  } 
  
  if (is.null(z.matrix)) {
    if (is.null(annot.z.col) | !(annot.z.col %in% colnames(annot.data))) stop("annot.z.col should be a column in the annotations slot of meth.data")
    
  } else {
    
    #check common cpgs
    cpgs.c<-intersect(rownames(cpg.data),rownames(z.matrix))
    
    if (length(cpgs.c)==0) stop ("z.matrix must contain rownames matching featureNames in meth.data")
    
    else
      
      cat("z.matrix will be used to construct the hierarchical model\n")
      z.matrix <- z.matrix[cpgs.c,]
      cpg.data <- cpg.data[cpgs.c,]
      annot.data <- annot.data %>% filter(cpg %in% cpgs.c)
      
  }
    
  cond.v <- pheno.data %>% select (!!as.name(pheno.cond.col)) %>% pull()
  
  if (is.factor(cond.v)) {
    
    if (nlevels(cond.v) != 2 ) {
      
      stop(paste0("vector ", pheno.cond.col, " is a factor but must have 2 levels" )) 
      
    } else {
      cond.min <- levels(cond.v)[1]
      cat(paste0(pheno.cond.col," is a factor, it has been converted to a numerical vector of 0s and 1s with ",cond.min," as the reference level. A hierarchical logistic model will be fitted\n" ))
      cond <- as.numeric(cond.v)-1
      print(cond)
      
    }
  } else if (is.character(cond.v)) {
    cond.min <- min(cond.v)
    cat(paste0(pheno.cond.col," is a factor, it has been converted to a numerical vector of 0s and 1s with ",cond.min," as the reference level. A hierarchical logistic model will be fitted\n" ))
    cond <- as.numeric(as.factor(cond.v)) -1
    print(cond)
    
  } else if (is.numeric(cond.v)) {
    
    cat(paste0(pheno.cond.col," is a numerical vector, a hierarchical linear model will be constructed\n" ))
    cond <- cond.v
    cont <- TRUE
    
  } else stop ("Unrecognized type of condition")  
 

  rm(cond.v)
  #########################
  ## Covariates  ####
  #########################
  #flag covar
  if (!is.null(pheno.covar.col)){
    
    if (!(pheno.covar.col %in% colnames(pheno.data))) stop("covariate ", pheno.covar.col, " should be a variable in phenoData slot") 
    else   covar.v <- pheno.data %>% select (!!as.name(pheno.covar.col)) %>% pull()
    
    if(is.factor(covar.v) | is.character(covar.v)) {
      
      covar.v <- as.numeric(as.factor(covar.v)) -1
      cat(paste0("covariate ", pheno.covar.col, " has been converted to a numerical vector\n" ))
      
    } else cat(paste0("covariate ", pheno.covar.col, " will be included in the hierarchical model, as a continuous covariate\n")) 
    
  } else covar.v=NULL
    
  
  #########################
  #######  Cores  #########
  #########################
  
  if (!is.numeric(cores)) stop("cores should be numeric")

  #########################
  #####    LOOP     #######
  #########################
  
N = length(gene.list)
cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)


cpg.models<- foreach(i=1:N,.packages=c("rjags","tidyverse","MCMCvis"),.export="hmodel") %dopar% {
    gen <- gene.list[i]
    #select annotations for that gene in specified col
    if (is.null(z.matrix)){
      annot.gen <- annot.data %>% 
                  filter(grepl(gen,!!as.name(annot.gene.col))) %>% 
                  select(cpg,!!as.name(annot.gene.col),!!as.name(annot.z.col))
      
      if(nrow(annot.gen)>0){ 
        #number of genes max in a field
        nmax <- max(sapply(pull(annot.gen,!!as.name(annot.gene.col)),
                           FUN = function(x) length(unlist(strsplit(x,split=annot.mult.sep)))))
        
       #create 2 tibbles, one for genes and one for z vars
       annot.gen.g <- annot.gen %>% 
                      separate(!!as.name(annot.gene.col), paste("Gene",1:nmax,sep="_"),sep=annot.mult.sep) %>% 
                      select(cpg,Gene_1:paste("Gene",nmax,sep="_")) %>%
                      gather("GeneVar","Gene",Gene_1:paste("Gene",nmax,sep="_"))
        
        annot.gen.z <- annot.gen %>% 
                      separate(!!as.name(annot.z.col), paste("GeneGroup",1:nmax,sep="_"),sep=annot.mult.sep) %>%
                      select(cpg,GeneGroup_1:paste("GeneGroup",nmax,sep="_")) %>%
                      gather("GeneGroupVar","GeneGroup",GeneGroup_1:paste("GeneGroup",nmax,sep="_"))
        
        annot.gen <- inner_join(annot.gen.g,annot.gen.z,by=c("cpg"="cpg")) %>%
                      filter(Gene==gen, !is.na(GeneGroup)) %>%
                      select(cpg,Gene,GeneGroup) %>%
                      distinct()
        
        #genegroup to dummies
        if(nrow(annot.gen)>0){
          z.n <- annot.gen   %>% 
                select(-Gene) %>% 
                mutate(var = 1) %>% 
                spread(GeneGroup, var, fill = 0, sep = "_") 
          
          cpgs <- z.n %>% pull(cpg)
          
          z.mat<-as.data.frame(select(z.n,-cpg)) 
          rownames(z.mat) <- cpgs
          g.matrix <- cpg.data[cpgs,] 
          hmodel(g.matrix = g.matrix, 
                 z.matrix = z.mat,
                 cond = cond, 
                 cont = cont, 
                 covar.matrix = covar.v, 
                 seed = seed)  
        }
      }  
    } else {
      annot.gen <- annot.data %>% 
                    filter(grepl(gen,!!as.name(annot.gene.col))) %>% 
                    select(cpg,!!as.name(annot.gene.col))
      if(nrow(annot.gen)>0){ 
        #number of genes max in a field
        nmax <- max(sapply(pull(annot.gen,!!as.name(annot.gene.col)),
                           FUN = function(x) length(unlist(strsplit(x,split=annot.mult.sep)))))
        annot.gen.g <- annot.gen %>% 
                      separate(!!as.name(annot.gene.col), paste("Gene",1:nmax,sep="_"),sep=annot.mult.sep) %>% 
                      select(cpg,Gene_1:paste("Gene",nmax,sep="_")) %>%
                      gather("GeneVar","Gene",Gene_1:paste("Gene",nmax,sep="_")) %>%
                      filter(Gene==gen) %>% 
                      select(cpg,Gene) %>%
                      distinct()
        if(nrow(annot.gen.g)>0){
          cpgs <-  annot.gen.g %>% pull(cpg)
          z.mat <- z.matrix[cpgs,]
          if (is.vector(z.mat)) z.mat <-as.matrix(t(z.mat)) 
          rownames(z.mat) <- cpgs
          g.matrix <- cpg.data[cpgs,] 
          hmodel(g.matrix,z.mat,cond,cont=cont) #function that performs model
        }  
      }  
    }
 }
stopCluster(cl)

names(cpg.models) <- gene.list[1:N]

results<-cpg.models
#results$call <- call
class(results) <- "HOmics"
return(results)  
}
# ===============================================================================