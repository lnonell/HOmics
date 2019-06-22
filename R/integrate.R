#' Integrates transcritpomics and methylation
#'
#' @param meth.data GenomicRatioSet o multidataset o expressionset
#' @param pheno.cond.col response variable, column of annotation that contains condition (one of the variable names contained in pData of meth.data)
#' @param annot.gene.col column of annotation that contains gene names (same annotations as gene.list) 
#' @param annot.z.cols prior information. vector of  Z variables(for GenomicRatioSet functions getAnnotation(brge_methy) or granges(brge_methy) could be used). 
#' @param annot.mult.sep separartor in case annot.gene.col and maybe annot.z.cols contain multiple values such as in methylation data. Default=NULL
#' @param zvars prior information. Alternative parameter to annot.z.cols as matrix with rows the featureNames of meth.data and columns the Zvars can be given. Default=NULL
#' @param gene.list genes to analyze as a vector with symbols (HUGO)
#' @param cores cores in case of parallelization. Default=1

#' @return integrate

#' @import tidyverse
#' @import MultiDataSet
#' @import minfi
#' 
#' @return 
#' @examples to be built

#' @export integrate

#from the multiDataSet define the G and the Z for specific samples
#Z should be in the phenoData

integrate <- function(meth.data=GRSet,pheno.cond.col="tissue.ch1",annot.gene.col="UCSC_RefGene_Name",
                      annot.z.cols=c("UCSC_RefGene_Group"),annot.mult.sep=NULL,zvars=NULL,gene.list=genes.u,
                      cores=1,
                      ...)
{
 source(file="D:/Doctorat/Hierarchical/Package/HOmics/R/hmodel.R")
  #data should be in a MultiDataSet format!!! 
  #for the moment only methylation info is needed
  #annotations of methylation should contain gene info -> Z matrix
  #pheno should contain the cond vector
  #cond should be numeric, but can be converted to number
  #the annotations should contain Z column/s, if the separator string is specified i wll
  #assume that all categories are in the same field
  
  #########################
  ##Argument checking #####
  #########################
  call<-match.call() #to be returned at the end
  # de MEAL createmodel.R  
  # if (is(set, "eSet")){
  #   pFun <- Biobase::pData
  # }
  # else if (is(set, "SummarizedExperiment")){
  #   pFun <- SummarizedExperiment::colData
  # } else{
  #   stop("set must be an eSet or a SummarizedExperiment derived object")
  # }
  
  #matrices, the same for meth.data, tibble with id?
  if (!is(meth.data, "GenomicRatioSet")) stop("meth.data should be a GenomicRatioSet or an ExpressionSet or a SummarizedExperiment") else{
    pheno.data <- pData(meth.data)
 #   col.nm <- colnames(pheno.data) #to prevent from name changing
    pheno.data$id <- rownames(pheno.data)
    pheno.data <- as_tibble(pheno.data)
 #    colnames(pheno.data) <-col.nm
    
    annot.data <- getAnnotation(meth.data)
    annot.data$cpg <- rownames(annot.data)
    annot.data <- as_tibble(annot.data)
    
    cpg.data <- getBeta(meth.data)
  }
  
  if (is.null(zvars)) {
    if (is.null(pheno.cond.col) | !(pheno.cond.col %in% colnames(pheno.data))) stop("cond variable should be a variable in phenoData slot")
    #check levels in cond
    if (is.null(annot.gene.col) | !(annot.gene.col %in% colnames(annot.data))) stop("annot.gene.col should be a column in the annotations slot")
    if (is.null(annot.z.cols) | !(annot.z.cols %in% colnames(annot.data))) stop("annot.z.cols should be a list of columns in the annotations slot")
    
  } #else...comprovacions matriu Z
  
  #ARGCHECK: phenotype
  #variable selection: SOLVE NAMING
  cond.v <- pheno.data %>% select (!!as.name(pheno.cond.col)) %>% pull()
  if (is.factor(cond.v)) {
    if (nlevels(cond.v) < 2) stop(paste0("vector ", pheno.cond.col, " is a factor but must contain at least 2 levels" )) 
    else if (nlevels(cond.v) > 2) {
      warning (paste0("vector ", pheno.cond.col," is a factor and contains more than three levels, it has been converted to a numerical vector of 0s and 1s" ))  
      cond <- as.numeric(cond.v) %% 2
      cond
   } else {
      warning (paste0("vector ",pheno.cond.col," is a factor and has been converted to a numerical vector of 0s and 1s" ))
      cond <- as.numeric(cond.v) %% 2
      cond
    }
  } else if (is.character(cond.v)) {
    warning (paste0("vector ", pheno.cond.col," is a character vector, it has been converted to a numerical vector of 0s and 1s" ))  
    cond <- as.numeric(as.factor(cond.v)) %% 2
    cond
  } else if (is.numeric(cond.v)) {
    warning (paste0("vector ", pheno.cond.col," is a numerical vector, a hierarchical linear model will be constructed" ))  
    cond <- cond.v
  } else stop ("Unrecognized type of condition")  

  rm(cond.v)

  #covariates
  
  #ARGCHECK: cores numeric
  if (!is.numeric(cores)) stop("cores should be numeric")

  #########################
  ##NAs and imputation? ####
  #########################
  
  #########################
  ##### check annot #######
  #########################
 
  #########################
  #####    LOOP     #######
  #########################
  
N = length(gene.list)
cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)

# depending on cond and covariates a different function should be used 
cpg.models<- foreach(i=1:N,.packages=c("rjags","tidyverse","MCMCvis"),.export="hmodel") %dopar% {
    print(i)
    gen <- gene.list[i]
    #select annotations for that gene in specified col
    annot.gen <- annot.data %>% 
      filter(grepl(gen,!!as.name(annot.gene.col))) %>% 
      select(cpg,!!as.name(annot.gene.col),!!as.name(annot.z.cols))
    
    if(nrow(annot.gen)>0){ 
    
         #if (!is.null(annot.mult.sep))  # AKI: separate genes and z variables
          
        #number of genes max in a field
        nmax <- max(sapply(pull(annot.gen,!!as.name(annot.gene.col)),FUN = function(x) length(unlist(strsplit(x,split=annot.mult.sep)))))
        
        annot.gen.s <- annot.gen %>% separate(UCSC_RefGene_Name, paste("Gene",1:nmax,sep="_"),sep=annot.mult.sep) %>% 
          separate(UCSC_RefGene_Group, paste("GeneGroup",1:nmax,sep="_"),sep=annot.mult.sep)
        
        #create 2 tibbles, one for genes and one for z vars
        annot.gen.g <- annot.gen.s %>% select(cpg,Gene_1:paste("Gene",nmax,sep="_")) %>%
          gather("GeneVar","Gene",Gene_1:paste("Gene",nmax,sep="_"))
        annot.gen.z <- annot.gen.s %>% select(cpg,GeneGroup_1:paste("GeneGroup",nmax,sep="_")) %>%
          gather("GeneGroupVar","GeneGroup",GeneGroup_1:paste("GeneGroup",nmax,sep="_"))
        
        annot.gen <- inner_join(annot.gen.g,annot.gen.z,by=c("cpg"="cpg")) %>%
          filter(Gene==gen, !is.na(GeneGroup)) %>%
          select(cpg,Gene,GeneGroup) %>%
          distinct()
        #genegroup to dummies
        if(nrow(annot.gen)>0){
          z.n <- annot.gen   %>% select(-Gene) %>% 
            mutate(var = 1) %>% 
            spread(GeneGroup, var, fill = 0, sep = "_") 
          
          cpgs <- z.n %>% pull(cpg)
          #check if there are no cpgs
          z.matrix<-as.data.frame(select(z.n,-cpg)) #in fact it is a data.frame
          rownames(z.matrix) <- cpgs
          g.matrix <- cpg.data[cpgs,] #in fact it is a data.frame
          #muntar matriu, cpgs i pheno
      
          hmodel(g.matrix,z.matrix,cond) #function that performs model
      }  
    }
  }
stopCluster(cl)

names(cpg.models) <- gene.list[1:N]
return(cpg.models)

  # results$call <- call
  # results$set.def <- data_m
  # results$vars_df.def <- vars_df
  # class(results)<-"nlAssoc"
  # return(results)  
}
# ===============================================================================