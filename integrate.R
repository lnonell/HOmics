#' Integrates transcritpomics and methylation
#'
#' @param meth.data GenomicRatioSet o multidataset o expressionset
#' @param gene.list	list of genes to study
#' @param cond response variable, one of the variable names contained in pData of meth.data
#' @param zvars	prior information. Either a vector of variable names contained in featureData (for GenomicRatioSet functions getAnnotation(brge_methy) or granges(brge_methy) could be used). Alternatively a matrix with rows the featureNames of meth.set and columns the Zvars can be given
#' @param cores cores in case of parallelization

#' @return integrate

#' @import tidyverse
#' @import MultiDataSet
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

# integrate <- function(expr.data,
#                       meth.data,
#                       pheno.data,
#                       
#                       gene.cpg.rel,
#                       hg.build,
#                       cores,
#                       ...)
# {
  
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
  
  #matrices, the same for meth.data, tibble with id?
  if (!is(meth.data, "GenomicRatioSet")) stop("meth.data should be a GenomicRatioSet or an ExpressionSet or a SummarizedExperiment") else{
    pheno.data <- pData(meth.data)  
    pheno.data$id <- rownames(pheno.data)
    pheno.data <- as_tibble(pheno.data)
    
    annot.data <- getAnnotation(meth.data)
    annot.data$cpg <- rownames(annot.data)
    annot.data <- as_tibble(annot.data)
    
    cpg.data <- getBeta(meth.data)
  }
  
  if (is.null(zvars)) {
    if (is.null(cond) | !(pheno.cond.col %in% colnames(pheno.data))) stop("cond variable should be a variable in phenoData slot")
    #check levels in cond
    if (is.null(annot.gene.col) | !(annot.gene.col %in% colnames(annot.data))) stop("annot.gene.col should be a column in the annotations slot")
    if (is.null(annot.z.cols) | !(annot.z.cols %in% colnames(annot.data))) stop("annot.z.cols should be a list of columns in the annotations slot")
    
  } #else...comprovacions matriu Z

  #ARGCHECK: cores numeric
  if (!is.numeric(cores)) stop("cores should be numeric")

  #########################
  ##NAs and imputation? ####
  #########################
  
  #########################
  ##### check annot #######
  #########################
  #start with loop
  
  #no es viable crear una super matriu d'anotacions amb tots els gens: cannot allocate vector
  #s'haurà de fer a cada pas del loop
 lg<-length(genes.u)
 lg=10
  
t1=Sys.time()
  for (i in 1:lg){
    print(i)
    gen <- genes.u[i]
    #select annotations for that gene in specified col
    annot.gen <- annot.data %>% 
                 filter(grepl(gen,!!as.name(annot.gene.col))) %>% 
                 select(cpg,!!as.name(annot.gene.col),!!as.name(annot.z.cols))
    
    if (!is.null(annot.mult.sep)) #separate genes and z variables
    
    #number of genes max in a field
    nmax <- max(sapply(pull(annot.gen,!!as.name(annot.gene.col)),FUN = function(x) length(unlist(strsplit(x,split=annot.mult.sep)))))
      
    annot.gen.s <- annot.gen %>% separate(UCSC_RefGene_Name, paste("Gene",1:nmax,sep="_"),sep=annot.mult.sep) %>% 
                  separate(UCSC_RefGene_Group, paste("GeneGroup",1:nmax,sep="_"),sep=annot.mult.sep)
    
    #create 2 tibbles, one for genes and one for 
    annot.gen.g <- annot.gen.s %>% select(cpg,Gene_1:paste("Gene",nmax,sep="_")) %>%
                  gather("GeneVar","Gene",Gene_1:paste("Gene",nmax,sep="_"))
    annot.gen.z <- annot.gen.s %>% select(cpg,GeneGroup_1:paste("GeneGroup",nmax,sep="_")) %>%
                  gather("GeneGroupVar","GeneGroup",GeneGroup_1:paste("GeneGroup",nmax,sep="_"))
      
    annot.gen <- inner_join(annot.gen.g,annot.gen.z,by=c("cpg"="cpg")) %>%
                 filter(Gene==gen, !is.na(GeneGroup)) %>%
                 select(cpg,Gene,GeneGroup) %>%
                 distinct()
     #genegroup to dummies
    z.n <- annot.gen   %>% select(-Gene) %>% 
                          mutate(var = 1) %>% 
                          spread(GeneGroup, var, fill = 0, sep = "_") 
   
    cpgs <- z.n %>% pull(cpg)
    z.matrix<-as.data.frame(select(z.n,-cpg)) #in fact it is a data.frame
    rownames(z.matrix) <-cpgs
    g.matrix <- cpg.data[cpgs,] #in fact it is a data.frame
    #muntar matriu, cpgs i pheno
    cond.v <- pheno.data %>% select (!!as.name(pheno.cond.col)) %>% pull()
    
    hmodel(g.matrix,z.matrix,cond.v) #function that performs model
}
t2<-Sys.time()
t2-t1 #Time difference of 8.127668 mins per fer-ne 10
# si trec el codi de jags fora Time difference of 8.021928 mins
#  2.17984 mins si ajusto n.iter y n.adapt -> 70h per als 21231 gens...mmmm
# buscar quins haurien de ser els valors y posar-ho potser com a param

#############################AMB FOREACH: MILLORA :COMBINAR ITER I foreach
#convertir a foreach!
library(doParallel)
library(foreach)
  
N = 10
cores=4 
cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
registerDoParallel(cl)
t1=Sys.time()
cpgs.models<- foreach(i=1:N,.combine=rbind, .packages=c("rjags","tidyverse")) %dopar% {
    print(i)
    gen <- genes.u[i]
    #select annotations for that gene in specified col
    annot.gen <- annot.data %>% 
      filter(grepl(gen,!!as.name(annot.gene.col))) %>% 
      select(cpg,!!as.name(annot.gene.col),!!as.name(annot.z.cols))
    
    if (!is.null(annot.mult.sep)) #separate genes and z variables
      
      #number of genes max in a field
     nmax <- max(sapply(pull(annot.gen,!!as.name(annot.gene.col)),FUN = function(x) length(unlist(strsplit(x,split=annot.mult.sep)))))
    
    annot.gen.s <- annot.gen %>% separate(UCSC_RefGene_Name, paste("Gene",1:nmax,sep="_"),sep=annot.mult.sep) %>% 
      separate(UCSC_RefGene_Group, paste("GeneGroup",1:nmax,sep="_"),sep=annot.mult.sep)
    
    #create 2 tibbles, one for genes and one for 
    annot.gen.g <- annot.gen.s %>% select(cpg,Gene_1:paste("Gene",nmax,sep="_")) %>%
      gather("GeneVar","Gene",Gene_1:paste("Gene",nmax,sep="_"))
    annot.gen.z <- annot.gen.s %>% select(cpg,GeneGroup_1:paste("GeneGroup",nmax,sep="_")) %>%
      gather("GeneGroupVar","GeneGroup",GeneGroup_1:paste("GeneGroup",nmax,sep="_"))
    
    annot.gen <- inner_join(annot.gen.g,annot.gen.z,by=c("cpg"="cpg")) %>%
      filter(Gene==gen, !is.na(GeneGroup)) %>%
      select(cpg,Gene,GeneGroup) %>%
      distinct()
    #genegroup to dummies
    z.n <- annot.gen   %>% select(-Gene) %>% 
      mutate(var = 1) %>% 
      spread(GeneGroup, var, fill = 0, sep = "_") 
    
    cpgs <- z.n %>% pull(cpg)
    z.matrix<-as.data.frame(select(z.n,-cpg)) #in fact it is a data.frame
    rownames(z.matrix) <-cpgs
    g.matrix <- cpg.data[cpgs,] #in fact it is a data.frame
    #muntar matriu, cpgs i pheno
    cond.v <- pheno.data %>% select (!!as.name(pheno.cond.col)) %>% pull()
    
    hmodel(g.matrix,z.matrix,cond.v) #function that performs model
    
    c("YES")
  }
stopCluster(cl)
t2<-Sys.time()
t2-t1 #Time dif 1.380336 mins 6 cores  2.157816 mins 1 core
# Time difference of 1.743935 mins (8  cores?? dec estar saturant)
# Time difference of 1.392021 mins 4 cores...
 
  # results$call <- call
  # results$set.def <- data_m
  # results$vars_df.def <- vars_df
  # class(results)<-"nlAssoc"
  # return(results)  

# ===============================================================================
