#' Integrates two omic data through hierarchical modeling
#'
#' @param data.matrix matrix with features as rownames and samples as columns
#' @param agg.matrix matrix with colnames the features and rows the groups according to some feature aggregation criteria, 0 for non pertenance
#' @param cond response variable, usually a numerical factor with two levels representing the conditions to compare. If cond is a numerical vector (continuous response), a hiearchical linear regression model will be fit instead of the default hierarchical logistic regression model
#' @param z.matrix matrix with prior information related to features, with rownames the features and columns the samples
#' @param cont.covar.matrix vector or matrix of continuous covariates, with samples as rownames (in the same order as cond) and covariates as columns. Default = NULL
#' @param seed numerical seed for the use of function set.seed in the generation of the model, for reproducibility
#' @param cores cores in case of parallelization. Default = 1 (no parallelization)

#' @import dplyr
#' @import parallel
#' @import doParallel
#' @import foreach
#' 
#' @return an object of class HOmics
#'  
#' @examples to be built

#' @export HOmics


HOmics <- function(data.matrix, agg.matrix, cond, z.matrix, cont.covar.matrix = NULL, seed=NULL, cores=1, ...)
{
  
 # source(file="D:/Doctorat/Hierarchical/Package/HOmics/R/hmodel.R")
  
  call<-match.call() #to be returned at the end
  
  cont = FALSE
  
  ############################
  ###### Data relations ######
  ############################
  if (is.null(data.matrix) | !(is.matrix(data.matrix))) stop("data.matrix should be a matrix")
  if (is.null(agg.matrix) | !(is.matrix(agg.matrix))) stop("agg.matrix should be a matrix")
 
  #############################
  ####### Completeness  #######
  #############################
  if (anyNA(data.matrix)) {
    
    data.matrix<- data.matrix[complete.cases(data.matrix),]
    z.matrix <- z.matrix[rownames(data.matrix)]
    cat("some missing values were detected in data.matrix, only complete features will be selected\n" )
  } 
  
  if (!all(agg.matrix %in% c(0,1))) stop("agg.matrix should be a matrix of 0s and 1s, with 1 indicating association and 0 non association")
  if (!all(z.matrix %in% c(0,1))) stop("z.matrix should be a matrix of 0s and 1s, with 1 indicating prior information and 0 no prior information")
  
  #############################
  ##### Matrices matching #####
  #############################
  if (!all(rownames(z.matrix) %in% rownames(data.matrix))) stop("all rownames in z.matrix should match some rownames in data.matrix")
  if (!all(colnames(agg.matrix) %in% rownames(data.matrix))) stop("all colnames in agg.matrix should match some rownames in data.matrix")
  if (length(cond)!= ncol(data.matrix)) stop("cond vector length should match ncol and be in the same order")

  #############################
  ######### Condition #########
  #############################
  if (is.factor(cond)) {
    
    if (nlevels(cond) < 2) stop(paste0("cond is a factor but must contain at least 2 levels" )) 
    
    else if (nlevels(cond) > 2) {
      
      cat(paste0("cond is a factor and contains more than three levels, it has been converted to a numerical vector of 0s and 1s\n" ))
      cond <- as.numeric(cond) %% 2
      print(cond)
      
    } else {
      
      cat(paste0("cond is a factor, it has been converted to a numerical vector of 0s and 1s to fit a hierarchical logistic model\n" ))
      cond <- as.numeric(cond) %% 2
      print(cond)
      
    }
  } else if (is.character(cond)) {
    
    cat(paste0("cond is a character vector, it has been converted to a numerical vector of 0s and 1s\n" ))
    cond <- as.numeric(as.factor(cond)) %% 2
    print(cond)
    
  } else if (is.numeric(cond)) {
    
    cat(paste0("cond is a numerical vector, a hierarchical linear model will be constructed\n" ))
    cond <- cond
    cont <- TRUE
    
  } else stop ("Unrecognized type of condition")  
  
  #############################
  #######  Covariates  ########
  #############################
  #only continuous for the moment
  if(!is.null(cont.covar.matrix)){
    if (is.vector(cont.covar.matrix)){
      if (length(cont.covar.matrix) != length(cond) | !is.numeric(cont.covar.matrix))
       stop("cont.covar.matrix is a vector but must be continuous and have the same length as cond")
    } else if (nrow(cont.covar.matrix)!= length(cond) | !is.numeric(cont.covar.matrix)) 
      stop("cont.covar.matrix is a matrix but must be continuous and have the same number of rows as length of cond")
  cat(paste0("cont.covar.matrix will be included in the hierarchical model, as continuous covariates\n" ))
  } 
    
 
  #############################
  ############ Cores ##########
  #############################
  
  if (!is.numeric(cores)) stop("cores should be numeric")
  
  #############################
  ########### LOOP ############
  #############################
  
  N = nrow(agg.matrix)
  cl <- makeCluster(cores,type="PSOCK",outfile="output.txt")
  registerDoParallel(cl)
  
  # depending on cond and covariates a different function should be used 
  g.models<- foreach(i=1:N,.packages=c("rjags","tidyverse","MCMCvis"),.export="hmodel") %dopar% {
    #select feature belonging to g[i]
    gi.f <- colnames(agg.matrix)[as.logical(agg.matrix[i,])]
    g.matrix <- data.matrix[gi.f,] #un lio de terminologia amb la g.matrix que es la data matrix ara, canviar!
    z.mat <- z.matrix[gi.f,]
    hmodel(g.matrix = g.matrix, 
           z.matrix = z.mat,
           cond = cond, 
           cont = cont, 
           covar.matrix = cont.covar.matrix, 
           seed = seed)  
  }
  stopCluster(cl)
  
  names(g.models) <- rownames(agg.matrix)
  
  results<-g.models
 # results$call <- call
  class(results)<-"HOmics"
  return(results)  
}
# ===============================================================================