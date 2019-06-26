#' Integrates transcritpomics and methylation
#'
#' @param data.matrix matrix with rownames the features and columns the samples
#' @param agg.matrix matrix with colnames the features and columns the aggregation criteria, 0 for non pertenance
#' @param cond response variable, usually a numerical factor with two levels representing the conditions to compare. If cond is a numerical vector (continuous response), a hiearchical linear regression model will be fit instead of the default hierarchical logistic regression model
#' @param z.matrix column of annotation that contains gene names (same annotations as gene.list) 
#' @param cores cores in case of parallelization. Default=1

#' @import parallel
#' @import foreach
#' 
#' @return an object of class HOmics
#'  
#' @examples to be built

#' @export HOmics


HOmics <- function(data.matrix,agg.matrix,cond,z.matrix,cores=1,
                      ...)
{
  
 # source(file="D:/Doctorat/Hierarchical/Package/HOmics/R/hmodel.R")
  
  call<-match.call() #to be returned at the end
  
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
  ############ Cores ##########
  #############################
  
  if (!is.numeric(cores)) stop("cores should be numeric")
  
  cont = FALSE
  covar = FALSE
  
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
    hmodel(g.matrix,z.mat,cond,cont=cont)  
  }
  stopCluster(cl)
  
  names(g.models) <- rownames(agg.matrix)
  
  results<-g.models
 # results$call <- call
  class(results)<-"HOmics"
  return(results)  
}
# ===============================================================================