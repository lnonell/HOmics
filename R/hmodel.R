#' Performs rjags models
#'
#' @param g.matrix matrix of cpgs, with samples in rows and cpgs in columns
#' @param z.matrix matrix of functional annotations for g.matrix cpgs, with cpgs in rows and functions in columns
#' @param cond response variable, as numerical vector
#' @param cont flag to consider a continuous response
#' @param covar.matrix matrix of continuous covariates (in columns)

#' @import rjags
#' @import MCMCvis
#' 
#' @return results for a gene model

#' @examples to be built

#' @export hmodel

hmodel <- function(g.matrix, z.matrix, cond, cont, covar.matrix=NULL, seed){
  #mirar on posar el model, si es pot posar fora o què
    ### G matrix
    if(is.matrix(g.matrix)) G<- t(g.matrix) else if (is.vector(g.matrix)) G <-as.matrix(g.matrix) else stop("error")
    
    ### Z matrix
    Z <- z.matrix
     
      ### phenotype (condition)
    y <-cond

    N = nrow(G)
    nG = ncol(G)
    nZ = ncol(Z)
   
    
    if(!is.null(seed)) set.seed(seed)
    
    if(is.null(covar.matrix)){
      jags.file <- ifelse(cont,
                          system.file("JAGSmodels", "hmodel.cont.txt", package="HOmics"),
                          system.file("JAGSmodels", "hmodel.binary.txt", package="HOmics"))
      dat <- list("n" = N, "y" = y,  "G" = G, "Z" = Z, "nG" = nG,"nZ" = nZ)
    } else {
      jags.file <- ifelse(cont,
                          system.file("JAGSmodels", "hmodel.cont.covariates.txt", package="HOmics"),
                          system.file("JAGSmodels", "hmodel.binary.covariates.txt", package="HOmics"))
      
      if(is.matrix(covar.matrix)) W<- covar.matrix else if (is.vector(covar.matrix)) W <-as.matrix(covar.matrix) else stop("error")
      nW <-ncol(W)
      
      dat <- list("n" = N, "y" = y,  "G" = G, "Z" = Z, "nG" = nG,"nZ" = nZ, W = W, nW =nW)
      
    }
    
    jags.m <- jags.model( file = jags.file, data = dat, n.chains = 3, n.adapt = 1000 ) 
   
    # burn-in
    update(jags.m, n.iter = 2000) 
    
    # estimate
    samps <-coda.samples(jags.m, variable.names = c("beta"), n.iter = 1000, thin = 1) 

    # pvals
    sampmat<-as.matrix(samps)
    
    p.pos <- apply(sampmat, 2, function(x) round(mean(x > 0), 4)) #es realment una frequencia dels pos i dels negs
    # apply(sampmat, 2, function(x) length(x[x < 0]))
    p.neg <- apply(sampmat, 2, function(x) round(mean(x < 0), 4))
    
    summ <- MCMCsummary(samps, Rhat = TRUE, n.eff = TRUE, round = 2)
 
    summ <- as_tibble(cbind(summ, p.pos, p.neg))
    summ <- summ %>% mutate(feature = rownames(Z)) %>% select(feature, mean:p.neg)
    
    return(summ)
}
# ===============================================================================