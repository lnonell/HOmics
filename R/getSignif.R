#' Performs rjags models
#'
#' @param g.matrix matrix of cpgs, with samples in rows and cpgs in columns
#' @param z.matrix matrix of functional annotations for g.matrix cpgs, with cpgs in rows and functions in columns
#' @param cond response variable, as numerical vector

#' @import MCMCvis
#' @return hmodel

#' @examples to be built

#' @export hmodel

hmodel <- function(g.matrix,z.matrix,cond){
  #mirar on posar el model, si es pot posar fora o què
    ### G matrix
    G <- t(g.matrix)
    
    ### Z matrix
    Z <- z.matrix
    
    ### phenotype (condition)
    y <-cond

    N = nrow(G)
    nG = ncol(G)
    nZ = ncol(Z)
    dat <- list("n" = N, "y" = y,  "G" = G, "Z" = Z, "nG" = nG,"nZ" = nZ)
    
   # hmodel <- system.file("JAGSmodels", "hmodel.binary.txt", package="HOmics") 
    #depending on variable y, and if it has covariates, model should change
    jags.m <- jags.model( file = "D:/Doctorat/Hierarchical/Package/HOmics/JAGSmodels/hmodel.binary.txt", data=dat, n.chains=3,  
                          n.adapt=1000 ) 
   
    # burn-in
    update(jags.m, n.iter = 2000) 
    
    # estimate
    samps <-coda.samples(jags.m, variable.names = c("beta"),n.iter=1000, thin=1) #thin has no effect on time

    # pvals
    sampmat<-as.matrix(samps)
    
    p.pos <- apply(sampmat[, 1:nG], 2, function(x) round(mean(x > 0),4)) # Get p(beta > 0)
    p.neg <- apply(sampmat[, 1:nG], 2, function(x) round(mean(x < 0),4)) # Get p(beta < 0)
    
    summ <- MCMCsummary(samps, Rhat=TRUE, n.eff=TRUE, round = 2)
    
    return(cbind(summ,p.pos,p.neg))
}
# ===============================================================================