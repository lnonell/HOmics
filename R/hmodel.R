#' Performs rjags models
#'
#' @param g.matrix matrix of cpgs, with samples in rows and cpgs in columns
#' @param z.matrix matrix of functional annotations for g.matrix cpgs, with cpgs in rows and functions in columns
#' @param cond response variable, as numerical vector
#' @param cont flag to consider a continuous response

#' @import rjags
#' @import MCMCvis
#' 
#' @return results for a gene model

#' @examples to be built

#' @export hmodel

hmodel <- function(g.matrix,z.matrix,cond,cont){
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
    dat <- list("n" = N, "y" = y,  "G" = G, "Z" = Z, "nG" = nG,"nZ" = nZ)
    
    jags.file <- ifelse(cont,
                        system.file("JAGSmodels", "hmodel.cont.txt", package="HOmics"),
                        system.file("JAGSmodels", "hmodel.binary.txt", package="HOmics"))
    # #depending on variable y, and if it has covariates, model should change
    # jags.file <- ifelse(cont,"D:/Doctorat/Hierarchical/Package/HOmics/JAGSmodels/hmodel.cont.txt",
    #                     "D:/Doctorat/Hierarchical/Package/HOmics/JAGSmodels/hmodel.binary.txt")
    
    jags.m <- jags.model( file = jags.file, data=dat, n.chains=3,  
                          n.adapt=1000 ) 
   
    # burn-in
    update(jags.m, n.iter = 2000) 
    
    # estimate
    samps <-coda.samples(jags.m, variable.names = c("beta"),n.iter=1000, thin=1) #thin has no effect on time

    # pvals
    sampmat<-as.matrix(samps)
    
    p.pos <- apply(sampmat, 2, function(x) round(mean(x > 0),4)) # Get p(beta > 0)
    p.neg <- apply(sampmat, 2, function(x) round(mean(x < 0),4)) # Get p(beta < 0)
    
    summ <- MCMCsummary(samps, Rhat=TRUE, n.eff=TRUE, round = 2)
   # summ <-paste(rownames(Z),rownames(summ),sep="-")
    summ <- as_tibble(cbind(summ,p.pos,p.neg))
    summ <- summ %>% mutate(cpg=rownames(Z)) %>% select(cpg,mean:p.neg)
    
    return(summ)
}
# ===============================================================================