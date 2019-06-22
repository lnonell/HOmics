#19/6/19 model def amb RJags afegint de moment info posicional

#el trec fora per no haver-lo d'executar cada vegada, ja veuria on posar-lo
cat("model
        {
    for(i in 1:n){
    y[i] ~ dbern(q[i])
    logit(q[i]) <- beta0 + inprod(G[i,],beta[]) 
    }
    
    for(h in 1:nG){
    Ebeta[h] <- inprod(Z[h,],pi[])
    beta[h] ~ dnorm(Ebeta[h],prec.beta)
    }
    
    for(s in 1:nZ){
    pi[s] ~ dnorm(0,prec.pi)
    }
    # for(h in 1:p){
    #   beta[h] ~ dnorm(0,0.001)
    #  }
    beta0 ~ dnorm(0,0.001) #aquests valors quèeee
    
    prec.pi ~ dgamma(10,1)
    prec.beta ~ dgamma(1,1)
    }", file="Fnal.matrix.model.txt")

hmodel <- function(g.matrix,z.matrix,cond.v){
    ### G matrix
    G <- t(g.matrix)
    
    ### Z matrix
    Z <- z.matrix
    
    ### phenotype (condition)
    y <-as.numeric(as.factor(cond.v))
    #trick to convert to 2 categories: use module! (%%) 
    y <- y %% 2

    N = nrow(G)
    nG = ncol(G)
    nZ = ncol(Z)
    dat <- list("n" = N, "y" = y,  "G" = G, "Z" = Z, "nG" = nG,"nZ" = nZ) 
    jags.m <- jags.model( file = "Fnal.matrix.model.txt", data=dat, n.chains=3,  n.adapt=1000 ) #abans 500
    #SI, funciona !! a veure els resultats
    
    # burn-in
    update(jags.m, n.iter = 2000) #abans 10000
    summary(jags.m)
    
    # estimate
    samps <-coda.samples(jags.m, variable.names = c("beta"),n.iter=1000, thin=5) #abans 5000
    quants <-round(summary(samps)$quantiles,2)
    
    #pvals
    sampmat<-as.matrix(samps)
    
    apply(sampmat[, 1:nG], 2, function(x) mean(x > 0)) # Get p(beta > 0)
    apply(sampmat[, 1:nG], 2, function(x) mean(x < 0)) # Get p(beta < 0)

}