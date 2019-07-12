# HOmics

Hierarchical model to integrate omics data sets by accounting for prior knowledge.

It constructs a bayesian hierarchical model using the rjags interface to JAGS (Just Another Gibbs Sampler), a program for analysis of Bayesian hierarchical models using Markov Chain Monte Carlo (MCMC) simulation) .

The installation of JAGS is straightforward using this link [http://mcmc-jags.sourceforge.net/](http://mcmc-jags.sourceforge.net/).

With **HOmics** you perform a variety of different *omics* data integration. For instance, CpGs beta-values in a gene with prior information related to their gene position, or genes in a functional pathway with prior information on predicted genes of eQTL SNPs. Check vignettes for details an examples. This is an easy to work, user oriented available R package, adapted to common omic R structures including ExpressionSet and GenomicRaioSet or tibble.





