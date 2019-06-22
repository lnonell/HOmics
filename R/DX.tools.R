DX rjags model


dic.samples(jags.m,n.iter=1000)
# Mean deviance:  70.95 
# penalty 6.262 
# Penalized deviance: 77.21 

gelman.diag(samps, multivariate = F) #Rhat!!!

Point est. Upper C.I.
beta[1]        1.15       1.25
beta[2]        1.75       6.52
beta[3]        1.50       4.09
beta[4]        1.15       1.18
beta[5]        1.45       3.18
beta[6]        1.51       3.68
beta[7]        1.18       1.33
beta[8]        1.31       2.25
beta[9]        1.57       3.71
beta[10]       1.29       2.26
beta[11]       1.26       1.96
beta[12]       1.66       5.05
beta[13]       1.84       4.95
beta[14]       1.66       5.08
beta[15]       1.27       2.00
beta[16]       1.43       3.25
beta[17]       1.10       1.15
beta[18]       1.28       2.10
beta[19]       1.29       2.14
beta[20]       1.35       2.61
beta[21]       1.20       1.61
beta[22]       1.31       2.40
beta[23]       1.20       1.34
beta[24]       1.96       6.28
beta[25]       1.34       2.64
beta[26]       1.19       1.45
beta[27]       1.60       4.53
beta[28]       1.22       1.67
beta[29]       1.49       3.60

plot(samps,ask=T) #si cada color es una chain son molt variables ise'n van quantes mes iteracions

round(summary(samps)$quantiles,2)

library(MCMCvis)
MCMCsummary(samps, Rhat=TRUE, n.eff=TRUE, round = 2)
