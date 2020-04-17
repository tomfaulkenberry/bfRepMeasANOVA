


##########################
# functions for simulation 

simData = function(N,k,ES,rho){
  # note: for all sims, keep sigmaP^2+sigmaE^2=100000. This leads to following:
  # if rho = 0.2, set sigmaP^2 = 20,000 --> sigmaP=141 and sigmaE = 283
  # if rho = 0.8, set sigmaP = 283 and sigmaE = 141
  sigmaP = sqrt(rho*1e5)
  sigmaE = sqrt((1-rho)*1e5) 
  diff = ES*sigmaP/sqrt(rho)
  mid = runif(n=1, min=0, max=1)
  mu = c(1000, 1000+diff*mid, 1000+diff)
  
  p = rnorm(n=N, mean=0, sd=sigmaP)
  epsilon = matrix(rnorm(n=N*k, mean=0, sd=sigmaE), nrow=N, ncol=k)
  
  data = matrix(0, nrow=N, ncol=k)
  for (i in 1:N){
    for (j in 1:k){
      data[i,j] = mu[j] + p[i] + epsilon[i,j]
    }
  }
  return(data)
}

bf = function(F, n, k){
  return(sqrt((n*k-n)^(k-1)*(1+F/(n-1))^(n-n*k)))
}

anovas = function(data, N, k){
  # BIC computations
  tmtMeans = apply(data, 2, mean)
  subjMeans = apply(data, 1, mean)
  
  SStot = sum(data^2)-sum(data)^2/length(data)
  SSsubj = length(tmtMeans)*sum((subjMeans-rep(mean(data), length(subjMeans)))^2)
  SStmt = length(subjMeans)*sum((tmtMeans-rep(mean(data), length(tmtMeans)))^2)
  SSresid = SStot-SSsubj-SStmt
  
  F = (SStmt/SSresid)*(N-1)

  #BForiginal = exp(0.5*(N*(k-1)*log((SStot-SStmt-SSsubj)/(SStot-SSsubj))+(k-1)*log(N*k-N)))
  BForiginal = bf(F, n=N, k=k)
  BFimproved = exp(0.5*(N*(k-1)*log((SStot-SStmt-SSsubj)/(SStot-SSsubj))+(k+2)*log((N*(SStot-SStmt))/SSsubj)-3*log(N*SStot/SSsubj)))
  
  return(list(PMPoriginal = BForiginal/(1+BForiginal), # original PMP from Masson 2011
              PMPimproved = BFimproved/(1+BFimproved)  # "improved" PMP from Nathoo & Masson (2016)
  ))
}


# do simulations

# make wrapper function to generate sims into CSV
doSim = function(ES,rho,numSims=100){
  # set up empty matrices for model probs
  pNullOriginal = matrix(0, nrow=numSims, ncol=3)
  pNullImproved = matrix(0, nrow=numSims, ncol=3)

  for (i in 1:numSims){
    # N=20
    data = simData(N=20, k=3, ES, rho)
    results = anovas(data, N=20, k=3)
    pNullOriginal[i,1] = results$PMPoriginal
    pNullImproved[i,1] = results$PMPimproved

    # N=50
    data = simData(N=50, k=3, ES, rho)
    results = anovas(data, N=50, k=3)
    pNullOriginal[i,2] = results$PMPoriginal
    pNullImproved[i,2] = results$PMPimproved

    # N=80
    data = simData(N=80, k=3, ES, rho)
    results = anovas(data, N=80, k=3)
    pNullOriginal[i,3] = results$PMPoriginal
    pNullImproved[i,3] = results$PMPimproved
  }

  # export to CSV
  probs = cbind(pNullOriginal, pNullImproved)
  setwd("~/Dropbox/papers/2019/repMeasANOVAcalc/sim")
  filename = sprintf("sim-%.1f-%.1f.csv",ES,rho)
  write.csv(probs, file=filename, row.names=FALSE)
}

# do simulations

doSim(ES=0, rho=0.2, numSims = 1000)
doSim(ES=0.2, rho=0.2, numSims = 1000)
doSim(ES=0.5, rho=0.2, numSims = 1000)
doSim(ES=0, rho=0.8, numSims = 1000)
doSim(ES=0.2, rho=0.8, numSims = 1000)
doSim(ES=0.5, rho=0.8, numSims = 1000)



# make boxplots of posterior model probabilities

sims02 = read.csv("sim-0.0-0.2.csv")
sims22 = read.csv("sim-0.2-0.2.csv")
sims52 = read.csv("sim-0.5-0.2.csv")
sims08 = read.csv("sim-0.0-0.8.csv")
sims28 = read.csv("sim-0.2-0.8.csv")
sims58 = read.csv("sim-0.5-0.8.csv")


# set plot layout
layout(mat = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2))

# change arrangement for plotting (requested by both reviewers)
sims02r = cbind(sims02[,1], sims02[,4], sims02[,2], sims02[,5], sims02[,3], sims02[,6])
sims22r = cbind(sims22[,1], sims22[,4], sims22[,2], sims22[,5], sims22[,3], sims22[,6])
sims52r = cbind(sims52[,1], sims52[,4], sims52[,2], sims52[,5], sims52[,3], sims52[,6])
sims08r = cbind(sims08[,1], sims08[,4], sims08[,2], sims08[,5], sims08[,3], sims08[,6])
sims28r = cbind(sims28[,1], sims28[,4], sims28[,2], sims28[,5], sims28[,3], sims28[,6])
sims58r = cbind(sims58[,1], sims58[,4], sims58[,2], sims58[,5], sims58[,3], sims58[,6])

boxplot(sims02r,
        ylim=c(0,1),
        outline=FALSE,
        col = rep(c("white","dark gray"),3),
        staplewex = 0.2,
        names = c(20,20,50,50,80,80),
        xlab = "Number of observations",
        ylab = "Posterior probability of H0",
        main = "Null effect with correlation = 0.2")

boxplot(sims22r,
        ylim=c(0,1),
        outline=FALSE,
        col = rep(c("white","dark gray"),3),
        staplewex = 0.2,
        names = c(20,20,50,50,80,80),
        xlab = "Number of observations",
        ylab = "Posterior probability of H0",
        main = "Small effect with correlation = 0.2")

boxplot(sims52r,
        ylim=c(0,1),
        outline=FALSE,
        col = rep(c("white","dark gray"),3),
        staplewex = 0.2,
        names = c(20,20,50,50,80,80),
        xlab = "Number of observations",
        ylab = "Posterior probability of H0",
        main = "Medium effect with correlation = 0.2")

boxplot(sims08r,
        ylim=c(0,1),
        outline=FALSE,
        col = rep(c("white","dark gray"),3),
        staplewex = 0.2,
        names = c(20,20,50,50,80,80),
        xlab = "Number of observations",
        ylab = "Posterior probability of H0",
        main = "Null effect with correlation = 0.8")

rect(3,0.3,3.5,0.4, col="white")
text(3.5,0.35, "Minimal BIC method", pos=4, cex=0.8)
rect(3,0.1,3.5,0.2, col="dark gray")
text(3.5,0.15, "Nathoo & Masson method", pos=4, cex=0.8)


boxplot(sims28r,
        ylim=c(0,1),
        outline=FALSE,
        col = rep(c("white","dark gray"),3),
        staplewex = 0.2,
        names = c(20,20,50,50,80,80),
        xlab = "Number of observations",
        ylab = "Posterior probability of H0",
        main = "Small effect with correlation = 0.8")

boxplot(sims58r,
        ylim=c(0,1),
        outline=FALSE,
        col = rep(c("white","dark gray"),3),
        staplewex = 0.2,
        names = c(20,20,50,50,80,80),
        xlab = "Number of observations",
        ylab = "Posterior probability of H0",
        main = "Medium effect with correlation = 0.8")

# export as PDF, 6 (width) x 8.5 (height) in



# construct table of correct model choice proportions

# convert posterior model probabilities to log Bayes factors (null/alternative)
bf02 = log(sims02/(1-sims02))
bf08 = log(sims08/(1-sims08))
bf22 = log(sims22/(1-sims22))
bf28 = log(sims28/(1-sims28))
bf52 = log(sims52/(1-sims52))
bf58 = log(sims58/(1-sims58))

# compute proportion choose "correct" model
# log BF > 0 if null effect
# log BF < 0 if non-null effect
# key: original20, original50, original80, improved20, improved50, improved80

props02 = apply(bf02 > 0, 2, sum)/1000; props02
props08 = apply(bf08 > 0, 2, sum)/1000; props08
props22 = apply(bf22 < 0, 2, sum)/1000; props22
props28 = apply(bf28 < 0, 2, sum)/1000; props28
props52 = apply(bf52 < 0, 2, sum)/1000; props52
props58 = apply(bf58 < 0, 2, sum)/1000; props58

# compute consistency between methods

# effect=0, rho = 2
sum(bf02[,1]*bf02[,4]>0)/1000 # N=20
sum(bf02[,2]*bf02[,5]>0)/1000 # N=50
sum(bf02[,3]*bf02[,6]>0)/1000 # N=80

# effect=0, rho = 8
sum(bf08[,1]*bf08[,4]>0)/1000 # N=20
sum(bf08[,2]*bf08[,5]>0)/1000 # N=50
sum(bf08[,3]*bf08[,6]>0)/1000 # N=80

# effect=2, rho=2
sum(bf22[,1]*bf22[,4]>0)/1000 # N=20
sum(bf22[,2]*bf22[,5]>0)/1000 # N=50
sum(bf22[,3]*bf22[,6]>0)/1000 # N=80

# effect=2, rho=8
sum(bf28[,1]*bf28[,4]>0)/1000 # N=20
sum(bf28[,2]*bf28[,5]>0)/1000 # N=50
sum(bf28[,3]*bf28[,6]>0)/1000 # N=80

# effect=5, rho=2
sum(bf52[,1]*bf52[,4]>0)/1000 # N=20
sum(bf52[,2]*bf52[,5]>0)/1000 # N=50
sum(bf52[,3]*bf52[,6]>0)/1000 # N=80

# effect=5, rho=8
sum(bf58[,1]*bf58[,4]>0)/1000 # N=20
sum(bf58[,2]*bf58[,5]>0)/1000 # N=50
sum(bf58[,3]*bf58[,6]>0)/1000 # N=80


# correlations of PMP between methods
# and scatterplots (use N=50, other plots have similar shape)


plot(sims02[,2], sims02[,5], pch=3, lwd=1.5,
     bty="n",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab = "Minimal BIC method",
     ylab = "Nathoo & Masson method",
     main = "Null effect with correlation = 0.2"
); abline(a=0, b=1, lty=2)

plot(sims22[,2], sims22[,5], pch=3, lwd=1.5,
     bty="n",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab = "Minimal BIC method",
     ylab = "Nathoo & Masson method",
     main = "Small effect with correlation = 0.2"
); abline(a=0, b=1, lty=2)

plot(sims52[,2], sims52[,5], pch=3, lwd=1.5,
     bty="n",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab = "Minimal BIC method",
     ylab = "Nathoo & Masson method",
     main = "Medium effect with correlation = 0.2"
); abline(a=0, b=1, lty=2)

plot(sims08[,2], sims08[,5], pch=3, lwd=1.5,
     bty="n",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab = "Minimal BIC method",
     ylab = "Nathoo & Masson method",
     main = "Null effect with correlation = 0.8"
); abline(a=0, b=1, lty=2)

plot(sims28[,2], sims28[,5], pch=3, lwd=1.5,
     bty="n",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab = "Minimal BIC method",
     ylab = "Nathoo & Masson method",
     main = "Small effect with correlation = 0.8"
); abline(a=0, b=1, lty=2)

plot(sims58[,2], sims58[,5], pch=3, lwd=1.5,
     bty="n",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab = "Minimal BIC method",
     ylab = "Nathoo & Masson method",
     main = "Medium effect with correlation = 0.8"
); abline(a=0, b=1, lty=2)


# compute correlation of PMP between methods

# effect=0, rho = 2
cor(sims02[,1], sims02[,4])
cor(sims02[,2], sims02[,5])
cor(sims02[,3], sims02[,6])

# effect=0, rho = 8
cor(sims08[,1], sims08[,4])
cor(sims08[,2], sims08[,5])
cor(sims08[,3], sims08[,6])

# effect=2, rho=2
cor(sims22[,1], sims22[,4])
cor(sims22[,2], sims22[,5])
cor(sims22[,3], sims22[,6])

# effect=2, rho=8
cor(sims28[,1], sims28[,4])s
cor(sims28[,2], sims28[,5])
cor(sims28[,3], sims28[,6])

# effect=5, rho=2
cor(sims52[,1], sims52[,4])
cor(sims52[,2], sims52[,5])
cor(sims52[,3], sims52[,6])

# effect=5, rho=8
cor(sims58[,1], sims58[,4])
cor(sims58[,2], sims58[,5])
cor(sims58[,3], sims58[,6])








