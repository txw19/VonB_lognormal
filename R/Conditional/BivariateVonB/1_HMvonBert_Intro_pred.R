#rm(list=ls())


#### Load libraries
library(mgcv)
library(MCMCpack) # rwish function
library(R2jags)


#read in data
# setwd("S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert")
dat <- read.csv("FHCGrowthMetaData.csv")

#remove old data from analysis
dat <- subset(dat, Year >= 2000) 
sort(unique(dat$Year)) #check


#only interested in intro pops
intro <-subset(dat, Status=="introduced")
#removing lakes from intro pops did not change relationship, so keeping them in model
# intro2<-intro[intro$River!="Lake Marion", ]
# intro2<-intro2[intro2$River!="Lake Moultrie", ]


#must drop levels before continuing
intro <- droplevels(intro)
# intro2 <- droplevels(intro2)


#river name column
intro <- intro[order(intro$River), ]
intro$River_name <- intro$River
intro$River <- as.numeric(intro$River)

#final check:
length(unique(intro$River)) #should be 12 for intro pops
unique(intro$River) #should be in a consecutive numeric order
unique(intro$River_name) #does not include removed pops above or NATIVES
sort(unique(intro$TimeSinceEstab))
mean(intro$TimeSinceEstab)
###################################ADD PREDICTORS
# Site-level predictors
estab <- as.numeric(by(intro$TimeSinceEstab, intro$River, mean)) 
z_estab <- as.numeric(scale(estab))

#################################################################
########## BUGS CODE ############################################
#################################################################
sink("MetaAnalysisHMvonBmodel_intro_pred.txt")
cat("
    model{
    for(i in 1:n){
    y[i] ~ dnorm(y.hat[i], tau.y)
    y.hat[i] <- Linf[g[i]] * (1-exp(-k[g[i]] * (age[i] - t0[g[i]] )))
    }
    
    tau.y <- pow(sigma.y,-2)
    sigma.y ~ dunif(0,100)
    
    # Linf and K modeled on log-scale
    for(j in 1:J){
    Linf[j] <- exp(BB[j,1])
    k[j] <- exp(BB[j,2])
    t0[j] ~ dnorm(mu.t0, tau.t0)
    BB[j,1:K] ~ dmnorm (BB.hat[j,], Tau.B[,]) # Multivariate normal dist'n;  Tau.B is a precision 
    BB.hat[j,1] <- mu.Linf + Lgamma1 * time[j] 
    BB.hat[j,2] <- mu.k + Kgamma1 * time[j] 
    }
    
    # Priors for population-average parameters
    tau.t0 <- pow(sigmat0,-2)
    sigmat0 ~ dunif(0.001,1)

    mu.Linf ~ dnorm(0,.0001)
    mu.k ~ dnorm(0,.0001)
    mu.t0 ~ dnorm(0,.0001)
    Lgamma1 ~ dnorm(0,.0001)
    Kgamma1 ~ dnorm(0,.0001)
    
    
    
    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    
    } # end model
    ",fill=TRUE)
sink()

################
J <- length(unique(intro$River))


# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 2

# Create identity matrix for Wishart dist'n
W <- diag(K)

# load data
data <- list(y = intro$TL, age = intro$Age, g = intro$River, n = dim(intro)[1],
             J = J, W=W, K=K, time=z_estab)

# Initial values
inits <- function(){list(mu.Linf = rnorm(1,3,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.7,0.001),Lgamma1=rnorm(1),Kgamma1=rnorm(1), 
                         sigma.y = runif(1,1,10),  
                         BB=array(c(rep(log(950) +rnorm(1,0.01,0.01),J),rep(log(0.04)+rnorm(1,0.001,0.1),J)),
                                  c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }

# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B",
             "Lgamma1","Kgamma1","sigmat0" )
rep(log(-2+10)+rnorm(1,0.01,0.1),J)

# MCMC settings
ni <- 100000
nt <- 3
nb <- 50000
nc <- 3


############################################
start.time = Sys.time()  

intro.pred.out <- jags(data = data, inits = inits, parameters.to.save = params1, 
            model.file = "MetaAnalysisHMvonBmodel_intro_pred.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)

#Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 


write.csv(intro.pred.out$BUGSoutput$summary,"BUGSoutBivariate.csv")
