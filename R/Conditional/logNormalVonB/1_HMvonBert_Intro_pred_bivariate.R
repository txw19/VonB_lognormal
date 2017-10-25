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
    y[i] ~ dlnorm(y.hat[i], tau.y)
    y.hat[i] <- log(Linf[g[i]] * (1-exp(-k[g[i]] * (age[i] - t0[g[i]] ))))
    }
    
    tau.y <- pow(sigma.y,-2)
    sigma.y ~ dunif(0,10)
    
    # Parameters modeled on log-scale
    for(j in 1:J){
    Linf[j] <- exp(BB[j,1])
    k[j] <- exp(BB[j,2])
    t0[j] ~ dnorm(mu.t0,tau.t0)T(-2,0)
    BB[j,1:K] ~ dmnorm (BB.hat[j,], Tau.B[,]) 
    BB.hat[j,1] <- mu.Linf + Lgamma1 * time[j]
    BB.hat[j,2] <- mu.k + Kgamma1 * time[j]
    
    }
    
    Lgamma1 ~ dnorm(0,.0001)
    Kgamma1 ~ dnorm(0,.0001)
    
    tau.t0 <- pow(sigmat0,-2)
    sigmat0 ~ dunif(0.001,1)
    # Priors for population-average parameters
    
    mu.Linf ~ dnorm(0,0.001)
    mu.k ~ dnorm(0.0,0.001)
    # mu.Linf ~ dunif(1300,1600)
    # mu.k ~ dunif(0.01,0.5)
    mu.t0 ~ dunif(-2,0)
    
    
    
    
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
inits <- function(){list(Tau.B=rwish(K+1,diag(K)) ) }

# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B",
             "Lgamma1","Kgamma1","Linf","k","t0" )


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
##############################

print(intro.pred.out)

write.csv(intro.pred.out$BUGSoutput$summary,"BUGSbivOut.csc")

