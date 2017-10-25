# rm(list=ls())

#### Load libraries
library(MCMCpack) # rwish function
library(R2jags)

# Read in data
dat<-read.csv('GrowthDataset.csv',na.strings='NA', header=T)
head(dat)
dim(dat)


# Quick and ugly plots
# By sex across all reservioirs
plot(length ~ jitter(age), data=dat, col=dat$sex,pch=16,xlab='Age (yrs)',ylab='Length (mm)')
legend(16,300,c("Female","Male"),pch=c(16,16),
              col=c(1,2))

# By reservoir
plot(length ~ jitter(age), data=dat, col=dat$reservoir,pch=16,xlab='Age (yrs)',ylab='Length (mm)')


# Grab female data
dat <- dat[dat$sex=='m',]
summary(dat)
dat <- droplevels(dat)
summary(dat)


#################################################################
########## BUGS CODE ############################################
#################################################################


sink("vonBmodel.txt")
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
	BB.hat[j,1] <- log(mu.Linf)
	BB.hat[j,2] <- log(mu.k)

}

 tau.t0 <- pow(sigmat0,-2)
 sigmat0 ~ dunif(0.001,1)
  # Priors for population-average parameters


  mu.Linf ~ dunif(200,800)
  mu.k ~ dunif(0.01,0.5)
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

# Number of sites
J <-length(unique(dat$reservoir))

# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 2

# Create identity matrix for Wishart dist'n
W <- diag(K)

# load data
data <- list(y = dat$length, age = dat$age, g = as.numeric(dat$reservoir), n = dim(dat)[1],
             J = J, W=W, K=K )


# Initial values
# inits <- function(){list(mu.Linf = rnorm(1,500,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.5,0.001),
#                          sigma.y = runif(1,1,10), 
#                     BB=array(c(rep(500,J),rep(0.4,J),rep(0.5,J)),
#                                  c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }

inits <- function(){list(Tau.B=rwish(K+1,diag(K)) ) }

# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B","Linf","k","t0")


# MCMC settings
ni <- 30000
nt <- 3
nb <- 20000
nc <- 3


############################################
###### DO analysis in JAGS 
start.time = Sys.time()         # Set timer (1.12 mins)

out <- jags(data = data, inits = inits, parameters.to.save = params1, 
            model.file = "vonBmodel.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb)

end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize the result
print(out, digits = 2)
# str(out)

write.csv(out$BUGSoutput$summary, "BUGSout.csv",row.names = T)

# Find which parameters, if any, have Rhat > 1.1
which(out$BUGSoutput$summary[, c("Rhat")] > 1.1)


# BB[samples(1:nsim),sites(1:10), parameters(1:3)]
# mean(out$BUGSoutput$sims.list$BB[,1,1] - out$BUGSoutput$sims.list$BB[,2,1])
# quantile(out$BUGSoutput$sims.list$BB[,1,1] - out$BUGSoutput$sims.list$BB[,3,1],c(0.025,0.975))


# Multiple comparisons Linf among 10 reservoirs
Linf_Diff <- matrix(NA,nrow=J,ncol=J)
for (j in 1:J){
  for (jj in 1:J){
    
    diff <- out$BUGSoutput$sims.list$BB[,j,1] - out$BUGSoutput$sims.list$BB[,jj,1]
    diffCI <- quantile(diff,c(0.025,0.975))
    Linf_Diff[j,jj] <- (diffCI[1] * diffCI[2]) > 0 # Does the 95% CI overlap zero? TRUE = signifcant difference
  }
}
Linf_Diff


# quantile(out$BUGSoutput$sims.list$BB[,1,2] - out$BUGSoutput$sims.list$BB[,7,2],c(0.025,0.975))

# Multiple comparisons k among 10 reservoirs
k_Diff <- matrix(NA,nrow=J,ncol=J)
for (j in 1:J){
  for (jj in 1:J){
    
    diff <- out$BUGSoutput$sims.list$BB[,j,2] - out$BUGSoutput$sims.list$BB[,jj,2]
    diffCI <- quantile(diff,c(0.025,0.975))
    k_Diff[j,jj] <- (diffCI[1] * diffCI[2]) > 0 # Does the 95% CI overlap zero? TRUE = signifcant difference
  }
}
k_Diff


##################################################
######## RUN JAGS ON PARALLEL PROCESSORS #########
# Number of sites
J <-length(unique(dat$reservoir))

# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 3

# Create identity matrix for Wishart dist'n
W <- diag(3)
y <- dat$length
age = dat$age
g = as.numeric(dat$reservoir)
n = dim(dat)[1]
Tau.B1 <- rwish(K+1,diag(K)) # This is for inits, jags did not like the old way (NOTE naming to Tau.B1, not Tau.B)
# load data
data <- list('y', 'age', 'g', 'n',
             'J', 'W', 'K','Tau.B1')


# Initial values
inits <- function(){list(mu.Linf = rnorm(1,3,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.7,0.001),
                         sigma.y = runif(1,1,10), 
                         BB=array(c(rep(log(500) +rnorm(1,0.01,0.01),J),rep(log(0.4)+rnorm(1,0.001,0.1),J),rep(log(0.5+10)+rnorm(1,0.01,0.1),J)),
                                  c(J,K)), Tau.B=Tau.B1 ) }



# Parameters monitored
parameters <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B")



start.time = Sys.time()         # Set timer (0.41 mins)
# Call JAGS from R
out2 <- jags.parallel(data, inits = inits, parameters.to.save = parameters, model.file = "vonBmodel.txt", n.chains = 3,  
                      n.thin = 3, n.iter = 30000, n.burnin = 20000)
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='')

print(out2)

