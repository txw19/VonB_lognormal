
#################################
### Group-specific
#################################
# fake data to predict
predX <- seq(min(dat$age), max(dat$age), length=100) 


# Container for predicted values
est.lineB <- array(NA, c(out$BUGSoutput$n.sims,length(predX),J) )

# Put each groups MCMC draws for all parameters in its own list
group.params <- list()
for(m in 1:J){
  group.params[[m]] <- out$BUGSoutput$sims.list$BB[,m,]
}

# Re-transform: look at - str(group.params)
for(j in 1:J){
  group.params[[j]][,1] <- exp(group.params[[j]][,1] )
  group.params[[j]][,2] <- exp(group.params[[j]][,2] )
  group.params[[j]][,3] <- exp(group.params[[j]][,3] )-10
}



for(k in 1:J){ # loop over groups (J)
  for(i in 1:out$BUGSoutput$n.sims ){  
    for(t in 1:length(predX)){
      # est.lineB[i,t,k] <-  group.params[[k]][i,1] + group.params[[k]][i,2] * predX[t]
      est.lineB[i,t,k] <-  group.params[[k]][i,1] * (1 -exp(-group.params[[k]][i,2] * (predX[t] - group.params[[k]][i,3])))
    }	  
  }
}

groupMean <- array(NA, c(1,length(predX),J) )
upper.CIB <- array(NA, c(1,length(predX),J) )
lower.CIB <- array(NA, c(1,length(predX),J) )

for(i in 1:J){
  # Means
  groupMean[,,i] <- apply(est.lineB[,,i], 2, mean )
  # 95% CIs for fitted values
  upper.CIB[,,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.975) )
  lower.CIB[,,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.025) )
}



################# PLOT ############################

res <- 6
name_figure <- "Group_Specific_vonB_male.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 1
x.label = 'Age (yrs)'
y.label = "Length (mm)"

nf <- layout(matrix(c(1:10),nrow=2,ncol=5,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Add transparency to points
source('AddTransFunction.R')


# Reservoir names
res_names <- unique(dat$reservoir)

# Group-specific plots

# Make variable to loop through reservoirs
dat$res2 <- as.numeric(dat$reservoir)

for(i in 1:J){
  
  plot(predX,z, ylim=c(min(dat$length,na.rm=T),max(dat$length,na.rm=T)),
       xlim=c(min(dat$age),max(dat$age)), axes=F, ylab='', xlab='', type='n')
  
  points(jitter(dat$age[dat$res2==i]), dat$length[dat$res2==i], cex=0.8, pch=16,col='black' )
  
  
  if( i <=5){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==1 | i==6 ){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1)
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  # Add credible region
  i.for <- order(predX)
  i.back <- order(predX, decreasing = TRUE )
  x.polygon <- c( predX[i.for] , predX[i.back] )
  y.polygon <- c( lower.CIB[,,i][i.for] , upper.CIB[,,i][i.back] )
  polygon( x.polygon , y.polygon , col = addTrans('black',100) , border = NA )
  
  # Add posterior means
  lines(predX, groupMean[,,i],lwd=1, col='black',lty=1)
  
  text(5,500,res_names[i],cex=1.2)
  box()
  
}

mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)


par(def.par)
dev.off()
