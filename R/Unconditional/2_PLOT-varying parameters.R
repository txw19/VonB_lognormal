###############################################################
##### Plot Popl'n average effect
# fake data to predict
predX <- seq(min(dat$age), max(dat$age), length=100) 

# Extract pop'n average coefficents

PopAve <- out$BUGSoutput$summary[c("mu.Linf", "mu.k", "mu.t0"),1]

# Re-transform
PopAve[1] <- exp(PopAve[1])
PopAve[2] <- exp(PopAve[2])
PopAve[3] <- exp(PopAve[3])-10

GroupCoef <- matrix(out$BUGSoutput$summary[1:(3*(length(unique(dat$reservoir)) )),1], c(length(unique(dat$reservoir)),3), byrow=F)

# Re-transform
GroupCoef[,1] <- exp(GroupCoef[,1])
GroupCoef[,2] <- exp(GroupCoef[,2])
GroupCoef[,3] <- exp(GroupCoef[,3])-10


# Generate fake y-axis for creating plot
z <- seq(min(dat$age),max(dat$age),length=100) 

# Create Popl'n average mean fitted line
y.pred <- PopAve[1] * (1-exp(-PopAve[2]  * (predX -  PopAve[3] )))


#####---------------- Create Figure
res <- 6
name_figure <- "VonB_fit_male.png"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)

par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

size.labels = 1
size.text = 1
x.label = 'Age (years)'
y.label = 'Length (mm)'


plot(predX,z, ylim=c(min(dat$length),max(dat$length)),xlim=c(min(dat$age),max(dat$age)), axes=F, ylab='', xlab='', type='n')
points(jitter(dat$age),dat$length, cex=0.8, pch=16)

# Add group-specific fits
for(i in 1:J){
  y.pred2 <- GroupCoef[i,1] * (1-exp(-GroupCoef[i,2]  * (predX -  GroupCoef[i,3] )))
	lines(predX, y.pred2,lwd=0.5, col='black',lty=1)
}

# Add popl'n average fit
lines(predX, y.pred, lwd = 5, col="blue", lty = 1)

axis(side=1, cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0) ) 
axis(side=2,cex.axis=size.text,font=1 ,tck=-0.01, mgp=c(0,0.5,0), las=1) 
mtext(x.label, line = 1.3, side = 1, cex = size.text,outer=T)
mtext(y.label, line = 1.8, side = 2, cex = size.text,outer=T)
box()


par(def.par)
dev.off()


