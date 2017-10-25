######################################################
########## SCATTER PLOT FOR LEVEL 2 MODEL FIT#########
# Select random slopes  
mean.beta <- exp(intro.pred.out$BUGSoutput$mean$BB[,1] )

# Fake data to predict
fake1 <- seq(min(z_estab), max(z_estab), length=20)


# Obtain 90% CIs for fitted line
est.lineA <- matrix(NA, ncol=length(fake1), nrow=intro.pred.out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:intro.pred.out$BUGSoutput$n.sims ){
  for(t in 1:length(fake1) ){
    est.lineA[i,t] <- exp(intro.pred.out$BUGSoutput$sims.list$mu.Linf[i] + intro.pred.out$BUGSoutput$sims.list$Lgamma1[i] * fake1[t] )
  }
}

# CIs for fitted values
upper.CIA <- apply(est.lineA, 2, quantile, probs=c(0.95) )
lower.CIA <- apply(est.lineA, 2, quantile, probs=c(0.05) )
post.mean <- apply(est.lineA, 2, mean )

## Grab 90% CIs for beta's
u.alpha <- numeric(length(mean.beta) )
l.alpha <- numeric(length(mean.beta) )
for(i in 1:length(mean.beta)) { 
  u.alpha[i] <- exp(quantile(intro.pred.out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.95) ) )
  l.alpha[i] <- exp(quantile(intro.pred.out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.05) ) )
}

###########################################
####### FIGURE WITH CRI's
res <- 6
name_figure <- "linf_established.png"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...
par(mfrow = c(1,1), mar=c(3,3,1,1)) 


size.labels = 1
size.text = 1

x.label <- 'Established'
# y.label <- expression(paste('Slope (',beta[j],')'))
y.label <- expression(paste('Lake-specific growth rate (', K[j],')' ))

#xlab1 <- c(-3,-2,-1,0,1)
#xlab2 <- xlab1* sd(lat1) + mean(lat1)

plot(mean.beta ~ z_estab,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c( lower.CIA[i.for] , upper.CIA[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(z_estab, mean.beta,pch=16,cex=0.8)

segments(x0=z_estab, x1=z_estab,
         y0=l.alpha, y1=u.alpha, col='black',lwd=1)


lines(fake1,post.mean, lwd = 3, col="black", lty = 1)

mtext(x.label, line = 1.5, side = 1, cex = size.text)
mtext(y.label, line = 1.7, side = 2, cex = size.text)

box()
par(def.par)
dev.off()

####### END