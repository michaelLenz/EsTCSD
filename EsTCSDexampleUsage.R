source("EsTCSDcalculations.R")
######################
# Plot real and observed densities for selected uniform, normal, and gamma distributions (Figure 2)
######################
# change code of measurementDensity to uniform distribution
distributionFunction = function(z, a, b){
	return(dunif(z,a,b))
}
meanDistributionFunction = function(a,b){
	return(0.5*(a+b))
}
sdDistributionFunction = function(a,b){
	return((b-a)/sqrt(12))
}
x = 0.1*(1:1050); zU = rep(NA,105); for (i in 1:length(x)){zU[i] = measurementDensity(x[i],8,25,10,100)}
zU2 = rep(NA,length(x)); for (i in 1:length(x)){zU2[i] = measurementDensity(x[i],8,70,10,100)}
zU3 = rep(NA,length(x)); for (i in 1:length(x)){zU3[i] = measurementDensity(x[i],40,25,10,100)}
library(grDevices)
tiff("Figure2.tiff",width=9,height=3.5,units="in",res=300,compression="lzw")
par(mfrow=c(1,3))
plot(2*x,dunif(x,10,100)/2,type="l",lwd=2,lty=2,xlab="diameter",ylab="density",ylim=c(0,max(zU)/3*2.5),main="uniform distribution")
lines(2*x,zU/2,lwd=2,lty=1)
lines(2*x,zU2/2,lwd=2,lty=1,col=2)
lines(2*x,zU3/2,lwd=2,lty=1,col=4)
legend("topleft",legend=c("real","measured"),lwd=2,lty=c(2,1))
legend("topright",legend=c(expression(list(paste("c"[0],"= 25"),paste("h"[s],"= 8"))),expression(list(paste("c"[0],"= 70"),paste("h"[s],"= 8"))),expression(list(paste("c"[0],"= 25"),paste("h"[s],"= 40")))),lwd=2,lty=c(1),col=c(1,2,4))
axis(side=1,at=70,labels="",col=2)
mtext(expression("c"[0]), at = 70, side = 1, line = 1, col = 2,cex=0.7)
axis(side=1,at=25,labels=expression("c"[0]))

# change code of measurementDensity to normal distribution
distributionFunction = function(z, mu, sigma){
	return(dnorm(z,mu,sigma))
}
meanDistributionFunction = function(mu, sigma){
	return(mu)
}
sdDistributionFunction = function(mu, sigma){
	return(sigma)
}
x = 0.1*(1:1500); zU = rep(NA,length(x)); for (i in 1:length(x)){zU[i] = measurementDensity(x[i],8,25,75,20)}
zU2 = rep(NA,length(x)); for (i in 1:length(x)){zU2[i] = measurementDensity(x[i],8,25,37.5,20)}
zU3 = rep(NA,length(x)); for (i in 1:length(x)){zU3[i] = measurementDensity(x[i],8,25,75,40)}
plot(2*x,dnorm(x,75,20)/2,type="l",lwd=2,lty=2,xlab="diameter",ylab="density",ylim=c(0,max(zU)/3*2.5),main="normal distribution")
lines(2*x,zU/2,lwd=2,lty=1)
lines(2*x,dnorm(x,37.5,20)/2,lwd=2,lty=2,col=2)
lines(2*x,zU2/2,lwd=2,lty=1,col=2)
lines(2*x,dnorm(x,75,40)/2,lwd=2,lty=2,col=4)
lines(2*x,zU3/2,lwd=2,lty=1,col=4)
legend("topleft",legend=c("real","measured"),lwd=2,lty=c(2,1))
legend("topright",legend=c(expression(list(paste(mu,"= 150"),paste(sigma,"= 40"))),expression(list(paste(mu,"= 75"),paste(sigma,"= 40"))),expression(list(paste(mu,"= 150"),paste(sigma,"= 80")))),lwd=2,lty=2,col=c(1,2,4))
axis(side=1,at=25,labels=expression("c"[0]))

# change code of measurementDensity to gamma distribution
distributionFunction = function(z, k, theta){
	return(dgamma(z,k,scale=theta))
}
meanDistributionFunction = function(k, theta){
	return(k*theta)
}
sdDistributionFunction = function(k, theta){
	return(theta*sqrt(k))
}
x = 0.1*(1:1500); zU = rep(NA,length(x)); for (i in 1:length(x)){zU[i] = measurementDensity(x[i],8,25,20,3)}
zU2 = rep(NA,length(x)); for (i in 1:length(x)){zU2[i] = measurementDensity(x[i],8,25,5,12)}
zU3 = rep(NA,length(x)); for (i in 1:length(x)){zU3[i] = measurementDensity(x[i],8,25,2,30)}
plot(2*x,dgamma(x,20,scale=3)/2,type="l",lwd=2,lty=2,xlab="diameter",ylab="density",ylim=c(0,max(zU)/3*2), main="gamma distribution")
lines(2*x,zU/2,lwd=2,lty=1)
lines(2*x,dgamma(x,5,scale=12)/2,lwd=2,lty=2,col=2)
lines(2*x,zU2/2,lwd=2,lty=1,col=2)
lines(2*x,dgamma(x,2,scale=30)/2,lwd=2,lty=2,col=4)
lines(2*x,zU3/2,lwd=2,lty=1,col=4)
legend("topright",legend=c("real","measured"),lwd=2,lty=c(2,1))
legend("right",legend=c(expression(list("k = 20",paste(theta,"= 6"))),expression(list("k = 5",paste(theta,"= 24"))),expression(list("k = 2",paste(theta,"= 60")))),lwd=2,lty=2,col=c(1,2,4))
axis(side=1,at=25,labels=expression("c"[0]))
dev.off()

################################
# Perform correction (maximum likelihood) on measured cell size data (data from one individual are provided in "cellDiameters.Rdata")
################################
# load example data:
load("cellDiameters.Rdata")
ML = maximumLikelihoodGamma(cellDiameters/2, pars_init=c(mean(cellDiameters/2)^2/var(cellDiameters/2)+2,var(cellDiameters/2)/mean(cellDiameters/2)), hs=8, c0=25)
library(grDevices)
tiff("MLexample.tiff",width=5,height=5,res=300,units="in",compression="lzw")
h1 = hist(cellDiameters,breaks=30,xlab=expression(paste("Diameter [", mu, "m]")),freq=F,main="Example correction",ylim=c(0,0.025),xlim=c(20,125),ylab="Density")
axis(side=4,at=c(0,10,20,30,40,50)/h1$counts[1]*h1$density[1],labels=c("0","10","20","30","40","50"))
mtext("Frequency",side=4,line=3)
x = seq(min(h1$mids),max(h1$mids),len=100); lines(x,dgamma(x/2,ML$solution[1],scale=ML$solution[2])/2,lwd=2)
z = rep(NA,length(x));
for (i in 1:length(x)){
	z[i] = measurementDensity(x[i]/2,8,25,ML$solution[1],ML$solution[2])
}
lines(x,z/2,lty=2,lwd=2)
legend("topright",legend=c("fitted","corrected"),lty=c(2,1),bty="n",lwd=2)
dev.off()
