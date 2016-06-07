############################################################################################################################
# Example of how to use the EstiTiCS (Estimate Tissue Cell Size/type distribution) R package
#
# This example consists of two parts:
# Part 1: Simulation study of real and measured (cross-sectional) cell size distributions (reproduces Figure 3 from the EstiTiCS paper)
# Part 2: Correction of adipose tissue cell size distribution using a maximum likelihood approach (using data from one individual and one clinical investigation day of the EstiTiCS paper)
# 
# For more details on EstiTiCS see:
# Lenz M, Roumans NJT, Vink RG, van Baak M, Mariman E, Arts ICW, de Kok TM, Ertaylan G
# "Estimating real cell size distribution from cross section microscopy imaging"
# 
# Contact: Michael Lenz (michael.lenz@maastrichtuniversity.nl) or 
#          Gökhan Ertaylan (gokhan.ertaylan@maastrichtuniversity.nl)
#		   Maastricht Centre for Systems Biology, Maastricht University, The Netherlands (macsbio.maastrichtuniversity.nl)
############################################################################################################################

source("EstiTiCScalculations.R") # load EstiTiCS functions

############################################
# Part 1:
# Plot real and observed densities for selected uniform, normal, and gamma distributions (Figure 3 of the EstiTiCS paper)
############################################
# change distribution function to uniform distribution
distributionFunction = function(z, a, b){
	return(dunif(z,a,b))
}
meanDistributionFunction = function(a,b){
	return(0.5*(a+b))
}
sdDistributionFunction = function(a,b){
	return((b-a)/sqrt(12))
}

# Calculate measured cell size distribution based on underlying true uniform distribution
x = 0.1*(1:1050); # range of cell radii for which the density should be calculated
zU = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU[i] = measurementDensity(x[i],8,25,10,100) # calculate density for slice thickness of 8, imaging cutoff of 25, and a uniform true cell radii distribution with support on the interval [10,100] - cell diameters have a support on [20,200]
}
zU2 = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU2[i] = measurementDensity(x[i],8,70,10,100) # same as before, but with different imaging cutoff (70 instead of 25)
}
zU3 = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU3[i] = measurementDensity(x[i],40,25,10,100) # same as for zU, but with different slice thickness (40 instead of 8)
}

# Plot probability densities (true and measured)
library(grDevices)
tiff("Figure3.tiff",width=9,height=3.5,units="in",res=300,compression="lzw")
par(mfrow=c(1,3))
# plot true distribution (uniform distribution)
plot(2*x,dunif(x,10,100)/2,type="l",lwd=2,lty=2,xlab="diameter",ylab="density",ylim=c(0,max(zU)/3*2.5),main="uniform distribution")
# plot measured distributions for the three different parameter settings
lines(2*x,zU/2,lwd=2,lty=1)
lines(2*x,zU2/2,lwd=2,lty=1,col=2)
lines(2*x,zU3/2,lwd=2,lty=1,col=4)
legend("topleft",legend=c("real","measured"),lwd=2,lty=c(2,1))
legend("topright",legend=c(expression(list(paste("c"[0],"= 25"),paste("h"[s],"= 8"))),expression(list(paste("c"[0],"= 70"),paste("h"[s],"= 8"))),expression(list(paste("c"[0],"= 25"),paste("h"[s],"= 40")))),lwd=2,lty=c(1),col=c(1,2,4))
axis(side=1,at=70,labels="",col=2)
mtext(expression("c"[0]), at = 70, side = 1, line = 1, col = 2,cex=0.7)
axis(side=1,at=25,labels=expression("c"[0]))

# change distribution function to normal distribution
distributionFunction = function(z, mu, sigma){
	return(dnorm(z,mu,sigma))
}
meanDistributionFunction = function(mu, sigma){
	return(mu)
}
sdDistributionFunction = function(mu, sigma){
	return(sigma)
}

# Calculate measured cell size distribution based on underlying true normal distribution
x = 0.1*(1:1500); # range of cell radii for which the density should be calculated
zU = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU[i] = measurementDensity(x[i],8,25,75,20); # calculate density for slice thickness of 8, imaging cutoff of 25, and a normal true cell radii distribution with mean 75 and standard deviation 20 - cell diameters have a mean of 150 and standard deviation of 40
}
zU2 = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU2[i] = measurementDensity(x[i],8,25,37.5,20); # same as before, but with mean of cell radii mean of 37.5 (75 for cell diameters)
}
zU3 = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU3[i] = measurementDensity(x[i],8,25,75,40); # same as for zU, but with standard deviation of 40 (80 for cell diameters)
}

# plot true distribution (normal distribution with mean = 75 (150), sd = 20 (40))
plot(2*x,dnorm(x,75,20)/2,type="l",lwd=2,lty=2,xlab="diameter",ylab="density",ylim=c(0,max(zU)/3*2.5),main="normal distribution")
# plot corresponding measured distribution
lines(2*x,zU/2,lwd=2,lty=1)
# plot true distribution (normal distribution with mean = 37.5 (75), sd = 20 (40))
lines(2*x,dnorm(x,37.5,20)/2,lwd=2,lty=2,col=2)
# plot corresponding measured distribution
lines(2*x,zU2/2,lwd=2,lty=1,col=2)
# plot true distribution (normal distribution with mean = 75 (150), sd = 40 (80))
lines(2*x,dnorm(x,75,40)/2,lwd=2,lty=2,col=4)
# plot corresponding measured distribution
lines(2*x,zU3/2,lwd=2,lty=1,col=4)
legend("topleft",legend=c("real","measured"),lwd=2,lty=c(2,1))
legend("topright",legend=c(expression(list(paste(mu,"= 150"),paste(sigma,"= 40"))),expression(list(paste(mu,"= 75"),paste(sigma,"= 40"))),expression(list(paste(mu,"= 150"),paste(sigma,"= 80")))),lwd=2,lty=2,col=c(1,2,4))
axis(side=1,at=25,labels=expression("c"[0]))

# change distribution function to gamma distribution
distributionFunction = function(z, k, theta){
	return(dgamma(z,k,scale=theta))
}
meanDistributionFunction = function(k, theta){
	return(k*theta)
}
sdDistributionFunction = function(k, theta){
	return(theta*sqrt(k))
}

# Calculate measured cell size distribution based on underlying true gamma distribution
x = 0.1*(1:1500); # range of cell radii for which the density should be calculated
zU = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU[i] = measurementDensity(x[i],8,25,20,3); # calculate density for slice thickness of 8, imaging cutoff of 25, and a gamma-shaped true cell radii distribution with shape parameter 20 and scale parameter 3 - cell diameters have a shape parameter of 20 and scale parameter of 6
}
zU2 = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU2[i] = measurementDensity(x[i],8,25,5,12); # same as before, but with shape parameter 5 and scale parameter 12 (24 for cell diameter)
}
zU3 = rep(NA,length(x)); # variable to store probability density values
for (i in 1:length(x)){
	zU3[i] = measurementDensity(x[i],8,25,2,30); # same as before, but with shape parameter of 2 and scale parameter of 30 (60 for cell diameter)
}

# plot true distribution (gamma distribution with shape k = 20 and scale theta = 3 (6))
plot(2*x,dgamma(x,20,scale=3)/2,type="l",lwd=2,lty=2,xlab="diameter",ylab="density",ylim=c(0,max(zU)/3*2), main="gamma distribution")
# plot corresponding measured distribution
lines(2*x,zU/2,lwd=2,lty=1)
# plot true distribution (gamma distribution with shape k = 5 and scale theta = 12 (24))
lines(2*x,dgamma(x,5,scale=12)/2,lwd=2,lty=2,col=2)
# plot corresponding measured distribution
lines(2*x,zU2/2,lwd=2,lty=1,col=2)
# plot true distribution (gamma distribution with shape k = 2 and scale theta = 30 (60))
lines(2*x,dgamma(x,2,scale=30)/2,lwd=2,lty=2,col=4)
# plot corresponding measured distribution
lines(2*x,zU3/2,lwd=2,lty=1,col=4)
legend("topright",legend=c("real","measured"),lwd=2,lty=c(2,1))
legend("right",legend=c(expression(list("k = 20",paste(theta,"= 6"))),expression(list("k = 5",paste(theta,"= 24"))),expression(list("k = 2",paste(theta,"= 60")))),lwd=2,lty=2,col=c(1,2,4))
axis(side=1,at=25,labels=expression("c"[0]))
dev.off()

############################################
# Part 2:
# Perform correction (maximum likelihood) on measured cell size data (data from one individual are provided in "cellDiameters.Rdata")
############################################
source("EstiTiCScalculations.R") # (re)load EstiTiCS functions
# load example data:
load("cellDiameters.Rdata")

# use maximum likelihood approach to infer parameters of corrected cell radii distribution - use gamma distribution assumption
ML = maximumLikelihoodGamma(cellDiameters/2, pars_init=c(mean(cellDiameters/2)^2/var(cellDiameters/2)+2,var(cellDiameters/2)/mean(cellDiameters/2)), hs=8, c0=25)

# visualization of measured and corrected distributions
library(grDevices)
tiff("MLexample.tiff",width=5,height=5,res=300,units="in",compression="lzw")
# plot histogram of measured cell diameters
h1 = hist(cellDiameters,breaks=30,xlab=expression(paste("Diameter [", mu, "m]")),freq=F,main="Example correction",ylim=c(0,0.025),xlim=c(20,125),ylab="Density")
axis(side=4,at=c(0,10,20,30,40,50)/h1$counts[1]*h1$density[1],labels=c("0","10","20","30","40","50"))
mtext("Frequency",side=4,line=3)
# plot corrected ("estimated true") cell diameter distribution
x = seq(min(h1$mids),max(h1$mids),len=100); lines(x,dgamma(x/2,ML$solution[1],scale=ML$solution[2])/2,lwd=2)
# determine and plot corresponding measurement distribution
z = rep(NA,length(x));
for (i in 1:length(x)){
	z[i] = measurementDensity(x[i]/2,8,25,ML$solution[1],ML$solution[2])
}
lines(x,z/2,lty=2,lwd=2)
legend("topright",legend=c("fitted","corrected"),lty=c(2,1),bty="n",lwd=2)
dev.off()

# calculate total number of adipocytes using total fat mass of the person (here: 31.37kg)
fm = 31.37 # kg
numCells = numCellsCorrected(fm,ML$solution[1],ML$solution[2])
