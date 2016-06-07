############################################################################################################################
# Specification of functions used to estimate tissue cell size distributions from cross-section images via EstiTiCS 
# (Estimate Tissue Cell Size/type distribution)
# For an example usage of these functions see "EstiTiCSexampleUsage.R"
# 
# For more details on EstiTiCS see:
# Lenz M, Roumans NJT, Vink RG, van Baak M, Mariman E, Arts ICW, de Kok TM, Ertaylan G
# "Estimating real cell size distribution from cross section microscopy imaging"
# 
# Contact: Michael Lenz (michael.lenz@maastrichtuniversity.nl) or 
#          Gökhan Ertaylan (gokhan.ertaylan@maastrichtuniversity.nl)
#		   Maastricht Centre for Systems Biology, Maastricht University, The Netherlands (macsbio.maastrichtuniversity.nl)
############################################################################################################################

###########################################
# The assumed underlying true distribution function f_z(z) has to be specified by the user - defaults to gamma distribution
# Three examples (gamma, normal, and uniform distribution) are provided in the following:
###########################################

# Example 1 (default): Gamma distribution
distributionFunction = function(z, k, theta){
	return(dgamma(z,k,scale=theta))
}
meanDistributionFunction = function(k, theta){
	return(k*theta)
}
sdDistributionFunction = function(k, theta){
	return(theta*sqrt(k))
}

# Example 2: Normal distribution
#distributionFunction = function(z, mu, sigma){
#return(dnorm(z,mu,sigma))
#}
#meanDistributionFunction = function(mu, sigma){
#return(mu)
#}
#sdDistributionFunction = function(mu, sigma){
#return(sigma)
#}

# Example 3: Uniform distribution
#distributionFunction = function(z, a, b){
#return(dunif(z,a,b))
#}
#meanDistributionFunction = function(a,b){
#return(0.5*(a+b))
#}
#sdDistributionFunction = function(a,b){
#return((b-a)/sqrt(12))
#}

###########################################
# Implementation of the functions to determine the measured distribution from the given true cell radii distribution (Equation 6 of the EstiTiCS paper)
###########################################

# Calculation of the normalization factor (denomiator) of Equation 6 of the EstiTiCS paper:
normFunToIntegrate = function(z,c0,hs,...){
	# helper function to calculate the normalization factor
	return((sqrt(z^2-c0^2/4)-hs/2)*distributionFunction(z,...))
}
normalizationFactor = function(c0,hs,...){
	# Function to calculate the normalization factor - helper function to calculate measured cell radius distribution
	return(integrate(normFunToIntegrate,sqrt(c0^2/4+hs^2/4),meanDistributionFunction(...)+10*sdDistributionFunction(...),c0,hs,...)$value)
}

# Calculation of the measured cell radii distribution (density) given the density of true cell radii
funToIntegrate = function(z,x,hs,c0,...){
	# helper function to calculate measured cell radius distribution
	return(ifelse(x > sqrt(z^2-hs^2/4) | x < c0/2,0,distributionFunction(z,...)*x/sqrt(z^2-x^2)))
}

measurementDensity = function(x,hs,c0,...){
	# function calculating the measured cell radius density given the true cell radius distribution
	return(1/normalizationFactor(c0,hs,...) * integrate(funToIntegrate,sqrt(x^2+hs^2/4),meanDistributionFunction(...)+10*sdDistributionFunction(...),x,hs,c0,...)$value)
}

measurementDensityNoNorm = function(x,hs,c0,...){
	# helper function to calculate logLikelihood - calculates measured cell radius density without dividing by the normalization factor
	return(integrate(funToIntegrate,sqrt(x^2+hs^2/4),meanDistributionFunction(...)+10*sdDistributionFunction(...),x,hs,c0,...)$value)
}

###########################################
# Implementation of the maximum likelihood approach assuming a true underlying Gamma-shaped distribution
###########################################

logLikelihoodGamma = function(pars,measurement,hs,c0){
	# function calculating the log-likelihood for given data "measurements" and parameters "pars" based on a gamma distribution 
	k = pars[1]
	theta = pars[2]
	logLike = 0
	normFac = normalizationFactor(c0,hs,k,theta)
	for (i in 1:length(measurement)){
		temp = measurementDensityNoNorm(measurement[i],hs,c0,k,theta)/normFac;
		if (temp<=.Machine$double.xmin){ # smallest non-zero normalized floating-point number - avoid infinite values (log(0) == -Inf)
			logLike= logLike + log(.Machine$double.xmin)
		}else{
			logLike = logLike + log(temp)
		} 
	}
	return(-logLike)
}

library(nloptr) # allows for non-linear constraints
# specify maximum likelihood approach:
maximumLikelihoodGamma = function(measurement, pars_init, hs=8, c0=25){
	opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
			"xtol_rel" = 1.0e-8,
			"maxeval" = 1000)
	inequalityConstraints = function(x,measurement,c0,hs){sqrt(c0^2/4+hs^2/4)-(meanDistributionFunction(x[1],x[2])+2*sdDistributionFunction(x[1],x[2]))} # This constraint is used to avoid numerical problems during the optimization. It assumes that mean+2*sd of the true distribution are above the imaging cutoff. If this assumption would not be fulfilled, the determination of true cell sizes would not be reliable anyways.
	return(nloptr(x0=pars_init,eval_f=logLikelihoodGamma,lb=c(5,1),ub=c(25,8),eval_g_ineq=inequalityConstraints,opts=opts,measurement=measurement,hs=hs,c0=c0)) # note: lower bounds (lb) and upper bounds (ub) may have to be adjusted. Specification of a narrow range can improve the numerics of the optimization, but the bounds should not be too narrow in order to not miss the optimal solution.
}

###########################################
# Implementation of the maximum likelihood approach assuming a true underlying Normal distribution
###########################################

# note: change "distributionFunction", "meanDistributionFunction" and "sdDistributionFunction" to the case of a normal distribution (example 2, above)
logLikelihoodNormal = function(pars,measurement,hs,c0){
	# function calculating the log-likelihood for given data "measurements" and parameters "pars" based on a normal distribution
	# Important: "distributionFunction", "meanDistributionFunction" and "sdDistributionFunction" have to be adjusted accordingly!
	mu = pars[1]
	sigma = pars[2]
	logLike = 0
	normFac = normalizationFactor(c0,hs,mu,sigma)
	for (i in 1:length(measurement)){
		temp = measurementDensityNoNorm(measurement[i],hs,c0,mu,sigma)/normFac;
		if (temp<=.Machine$double.xmin){ # smallest non-zero normalized floating-point number - avoid infinite values (log(0) == -Inf)
			logLike= logLike + log(.Machine$double.xmin)
		}else{
			logLike = logLike + log(temp)
		} 
	}
	return(-logLike)
}

library(nloptr) # allows for non-linear constraints
# specify maximum likelihood approach:
maximumLikelihoodNormal = function(measurement, pars_init, hs=8, c0=25){
	opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
			"xtol_rel" = 1.0e-8,
			"maxeval" = 1000)
	inequalityConstraints = function(x,measurement,c0,hs){sqrt(c0^2/4+hs^2/4)-(meanDistributionFunction(x[1],x[2])+2*sdDistributionFunction(x[1],x[2]))} # This constraint is used to avoid numerical problems during the optimization. It assumes that mean+2*sd of the true distribution are above the imaging cutoff. If this assumption would not be fulfilled, the determination of true cell sizes would not be reliable anyways.
	return(nloptr(x0=pars_init,eval_f=logLikelihoodNormal,lb=c(c0/2,1),ub=c(10*c0,80),eval_g_ineq=inequalityConstraints,opts=opts,measurement=measurement,hs=hs,c0=c0)) # note: lower bounds (lb) and upper bounds (ub) may have to be adjusted. Specification of a narrow range can improve the numerics of the optimization, but the bounds should not be too narrow in order to not miss the optimal solution.
}

###########################################
# Implementation of the maximum likelihood approach assuming a true underlying Uniform distribution
###########################################

# note: change "distributionFunction", "meanDistributionFunction" and "sdDistributionFunction" to the case of a uniform distribution (example 3, above)
logLikelihoodUniform = function(pars,measurement,hs,c0){
	# function calculating the log-likelihood for given data "measurements" and parameters "pars" based on a uniform distribution
	# Important: "distributionFunction", "meanDistributionFunction" and "sdDistributionFunction" have to be adjusted accordingly!
	a = pars[1]
	b = pars[2]
	logLike = 0
	normFac = normalizationFactor(c0,hs,a,b)
	for (i in 1:length(measurement)){
		temp = measurementDensityNoNorm(measurement[i],hs,c0,a,b)/normFac;
		if (temp<=.Machine$double.xmin){ # smallest non-zero normalized floating-point number - avoid infinite values (log(0) == -Inf)
			logLike= logLike + log(.Machine$double.xmin)
		}else{
			logLike = logLike + log(temp)
		} 
	}
	return(-logLike)
}

library(nloptr) # allows for non-linear constraints
# specify maximum likelihood approach:
maximumLikelihoodUniform = function(measurement, pars_init, hs=8, c0=25){
	opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
			"xtol_rel" = 1.0e-8,
			"maxeval" = 1000)
	inequalityConstraints = function(x,measurement,c0,hs){sqrt(c0^2/4+hs^2/4)-(meanDistributionFunction(x[1],x[2])+2*sdDistributionFunction(x[1],x[2]))} # This constraint is used to avoid numerical problems during the optimization. It assumes that mean+2*sd of the true distribution are above the imaging cutoff. If this assumption would not be fulfilled, the determination of true cell sizes would not be reliable anyways.
	return(nloptr(x0=pars_init,eval_f=logLikelihoodUniform,lb=c(0,sqrt(c0^2/4+hs^2/4)),ub=c(10*c0,100*c0),eval_g_ineq=inequalityConstraints,opts=opts,measurement=measurement,hs=hs,c0=c0)) # note: lower bounds (lb) and upper bounds (ub) may have to be adjusted. Specification of a narrow range can improve the numerics of the optimization, but the bounds should not be too narrow in order to not miss the optimal solution.
}


###########################################
# Calculation of the total number of fat cells in the body based on fat mass (fm) and adipocyte size distribution
###########################################

funToIntCellNum = function(x,...){
	# helper function to calculate mean fat cell volume
	return(distributionFunction(x,...)*x^3*pi*4/3)
}

numCellsCorrected = function(fm,...){
	# Function to calculate total number of fat cells based on fat mass (fm) and the parameters of the assumed true distribution of the cell radii
	Vc = integrate(funToIntCellNum,0,meanDistributionFunction(...)+10*sdDistributionFunction(...),...)$value*1e-15; # last factor is for conversion from micrometer^3 to decimeter^3
	rho = 0.9196 # density of fat in kg/l, i.e. kg/dm^3
	Wc=Vc*rho
	Nc=fm/Wc
	return(Nc)
}
