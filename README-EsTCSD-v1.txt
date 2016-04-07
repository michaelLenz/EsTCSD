README: EsTCSD - Estimate Tissue Cell Size/Type Distribution - Version 1.1

##########################################

General information:

Determination of tissue cell size distributions from cross-section images is subject to different types of bias due to underestimation of the number of small cells as well as the size of measured cells. Furthermore, imaging parameters and slice thickness do have an effect on the measurements (see the EsTCSDpublic ation for more details).
EsTCSD is an algorithm for the adjustment of these types of bias for cells that are (close to) spherical in shape, e.g. adipocytes. The algorithm can be used to either calculate the expected measured cell size distribution given the true distribution, or for solving the inverse problem of determining the adjusted/corrected distribution based on cross-sectional measurements. The inverse approach uses a maximum likelihood approach together with an assumption about the type of distribution (e.g. gamma distribution or normal distribution) and determines the corresponding optimal parameters of that distribution.

##########################################

How to use:

In order to perform the cell size correction, you need to provide the following information:
1) A vector of measured cell diameters (should be at least a few hundred, the unit is typically in micro meter)
2) The type of distribution you expect the corrected cell sizes to follow (examples for gamma, normal, and uniform distributions are provided)
3) The imaging cutoff parameter, i.e. the minimum diameter of particles on the microscopy image to be considered as cells of the desired type (e.g. adipocytes)
4) The slice thickness, i.e. the height of the slice that has been cut of the biopsy

Furthermore, you have to make sure that the cells can be considered spherical in shape (e.g. adipocytes)

In addition, in order to calculate the total number of fat cells (adipocytes) in a body, you have to provide a measurement of total fat mass (in kg 
if cell diameters are provided in micro meter).

How to perform the correction:
Load vector of measured cell sizes into R (should be a plain vector of numbers in R, can be read in from various files(.csv, .txt or .xlsx) into R)
Specify the type of distribution in the function "distributionFunction" and specify the functions "meanDistributionFunction" and "sdDistributionFunction" to calculate mean and standard deviations for this type of distribution (or use one of the provided examples - gamma distribution, normal distribution, uniform distribution) - see file EstCSDcalculations.R
Specify log-likelihood (e.g. as in "logLikelihoodGamma") and maximum-likelihood (e.g. as in maximumLikelihoodGamma) functions for the specific underlying true distribution function (or use one of the provided examples) - see file EstCSDcalculations.R
Determine the maximum likelihood estimates using the specified maximum-likelihood function with parameters:
	"measurement": the vector of measured cell sizes
	"pars_init": a vector of initial parameters for the "distributionFunction" function
	"hs": the slice thickness
	"c0": the imaging cutoff parameter 
Potentially (if desired) determine the total number of cells using the provided "numCellsCorrected" function with the fat mass (fm) and the maximum likelihood parameters as input.

Hint: The EsTCSDexampleUsage.R file provides helpful examples on how to use the EsTCSD code.

Note that the implementation in the current version is focused on the application to human adipocyte cell sizes using an underlying gamma distribution assumption.
Applications to other cases might need slight adjustments (e.g. specification of a different distribution function etc). Please feel free to contact us if you need help for a specific application (see contact information below).

##########################################

Description of files included:

EsTCSDexampleUsage.R - Example on how to use the EsTCSD R script
EstCSDcalculations.R - Functions performing the EsTCSD correction
cellDiameters.Rdata  - Example of measured cell diameters (one biopsy, 419 cells measured)

##########################################

Citation: Lenz M, Roumans NJT, Vink RG, van Baak M, Mariman E, Arts ICW, de Kok TM, Ertaylan G
	  "Estimating real cell size distribution from cross section microscopy imaging", submitted.

Contact: Michael Lenz (michael.lenz@maastrichtuniversity.nl) or 
         Gokhan Ertaylan (gokhan.ertaylan@maastrichtuniversity.nl)
	 Maastricht Centre for Systems Biology, Maastricht University, The Netherlands (macsbio.maastrichtuniversity.nl)
