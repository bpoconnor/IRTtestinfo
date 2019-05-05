


fluctuating <- function(plotdata1, plotdata2, rawdata1=NULL, rawdata2=NULL, rTruePop = .95, rObserved = NULL,
                        Nsets = 100, distribMN1 = 0, distribMN2 = 0, distribSD1 = 1, distribSD2 = 1 ) {


cat('\n\n\n\nIRT attenuation:\n')

# This code estimates the level of attenuation in correlation that occurs for two variables, 
# given the IRT test information values for the two measures. The degree of attenuation is determined
# empirically rather than by formulas. The program preserves and takes into account 
# the distributions of the user’s sample-specific raw data values. The program answers these two questions: 
# When there are two normally-distributed correlated variables in a population, how much attenuation 
# due to measurement error is caused by the fluctuating levels of test information for the variables 
# and by the sample-specific variable distributions? 
# What is the corrected-for-attenuation estimate of the correlation between the two variables?

# read in the test information & corresponding z values for Measure 1
# the z values should be the first column of scores, & 
# the information values should be the second column of scores

# plotdata1 = the test information & corresponding z values for Measure 1     

# plotdata2 = the test information & corresponding z values for Measure 2     

# rawdata1 = the raw data z scores for Measure 1    

# rawdata2 = the raw data z scores for Measure 2     

# rTruePop = set the true, population correlation to be used in the analyses

# rObserved = an observed population correlation

# Nsets = specify the number of random data sets to be used in the analyses

# distribMN1 = specify the Measure 1 mean to be used for the random data sets

# distribMN2 = specify the Measure 2 mean to be used for the random data sets

# distribSD1 = specify the Measure 1 standard deviation to be used for the random data sets

# distribSD2 = specify the Measure 2 standard deviation to be used for the random data sets



# check on the variable distributions
# v1MNdiff <- abs(distribMN1 - mean(rawdata1))
# v2MNdiff <- abs(distribMN2 - mean(rawdata2))
# v1SDdiff <- abs(distribSD1 - sd(rawdata1))
# v2SDdiff <- abs(distribSD2 - sd(rawdata2))
# ntest1 <- shapiro.test(rawdata1)
# ntest1$p.value > .05
# ntest2 <- shapiro.test(rawdata2)
# ntest2$p.value > .05
varDresults <- matrix(-9999,2,4)
varDresults[1,] <- c(distribMN1, mean(rawdata1), distribSD1, sd(rawdata1))
varDresults[2,] <- c(distribMN2, mean(rawdata2), distribSD2, sd(rawdata2))
colnames(varDresults) <- c('set distribMN','   Mean of rawdata','        set distribSD','    SD of rawdata')
rownames(varDresults) <- c('rawdata1','rawdata2')
cat('\n\nCheck the variable distributions:\n\n')
print(round(varDresults,2))
cat('\n\nThe means and standard deviations of the raw data should be very similar')
cat('\nto the set values for distribMN and distribSD. The default set values are 0 and 1.')
cat('\nThe rawdata scores should be standardized. If the distributions of the rawdata')
cat('\nare different from the set values, then errors may result.\n')


if (is.null(rawdata1) | is.null(rawdata2)) {
	
	if (is.null(rawdata1) == F) Ncases = length(rawdata1)
	if (is.null(rawdata2) == F) Ncases = length(rawdata2)
	if (is.null(rawdata1) & is.null(rawdata2)) Ncases = 1000
	
	if (is.null(rawdata1)) rawdata1 <- rnorm(n = Ncases, mean = distribMN1, sd = distribSD1) 
	if (is.null(rawdata2)) rawdata2 <- rnorm(n = Ncases, mean = distribMN2, sd = distribSD2) 
}


# setting any negative information values to 0
plotdata1[,2] <- pmax(plotdata1[,2],0.001)
plotdata2[,2] <- pmax(plotdata2[,2],0.001)

# computing SEM values	
SEMx <- 1 / sqrt(plotdata1[,2])
SEMy <- 1 / sqrt(plotdata2[,2])
		
# grand loop
results <- matrix(-9999,Nsets,4)
#pb = txtProgressBar(min = 0, max = Nsets, initial = 0, style = 3) 
for (lupe in 1:Nsets) {
	
#setTxtProgressBar(pb,lupe)

# generate raw data for 2 variables with a specific correlation
randdata <- data.frame(mvrnorm(n <- 1000,mu=c(distribMN1,distribMN2),
		    Sigma=matrix(c(distribSD1,rTruePop,rTruePop,distribSD2),nrow=2),empirical=TRUE))
Xr <- randdata[,1]
Yr <- randdata[,2]
	

# Go through the two lists of real-data values and extract from the random 
# variables the scores that are at the same z-value locations as the real, raw data scores.
indicesXz <- as.matrix( findInterval(rawdata1, sort(Xr), all.inside=TRUE) )
indicesYz <- as.matrix( findInterval(rawdata2, sort(Yr), all.inside=TRUE) )

extracted1 <- extracted2 <- matrix(-9999,length(rawdata2),2)

for (lupe2 in 1:length(indicesXz)) {	
	extracted1[lupe2,] <- as.matrix(randdata[indicesXz[lupe2],])
	extracted2[lupe2,] <- as.matrix(randdata[indicesYz[lupe2],])
}

extractedX1 <- extracted1[,1]
extractedY1 <- extracted1[,2]

extractedX2 <- extracted2[,1]
extractedY2 <- extracted2[,2]

results[lupe,1:2] <- c(cor(extractedX1, extractedY1), cor(extractedX2, extractedY2))


# For each extracted score value for each variable, identify the z value location on the
# latent trait continuum, find the corresponding value from the test information function,
# compute the location-specific SEM value, and then add this level of error to the extracted score.

indicesX1 <- as.matrix( findInterval(extractedX1, plotdata1[,1], all.inside=TRUE) )
indicesY1 <- as.matrix( findInterval(extractedY1, plotdata2[,1], all.inside=TRUE) )

indicesX2 <- as.matrix( findInterval(extractedX2, plotdata1[,1], all.inside=TRUE) )
indicesY2 <- as.matrix( findInterval(extractedY2, plotdata2[,1], all.inside=TRUE) )

SEMXextracted1 <- SEMYextracted1 <- matrix(-9999, length(extractedX1), 1) 
SEMXextracted2 <- SEMYextracted2 <- matrix(-9999, length(extractedX2), 1) 

for (luper in 1:length(extractedX1) ) {
	SEMXextracted1[luper,1] <- SEMx[indicesX1[luper,1]]
	SEMYextracted1[luper,1] <- SEMy[indicesY1[luper,1]]

	SEMXextracted2[luper,1] <- SEMx[indicesX2[luper,1]]
	SEMYextracted2[luper,1] <- SEMy[indicesY2[luper,1]]
}
		
extracted_error_X1 <- extractedX1 + as.matrix( rnorm(length(extractedX1), mean = 0, sd = SEMXextracted1), nrow=length(extractedX1) ) 
extracted_error_Y1 <- extractedY1 + as.matrix( rnorm(length(extractedY1), mean = 0, sd = SEMYextracted1), nrow=length(extractedY1) ) 
XYcorrESEM1 <- cor( extracted_error_X1, extracted_error_Y1 )

extracted_error_X2 <- extractedX2 + as.matrix( rnorm(length(extractedX2), mean = 0, sd = SEMXextracted2), nrow=length(extractedX2) ) 
extracted_error_Y2 <- extractedY2 + as.matrix( rnorm(length(extractedY2), mean = 0, sd = SEMYextracted2), nrow=length(extractedY2) ) 
XYcorrESEM2 <- cor( extracted_error_X2, extracted_error_Y2 )
	
results[lupe,3:4] <- c(XYcorrESEM1, XYcorrESEM2)

}	

cormeans <- colMeans(results)

propred1 <- 1 - (cormeans[3] / rTruePop)  # proportionate reduction in correlation
propred2 <- 1 - (cormeans[4] / rTruePop)  # proportionate reduction in correlation


cat('\n\n\nThe specified population correlation for the analyses = ', round(rTruePop,2),
    '\n\nThe specified number of random datasets for the analyses = ', Nsets,
    '\n',
    '\n\nThe mean random data correlation before the addition of measurement error when using the raw data distribution for Measure 1 =', round(cormeans[1],2),
    '\n\nThe mean random data correlation before the addition of measurement error when using the raw data distribution for Measure 2 =', round(cormeans[2],2),
    '\n',
    '\n\nThe mean random data correlation after the addition of measurement error when using the raw data distribution for Measure 1 =', round(cormeans[3],2),
    '\n\nThe mean random data correlation after the addition of measurement error when using the raw data distribution for Measure 2 =', round(cormeans[4],2),
    '\n',
    '\n\nThe proportionate reduction in correlation when using the raw data distribution for Measure 1 =', round(propred1,2),
    '\n\nThe proportionate reduction in correlation when using the raw data distribution for Measure 2 =', round(propred2,2) )    


if (is.null(rawdata1) == F & is.null(rawdata2) == F) rObserved <- cor(rawdata1,rawdata2)

if (is.null(rObserved) == F) {
	correctedR1 <- rObserved / (1 - propred1)
	correctedR2 <- rObserved / (1 - propred2)

	cat('\n\nThe computed, or provided, correlation between the two measures =', round(rObserved,2),
	    '\n\nThe corrected-for-attenuation correlation when using the raw data distribution for Measure 1 =', round(correctedR1,2),
	    '\n\nThe corrected-for-attenuation correlation when using the raw data distribution for Measure 2 =', round(correctedR2,2) )
}

cat('\n\n\nThe analyses were conducted while preserving the variable distributions in the raw data.',
    '\nThe program does so while first preserving the distribution for Measure 1, and then again',
    '\nwhile first preserving the distribution for Measure 2. The above results for Measures 1 & 2',
    '\nshould be highly similar. If they are not highly similar, then re-run the analyses',
    '\nusing a larger number of random data sets. The findings for the 2 measures will eventually converge.',
    '\n', 
    '\nTip: Set the true, population correlation to be used in the analyses (rTruePop) to a',
    '\nhigh value, e.g., .95. Fewer random data sets will be then required for the findings for',
    '\nthe two measures to converge. The results will be identical to those for lower rTruePop values,', 
    '\nwhich would require more random dataset processing before convergence.\n\n')


IRToutput <- list(rTruePop=rTruePop, propred1=propred1, propred2=propred2,
                  correctedR1=correctedR1, correctedR2=correctedR2)

return(invisible(IRToutput))

}



