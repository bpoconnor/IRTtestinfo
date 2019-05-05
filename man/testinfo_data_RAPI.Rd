\name{testinfo_data_RAPI}
\alias{testinfo_data_RAPI}
\docType{data}
\title{testinfo_data_RAPI}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the Rutgers Alcohol Problem Index, based on responses from 895 participants (Neal et al., 2006). 
	The z values are in the first column and the test information values are in the second column.
}
\usage{data(testinfo_data_RAPI)}
\source{Neal, D. J., Corbin, W. R., & Fromme, K. (2006). Measurement of alcohol-related consequences 
	among high school and college students: Application of item response models to the Rutgers 
	Alcohol Problem Index. Psychological Assessment, 18, 402_414.}
\examples{
head(testinfo_data_RAPI)

plot(testinfo_data_RAPI, type="l", 
     ylim=c(0,(max(testinfo_data_RAPI[,2])+1)), 
     xlim=c((min(testinfo_data_RAPI[,1])-1),(max(testinfo_data_RAPI[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}
