\name{testinfo_data_PICTo}
\alias{testinfo_data_PICTo}
\docType{data}
\title{testinfo_data_PICTo}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the Psychological Inventory of Criminal Thinking Styles, based on responses from 
	2872 participants (Walters et al., 2011). 
	The z values are in the first column and the test information values are in the second column.
}
\usage{data(testinfo_data_PICTo)}
\source{Walters, G. D., Hagman, B. T., & Cohn, A. M. (2011). Toward a hierarchical model of 
	criminal thinking: Evidence from item response theory and confirmatory factor analysis. 
	Psychological Assessment, 23, 925-936.}
\examples{
head(testinfo_data_PICTo)

plot(testinfo_data_PICTo, type="l", 
     ylim=c(0,(max(testinfo_data_PICTo[,2])+1)), 
     xlim=c((min(testinfo_data_PICTo[,1])-1),(max(testinfo_data_PICTo[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}
