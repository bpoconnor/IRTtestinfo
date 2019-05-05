\name{testinfo_data_SOGSRA}
\alias{testinfo_data_SOGSRA}
\docType{data}
\title{testinfo_data_SOGSRA}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the South Oaks Gambling Screen-Revised for Adolescents, based on responses from 
	981 participants (Chiesi et al., 2012). 
	The z values are in the first column and the test information values are in the second column.
}
\usage{data(testinfo_data_SOGSRA)}
\source{Chiesi, F, Donati, M. A., Galli, S., & Primi, C. (2013). The Suitability of the 
	SOGS-RA as Screening Tool: Item Response Theory-Based evidence. Psychology of 
	Addictive Behaviors, 27(1), 287-293.}
\examples{
head(testinfo_data_SOGSRA)

plot(testinfo_data_SOGSRA, type="l", 
     ylim=c(0,(max(testinfo_data_SOGSRA[,2])+1)), 
     xlim=c((min(testinfo_data_SOGSRA[,1])-1),(max(testinfo_data_SOGSRA[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}

