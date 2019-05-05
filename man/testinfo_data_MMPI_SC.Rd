\name{testinfo_data_MMPI_SC}
\alias{testinfo_data_MMPI_SC}
\docType{data}
\title{testinfo_data_MMPI_SC}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the MMPI-2-RF Somatic Complaints (RC1) Scale, based on responses from 399 
	participants (Thomas et al., 2010). The z values are in the first column and the test 
	information values are in the second column.
}
\usage{data(testinfo_data_MMPI_SC)}
\source{Thomas, M. L., & Locke, D. E. C. (2010). Psychometric properties of the MMPI-2-RF 
	Somatic Complaints (RC1) scale._Psychological Assessment,_22, 492-503.}
\examples{
head(testinfo_data_MMPI_SC)

plot(testinfo_data_MMPI_SC, type="l", 
     ylim=c(0,(max(testinfo_data_MMPI_SC[,2])+1)), 
     xlim=c((min(testinfo_data_MMPI_SC[,1])-1),(max(testinfo_data_MMPI_SC[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}
