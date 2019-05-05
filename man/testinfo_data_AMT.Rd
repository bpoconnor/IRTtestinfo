\name{testinfo_data_AMT}
\alias{testinfo_data_AMT}
\docType{data}
\title{testinfo_data_AMT}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the Autobiographical Memory Test, based on responses from 194 participants (Monahan et al., 2001). 
	The z values are in the first column and the test information values are in the second column.
}
\usage{data(testinfo_data_AMT)}
\source{Griffith, J. W., Kleim, B., Sumner, J. A., & Ehlers, A. (2012). The factor structure of the Autobiographical Memory Test in recent trauma survivors. Psychological Assessment, 24, 640- 646.}
\examples{
head(testinfo_data_AMT)

plot(testinfo_data_AMT, type="l", 
     ylim=c(0,(max(testinfo_data_AMT[,2])+1)), 
     xlim=c((min(testinfo_data_AMT[,1])-1),(max(testinfo_data_AMT[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}
