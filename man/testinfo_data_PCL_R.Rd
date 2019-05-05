\name{testinfo_data_PCL_R}
\alias{testinfo_data_PCL_R}
\docType{data}
\title{testinfo_data_PCL_R}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the Psychopathy Checklist-Revised (PCL-R), based on responses from 715 participants 
	(Cooke et al., 2001). The z values are in the first column and the test information values 
	are in the second column.
}
\usage{data(testinfo_data_PCL_R)}
\source{Cooke, D. J., Kosson, D. S. & Michie, C. (2001)_Psychopathy and ethnicity: 
	structural, item and test generalizability of the Psychopathy Checklist Revised 
	(PCL-R) in Caucasian and African-American participants._Psychological Assessment,_13,_531- 542.}
\examples{
head(testinfo_data_PCL_R)

plot(testinfo_data_PCL_R, type="l", 
     ylim=c(0,(max(testinfo_data_PCL_R[,2])+1)), 
     xlim=c((min(testinfo_data_PCL_R[,1])-1),(max(testinfo_data_PCL_R[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}
