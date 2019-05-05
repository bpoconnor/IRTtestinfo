\name{testinfo_data_NEO_FFI_A}
\alias{testinfo_data_NEO_FFI_A}
\docType{data}
\title{testinfo_data_NEO_FFI_A}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the NEO-FFI Antagonism-Agreeableness scale, based on responses from 1229 participants 
	in the MacArthur Violence Risk Assessment Study dataset (Monahan et al., 2001). 
	The z values are in the first column and the test information values are in the second column.
}
\usage{data(testinfo_data_NEO_FFI_A)}
\source{Monahan, J., Steadman, H., Silver, E., Appelbaum, P., Robbins, P., Mulvey, E., 
	Roth, L., Grisso, T., & Banks, S. (2001). Rethinking Risk Assessment: The MacArthur 
	Study of Mental Disorder and Violence. New York: Oxford University Press.}
\examples{
head(testinfo_data_NEO_FFI_A)

plot(testinfo_data_NEO_FFI_A, type="l", 
     ylim=c(0,(max(testinfo_data_NEO_FFI_A[,2])+1)), 
     xlim=c((min(testinfo_data_NEO_FFI_A[,1])-1),(max(testinfo_data_NEO_FFI_A[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}
