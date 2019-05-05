\name{testinfo_data_BDI}
\alias{testinfo_data_BDI}
\docType{data}
\title{testinfo_data_BDI}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the Beck Depression Inventory, based on responses from 1709 participants (Olino et al., 2012). 
	The z values are in the first column and the test information values are in the second column.
}
\usage{data(testinfo_data_BDI)}
\source{Olino, T.M., Yu, L., Klein, D.N., Rohde, P., Seeley, J.R., Pilkonis, P.A. & 
	Lewinsohn, P.M._ (2012). Measuring depression using item response theory: An 
	examination of three measures of depressive symptomatology. International Journal 
	of Methods in Psychiatric Research, 21, 76-85.}
\examples{
head(testinfo_data_BDI)

plot(testinfo_data_BDI, type="l", 
     ylim=c(0,(max(testinfo_data_BDI[,2])+1)), 
     xlim=c((min(testinfo_data_BDI[,1])-1),(max(testinfo_data_BDI[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}
