\name{testinfo_data_AATOCAR5}
\alias{testinfo_data_AATOCAR5}
\docType{data}
\title{testinfo_data_AATOCAR5}
\description{A data frame with IRT test information coefficients and the corresponding z values  
	for the Authority Acceptance scale from the Teacher Observation of Classroom Adaptation-Revised test, 
	based on responses from 8,820 participants (Wu et al., 2012). 
	The z values are in the first column and the test information values are in the second column.
}
\usage{data(testinfo_data_AATOCAR5)}
\source{Wu, J., King, K.M., Witkiewitz, K., Racz, S.J., McMahon, R.J., & Conduct Problems Prevention 
        Research Group. (2012). Item analysis and differential item functioning of a brief conduct problem 
        screen._Psychological Assessment, 24,_444-454.
}
\examples{
head(testinfo_data_AATOCAR5)

plot(testinfo_data_AATOCAR5, type="l", 
     ylim=c(0,(max(testinfo_data_AATOCAR5[,2])+1)), 
     xlim=c((min(testinfo_data_AATOCAR5[,1])-1),(max(testinfo_data_AATOCAR5[,1])+1)), 
     xlab="Latent Trait", ylab="Test Information", 
     main="Test Information Function for the NEO-FFI-A Scale", font.lab=2, cex.lab=1.3)
\dontrun{
}
}