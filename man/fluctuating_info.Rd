\name{fluctuating_info}
\alias{fluctuating_info}
\title{Effect size attenuation for fluctuating levels of IRT test information}
\description{This function estimates the level of attenuation in correlation that 
	occurs for two variables, given the IRT test information functions for the two 
	measures (following O'Connor, 2018a, 2018b). The degree of attenuation is 
	determined empirically rather than by formulas. The function preserves and takes 
	into account the distributions of the user's sample-specific raw data values,
	when such values are provided. The function answers these two questions: When 
	there are two normally-distributed correlated variables in a population, how much 
	attenuation due to measurement error is caused by the fluctuating levels of test 
	information for the variables? What is the corrected-for-attenuation estimate of 
	the correlation between the two variables?
}
\usage{fluctuating_info(plotdata1, plotdata2, rawdata1=NULL, rawdata2=NULL,
                   rTruePop = .95, rObserved = NULL, Nsets = 100,
                   distribMN1 = 0, distribMN2 = 0, distribSD1 = 1, distribSD2 = 1)
}
\arguments{	
  \item{plotdata1}{
  \code{}(required) A dataframe with the test information & corresponding z values for Measure 1.
         The z values should be the first column of scores & the information values
         should be the second column of scores.}

  \item{plotdata2}{
  \code{}(required) A dataframe with the test information & corresponding z values for Measure 2.
         The z values should be the first column of scores & the information values
         should be the second column of scores.}

  \item{rawdata1}{
  \code{}(optional) A dataframe with the respondents' z scores for Measure 1.}

  \item{rawdata2}{
  \code{}(optional) A dataframe with the respondents' z scores for Measure 2.}

  \item{rTruePop}{
  \code{}(optional) The true, population correlation to be used in the analyses. The default = .95.}

  \item{rObserved}{
  \code{}(optional) An observed sample correlation.}

  \item{Nsets}{
  \code{}(optional) The number of random data sets to be used in the analyses. The default = 100.}

  \item{distribMN1}{
  \code{}(optional) The Measure 1 mean to be used for the random data sets. The default = 0.}

  \item{distribMN2}{
  \code{}(optional) The Measure 2 mean to be used for the random data sets. The default = 0.}

  \item{distribSD1}{
  \code{}(optional) The Measure 1 standard deviation to be used for the random data sets. The default = 1.}

  \item{distribSD2}{
  \code{}(optional) The Measure 2 standard deviation to be used for the random data sets. The default = 1.}
}

\details{
The analyses are conducted while preserving the variable distributions in the raw data.
The function does so while first preserving the distribution for Measure 1, and then again
while first preserving the distribution for Measure 2. The results for Measures 1 & 2
should be highly similar. If they are not highly similar, then re-run the analyses
using a larger number of random data sets. The findings for the 2 measures will eventually converge.
Tip: Set the true, population correlation to be used in the analyses (rTruePop) to a
high value, e.g., .95. Fewer random data sets will be required for the findings for
the two measures to converge. The results will be identical to those for lower rTruePop values,
which would require more random dataset processing before convergence.

The analyses are conducted on the z values that correspond with the scale raw scores. 
The "rawdata" values should be standardized. When standardizing the scores, it is generally 
better to use the means and SDs from normative data for each 
measure, i.e., rather than using the (usually smaller) sample means and SDs.

If both rawdata1 and rawdata2 are provided, then the function will compute the observed correlation
between the two variables. The output will include an estimate of the true population correlation,
given the test information functions, the variable distributions, and the observed correlation.

If only one of rawdata1 or rawdata2 is provided, then the function will conduct the analyses
using the distribution of the provided data along with random normal generated data, for the
same number of cases, for the other variable. If a value for rObserved is not provided, then 
the function will return the proportionate reduction in effect size. If a value for rObserved 
is provided, then the output will include an estimate of the true population correlation, 
given the test information functions, the one provided variable distribution, and the provided 
observed correlation.

If neither rawdata1 or rawdata2 is provided, then the function will conduct the analyses
using random normal generated data for both variables, for 1000 cases.
If a value for rObserved is not provided, then the function will return the proportionate
reduction in effect size. If a value for rObserved is provided, then the
output will include an estimate of the true population correlation, given the test information 
functions and the provided observed correlation.
}

\value{The displayed output includes the set, population correlation for the analyses,
the mean random data correlation before the addition of measurement error when using the 
raw data distribution for Measure 1; the mean random data correlation before the addition of 
measurement error when using the raw data distribution for Measure 2; the mean random data 
correlation after the addition of measurement error when using the raw data distribution 
for Measure 1; the mean random data correlation after the addition of measurement error 
when using the raw data distribution for Measure 2; the proportionate reduction in correlation 
when using the raw data distribution for Measure 1; the proportionate reduction in correlation 
when using the raw data distribution for Measure 2; the raw data correlation between the two 
measures; the corrected-for-attenuation correlation when using the raw data distribution for 
Measure 1; and the corrected-for-attenuation correlation when using the raw data distribution 
for Measure 2. 
	   
	   The returned output is a list with the following elements:
  \item{rTruePop}{the set, population correlation for the analysess}
  \item{propred1}{the proportionate reduction in correlation when using the raw data distribution for Measure 1}
  \item{propred2}{the proportionate reduction in correlation when using the raw data distribution for Measure 2}
  \item{correctedR1}{the corrected-for-attenuation correlation when using the raw data distribution for Measure 1}
  \item{correctedR2}{the corrected-for-attenuation correlation when using the raw data distribution for Measure 2}
}
\references{
       {O'Connor, B. P. (2018a). Clarifications regarding test information and reliability, and new methods for estimating attenuation due to measurement error: Reply to Markon (2018). Psychological Assessment, 30(8), 1010-1012.}
\cr\cr {O'Connor, B. P. (2018b). An illustration of the effects of fluctuations in test information on measurement error, the attenuation of effect sizes, and diagnostic reliability. Psychological Assessment, 30(8), 991-1003.}
}
\author{Brian P. O'Connor }
\examples{
# when both test information functions and raw data are provided
fluctuating_info(
    plotdata1 = testinfo_data_PCL_SV, 
    plotdata2 = testinfo_data_NEO_FFI_A,
    rawdata1 = zscores_data_PCL_SV_NEO_FFI_A[,1], 
    rawdata2 = zscores_data_PCL_SV_NEO_FFI_A[,2], 
    rTruePop = .95, 
    Nsets = 100,
    distribMN1 = 0, 
    distribMN2 = 0, 
    distribSD1 = 1, 
    distribSD2 = 1 )

# when information functions but no raw data are provided
fluctuating_info(
    plotdata1 = testinfo_data_PCL_SV, 
    plotdata2 = testinfo_data_NEO_FFI_A,
    rTruePop = .95, 
    rObserved = .25, 
    distribMN1 = 0, 
    distribMN2 = 0, 
    distribSD1 = 1, 
    distribSD2 = 1 )


\dontrun{

# # the next, example demonstrates how to obtain and use test information values 
# # from the grm function (ltm package) and from the mirt function (mirt package)
# library(ltm)
# library(mirt)

# # obtain the test information and z values using the grm function (ltm package)
# # for the Science data in the ltm package
# # run the IRT analyses using grm
# test1results_<- grm(ltm::Science[c(1,3,4,7)])
# # display the test information curve
# plot(test1results, type = "IIC", z = seq(-4, 4, length = 100), items = 0) _
# # obtain the test information values & z scores
# test1info <- plot(test1results, type = "IIC", z = seq(-4, 4, length = 1000), items = 0)   
# head(test1info); dim(test1info); summary(test1info) 
# # The z values are in the first column & the test information 
# # values are in the second column.

         
# # obtain the test information and z values using the mirt function (mirt package)
# # for generated data as in the mirt package example
# set.seed(12345)
# a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
# d <- matrix(rnorm(15,0,.7),ncol=1)
# d <- cbind(d, d-1, d-2)
# itemtype <- rep('graded', nrow(a))
# N <- 1000
# dataset2 <- simdata(a, d, N, itemtype)
# # run the IRT analyses using mirt
# test2results <- mirt(dataset2, 1, itemtype = "graded", SE=TRUE)
# # display the test information curve
# plot(test2results, type ="info", theta_lim = c(-4,4))
# # obtain the test information values & z scores
# Theta <- matrix(seq(-4,4,.01))
# test2info <- cbind(Theta, testinfo(test2results, Theta))
# head(test2info); dim(test2info); summary(test2info) 
# # The z values are in the first column & the test information 
# # values are in the second column.


# fluctuating_info(
# 	plotdata1 = test1info, 
# 	plotdata2 = test2info,
# 	rTruePop = .95, 
# 	rObserved = .25, 
# 	Nsets = 100,
# 	distribMN1 = 0, 
# 	distribMN2 = 0, 
# 	distribSD1 = 1, 
# 	distribSD2 = 1 )
}

}


