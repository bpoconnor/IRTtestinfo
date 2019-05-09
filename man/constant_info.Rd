\name{constant_info}
\alias{constant_info}
\title{Effect size attenuation for constant levels of IRT test information}
\description{What constitutes an adequate level of IRT test information? For example, at what 
	levels of test information do effect size attenuation problems occur? This function 
	displays the proportionate reductions in the magnitudes of correlation coefficients that 
	occur for two measures with constant test information values (following 
	O'Connor, 2018a, 2018b). The plot data are produced 
	by converting the test information values into SEM values, adding these degrees of location-specific
	measurement error to two random normal variables (e.g., N = 1000) which has a specific, 
	previously set correlation, and then computing the observed correlation for the altered 
	variables. The procedure is repeated Nsets times for each value of test information. The 
	means of the correlations are used as the data points for the plot. The proportionate 
	reductions in the effect sizes are the same across levels of true population correlation. The function 
	permits users to specify alternative values for the sampling variances, such as the 
	values they may have for their own samples, and to observe the consequences for effect 
	size attenuation.
}
\usage{constant_info(test1info, 
              test2info, 
              Nsets = 100, 
              Ncases = 1000,
              rTruePop = .95,
              distribMN1 = 0, distribMN2 = 0, distribSD1 = 1, distribSD2 = 1)
}
\arguments{	

  \item{test1info}{
  \code{}(required) A set or range of information values for test 1.}

  \item{test2info}{
  \code{}(required) A set or range of information values for test 2.}

  \item{Nsets}{
  \code{}(optional) The number of random data sets to be used in the analyses. The default = 100.}

  \item{Ncases}{
  \code{}(optional) The number of cases for each random data set. The default = 1000.}

  \item{rTruePop}{
  \code{}(optional) The true, population correlation to be used in the analyses. The default = .95.}

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
It is best to first run the function using few values for test1info and test12info, in order
to become familiar with how the function works. For example, try seq(1:5) for both.
For smoother and more precise plots, use more numerous, finer values, such as seq(1, 10, by = .1),
which will take much longer to produce results.

Set the true, population correlation to be used in the analyses (rTruePop) to a
high value, e.g., .95. Fewer random data sets will be required for the findings
to converge. The results will be identical to those for lower rTruePop values.
}

\value{The displayed output is a heatmap. 
	   The returned output is a list with the plot data.
}
\references{
       {O'Connor, B. P. (2018a). Clarifications regarding test information and reliability, and new methods for estimating attenuation due to measurement error: Reply to Markon (2018). Psychological Assessment, 30(8), 1010-1012.}
\cr\cr {O'Connor, B. P. (2018b). An illustration of the effects of fluctuations in test information on measurement error, the attenuation of effect sizes, and diagnostic reliability. Psychological Assessment, 30(8), 991-1003.}
}
\author{Brian P. O'Connor}
\examples{
constant_info(
   test1info <- c(1,2,3,4,5), 
   test2info <- c(1,2,3,4,5), 
   Nsets = 100, 
   Ncases = 1000,
   rTruePop = .95, 
   distribMN1 = 0, 
   distribMN2 = 0, 
   distribSD1 = 1, 
   distribSD2 = 1 )
}
