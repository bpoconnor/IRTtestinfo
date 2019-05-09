

constant_info <- function(test1info, test2info, 
                        Nsets = 100, 
                        Ncases = 1000,
                        rTruePop = .95, 
                        distribMN1 = 0, distribMN2 = 0, distribSD1 = 1, distribSD2 = 1 ) {

# This code reveals the attenuation in effect sizes that occurs for combinations of
# constant test information values for two measures.

# test1info = provide a set or range of information values for test 1 

# test2info = provide a set or range of information values for test 2 

# Nsets = specify the number of random data sets to be used in the analyses

# Ncases = specify the # of cases for each random data set

# rTruePop = set the true, population correlation to be used in the analyses

# distribMN1 = specify the Measure 1 mean to be used for the random data sets

# distribMN2 = specify the Measure 2 mean to be used for the random data sets

# distribSD1 = specify the Measure 1 standard deviation to be used for the random data sets

# distribSD2 = specify the Measure 2 standard deviation to be used for the random data sets



##################### end of user specifications #########################################

pb = txtProgressBar(min = 0, max = length(test1info), initial = 0, style = 3) 

results2 <- matrix(-9999,length(test1info),length(test2info))

for (lupe1 in 1:length(test1info)) {
	
	setTxtProgressBar(pb,lupe1)

	for (lupe2 in 1:length(test2info)) {
		results1 <- matrix(-9999,Nsets,1)

		for (lupe3 in 1:Nsets) {

			# generate raw data for 2 variables with a specific correlation
			randdata <- data.frame(mvrnorm(n <- Ncases,mu=c(distribMN1,distribMN2),
		            Sigma=matrix(c(distribSD1,rTruePop,rTruePop,distribSD2),nrow=2),empirical=TRUE))
			Xr <- randdata[,1]
			Yr <- randdata[,2]

			XYcorrT = cor(Xr, Yr)  # XYcorrT

			# the mean SEM
			semMNx <- 1 / sqrt(test1info[lupe1])
			semMNy <- 1 / sqrt(test2info[lupe2])

			# add meas error: from  2007 Zimmerman p 924
			eX = t(t(( rnorm(1000, mean = 0, sd = semMNx) )))  
			eY = t(t(( rnorm(1000, mean = 0, sd = semMNy) ))) 

			Xmnsem = Xr + eX 
			Ymnsem = Yr + eY 

			# correlation between the scores with measurement error, based on one SEM per variable
			corrMNSEM = cor( Xmnsem, Ymnsem )

			PR <- 1 - (corrMNSEM / rTruePop)

			results1[lupe3,] = PR 
		}

	results2[lupe1,lupe2] <- mean(results1)
	}
}



# for colors -- from visreg
pal <- function(n, alpha=1) {
  if (n==2) {
    val <- hcl(seq(15,375,len=4), l=60, c=150, alpha=alpha)[c(1,3)]
  } else val <- hcl(seq(15,375,len=n+1), l=60, c=150, alpha=alpha)[1:n]
  val
}
color.palette=colorRampPalette(c(pal(3)[3],"gray90",pal(3)[1]),space="Lab")


# heat map
print(
levelplot( as.matrix(results2),row.values=test1info,column.values=test2info, cuts = 60, 
  main="Proportionate Reductions in Effect Size\nAcross Levels of Test Information",
  xlab="Test Information for Measure 1",
  ylab="Test Information for Measure 2",  
  col.regions=color.palette)
)


colnames(results2) <- c(test2info)
rownames(results2) <- c(test1info)

CONSTANToutput <- list(plotdata = results2)

return(invisible(CONSTANToutput))

}







