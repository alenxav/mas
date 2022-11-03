\name{mas-package}
\alias{mas-package}
\alias{mas}
\docType{package}
\title{\packageTitle{mas}}
\description{\packageDescription{mas}}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{mas}
\packageIndices{mas}
}
\author{
\packageAuthor{mas}
Maintainer: \packageMaintainer{mas}
}
\references{

Wei, J. and Xu, S., 2016. A random-model approach to QTL mapping in multiparent advanced generation intercross (MAGIC) populations. Genetics, 202(2), pp.471-486.

Xavier, A. and Habier, D., 2022. A new approach fits multivariate genomic prediction models efficiently. Genetics Selection Evolution, 54(1), pp.1-15.

}
\examples{

# load the toy dataset
data( soy ) 

# run gwas
fit = GWAS(y, Z, pop) 
# adjust variances
fit2 = CorrectBeavis( fit ) 

# compare plots
par(mfrow=c(1,2))
plot( fit,  h2=TRUE ) # display QTL h2
plot( fit2, h2=TRUE ) # display adjusted QTL h2


}
