\name{gwas}
\alias{GWAS}
\alias{MLM}
\alias{CorrectBeavis}
\alias{print.GenModel}
\alias{plot.GenModel}

\title{
 Genome-Wide Association Studies
}
\description{

Genome-wide association analysis that allows for integrating population membership into the analysis and information from correlated traits.

}

\usage{
GWAS(Y,GEN,M,maxit=500,logtol=-8,cores=1)
CorrectBeavis(fit)
}

\arguments{
  \item{Y}{Numeric matrix with phenotypic information. Columns are traits, rows are individuals. NA are allowed.}
  \item{GEN}{Numeric matrix of genotypic information. Columns are markers, rows are individuals. NA are not allowed.}
  \item{M}{Numeric matrix with membership information. Columns are population classes, rows are individuals. NA are not allowed.}
  \item{maxit}{Integer. Maximum number of iterations.}
  \item{logtol}{Scalar. Convergerce parameter of null-model in log scale.}
  \item{cores}{Integer. Number of CPUs to use OpenMP.}
  \item{fit}{Object of class GenModel, output of the GWAS function.}
}

\details{


The model for the genome-wide association model is set as follows:

\deqn{ y = Xb + Zu + g + e }

where \eqn{y} is the response variable, \eqn{Xb} corresponds to the fixed effect term set as the membership matrix, \eqn{Zu} corresponds to the marker-membership interaction term, \eqn{g} is the polygenic term defines by \eqn{g=Ma}, where \eqn{M} is the genotypic matrix and \eqn{g} are marker effects, and \eqn{e} is the residual term.

The null-hypothesis term cosists of a similar model, but without \eqn{Zu}. 

For populations without knwon membership, set \code{M=matrix(rep(1,nrow(Y)))}.

Theorical description of the model is provided by Wei and Xu (2016). Variance components are REML estimates under the univariate case. In the multivariate case, variance and covariance components are obtained through a general-purpose REML approximation (Van Raden and Jung 1987), which yields the exact same solution as REML when all traits are observed in all individuals. The multivariate generalization is described by Xavier and Habier (2022).

The package provides the function \code{MLM(Y,X,Z)} for users only interested in solving multivariate mixed models, which takes as input the phenotypic matrix (\code{Y}), design matrices of fixed effects (\code{X}) and random effect (\code{Z}), and the same controls (\code{logtol,cores}) as the function \code{GWAS}. Unlike the GWAS solved, this function uses an efficient solver without single-value decomposition, as described by Xavier and Habier (2022).

}

\value{
  
Returns a list of class "GenModel" with two sublists: GBLUP and GWAS. The GBLUP list contains fixed and random effect coefficients, variance components, genetic correlations, heritability. The GWAS list provides the Wald statistics and QTL variances for each marker-trait combination, as well as the regression coefficients for each marker-population combination for the last trait.

}


\references{

Van Raden, P.M. and Jung, Y.C., 1988. A general purpose approximation to restricted maximum likelihood: the tilde-hat approach. Journal of Dairy Science, 71(1), pp.187-194.

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

# Compare before and after correction
plot( fit,  h2=TRUE, pch=20, col=8 ) # display QTL h2
lines( fit2$GWAS$QTLh2, pch=20, type='o' ) # adjusted QTL h2
legend('topright',pch=16,col=c(8,1),c('Before Beavis correction','After Beavis correction'))

}
