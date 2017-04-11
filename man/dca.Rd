\name{dca}
\alias{dca}
\title{
Dynamic Correlation Analysis
}
\description{
The method finds a series of latent vectors, which serve as the LA scouting vectors for large numbers of variable pairs. 
}
\usage{
dca(array, top.pairs.prop = 0.95, max.pairs = 1e+06, n.fac = 10, 
sumabsv = sqrt(max.pairs)/10, normalization = "standardize", method = "PCA")
}
\arguments{
  \item{array}{
The data matrix, with variables in the rows and samples in the columns. 
}
  \item{top.pairs.prop}{
The method ranks all variable pairs from the most likely to have dynamic correlation relationship to the least likely. The top pairs are used for detection of latent signals. This parameter controls the percentage of pairs used in the computation.
}
  \item{max.pairs}{
The maximumn number of pairs to use. When the data contains too many variables, such as tens of thousands of variables in a gene expression matrix, this parameter limits the maximumn number of variable pairs to enter the calculation.
}
  \item{n.fac}{
The number of top latent factors to report. If the method "kmeans" is used, this parameter is used as the number of clusters. 
}
  \item{sumabsv}{
The sumabsv parameter to be passed on to the SPC() method.
}
  \item{normalization}{
The way the data matrix is to be row-normalized. The method requires each row to have mean 0 and SD 1. There are two options, "standardize", or "normal score".
}
  \item{method}{
The method for finding the latent factors. Current choices are "PCA", "SPCA", and "kmeans".}
}
\details{
	After finding the factors, the method attemps to rotate the factor using oblique rotation to achieve more interpretable results. 
}
\value{
The method returns a list.
 \item{fac}{The original factors found. This is the PC, SPC, or cluster mean vector depending on the method chosen.}
 \item{rotated}{The factors after rotation.}
 \item{ss.proj}{The sum of squared attributed to each rotated factor.}

}
\author{
Tianwei Yu <tianwei.yu@emory.edu>}
\seealso{
find.xy()
}
\examples{
x<-la.simu.gen(n=100,p=200,n.grp=2, n.noise.gene=100, rho=0.5, pwr=0.5)
z<-dca(x$dat, n.fac=2)
cor(z[[2]], x$z, method="spearman")
}