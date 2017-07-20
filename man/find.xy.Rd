\name{find.xy}
\alias{find.xy}
\title{
Find variable pairs for a given set of LA scouting vectors. 
}
\description{
After finding the latent vectors, this function can be used to find the pairs of variables (rows of the original data matrix) that are associated with each latent vector in terms of dynamic correlation. 
}
\usage{
find.xy(array, z, fdr.cut=0.05, normalization="standardize",
 center.z=FALSE, lac.percentile=0.8)
}
\arguments{
  \item{array}{
The data matrix with variables in the rows and samples in the columns.
}
  \item{z}{
The matrix of latent variables. Each column is a latent vector.
}
  \item{fdr.cut}{
The threshold of local fdr for the selection of variable pairs. 
}
  \item{normalization}{
The way the data matrix is to be row-normalized. The method requires each row to have mean 0 and SD 1. There are two options, "standardize", or "normal score".
}
  \item{center.z}{
Whether to remove mean from each z vector.
}
  \item{lac.percentile}{
The variable pairs that are considered potentially dymanically correlated will enter the computation. This is determined by the percentile of the LAC score among all variable pairs. If the percentile is higher than the provided threshold, then they are considered. 
}
}

\value{
A list is returned. Each numbered item is a matrix with three columns: row number 1, row number 2, and local fdr value. Only those selected by fdr threhsold are returned.
}
\author{
Tianwei Yu <tianwei.yu@emory.edu>
}
\seealso{
dca()
}
\examples{
x<-la.simu.gen(n=100,p=80,n.grp=2, n.noise.gene=100, rho=0.5, pwr=0.25)
z<-dca(x$dat, n.fac=2)
xy<-find.xy(x$dat, z[[2]], fdr.cut=0.01)
summary(xy)
xy[[1]][1:5,]
xy[[2]][1:5,]
}