\name{plot_la}
\alias{plot_la}
\title{
Visualizing the dynamic correlation between a pair of genes given the AL scouting vector.
}
\description{
Given three vectors, x, y and z, the function produces a color scatter plot of x and y, colored by the grouping of z. 
}
\usage{
plot_la(x, y, z, use.locfdr = FALSE, cols = c("red", "green", "blue"), cex = 0.5)
}
\arguments{
  \item{x}{
The x vector.
}
  \item{y}{
The y vector.
}
  \item{z}{
The z vector, conditioned on which X and Y have dynamic correlation.
}
  \item{use.locfdr}{
Whether to use local fdr to group the z values. If TRUE, locfdr() is used and the fdr (posterior probability of belonging to the more extreme groups) threshold is 0.5. If FALSE, the z values are cut at the 0.33 and 0.67 quantiles into three groups.
}
  \item{cols}{
The colors of the groups. 
}
  \item{cex}{
The point size to be passed to plot().
}
}
\details{
The locfdr approach is only to be used when it is clear z has heavy tails, and it is believed such points are driving the dynamic correlation.
}
\value{
A vector, each item is a string, recording the LA score and the correlation in each of the point groups as determined by z values.
}
\author{
Tianwei Yu <tianwei.yu@emory.edu>
}
\examples{
x<-la.simu.gen(n=100,p=20,n.grp=2, n.noise.gene=10, rho=0.25, pwr=1)
plot_la(x$dat[1,], x$dat[2,], x$z[,1],use.locfdr=FALSE)
}