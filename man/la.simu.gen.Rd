\name{la.simu.gen}
\alias{la.simu.gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a data matrix with underlying dynamic correlation signal in the Liquid Association (LA) framework
}
\description{
The simulation follows the LA framework, namely dynamic correlation is in the form of X~N(0,1), Y~N(0,1), Z~N(0,1), E(XY|Z) is a function of Z.
}
\usage{
la.simu.gen(n, p, n.grp, n.noise.gene, rho, pwr)
}
\arguments{
  \item{n}{
Sample size (number of columns of the data matrix).
}
  \item{p}{
The number of genes in each LA module, i.e. a group of genes regulated by the same latent dynamic correlation factor.
}
  \item{n.grp}{
The number of LA modules to simulate.
}
  \item{n.noise.gene}{
The number of pure noise genes to add to the matrix.
}
  \item{rho}{
The standard deviation of the Gaussian noise to be added to the simulated data in the modules.
}
  \item{pwr}{
The power for the transformation (see details)
}
}
\details{
Between modules, the latent LA factor z's are independent.

Within each module, 10 sub-modules are simulated. For each sub-module, we first generate a pair of X and Y vectors, which follows:

X~N(0,1), Y~N(0,1)
u=(pnorm(z)-0.5)*2
E(XY|z)=sign(u)*abs(u)^pwr

Then white noise with SD of rho is added to the hidden X, Y pair to generate pairs of observed X, Y vectors. 
}
\value{
A list is returned.
  \item{dat }{The data matrix.}
  \item{z }{The true z vectors.}
}
\author{
Tianwei Yu <tianwei.yu@emory.edu>
}
\seealso{
dca()
}
\examples{
x<-la.simu.gen(n=100,p=200,n.grp=3, n.noise.gene=100, rho=0.5, pwr=1)
x$dat[1:5,1:5]
x$z[1:5,]
}
