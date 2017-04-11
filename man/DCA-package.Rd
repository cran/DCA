\name{DCA-package}
\alias{DCA-package}
\alias{DCA}
\docType{package}
\title{
Dynamic Correlation Analysis for high dimensional data
}
\description{
Given a data matrix with variables in rows and samples in the columns, the method DCA finds dominant latent signals that regulate the dynamic correlation between many pairs of variables.
}
\details{
The subroutine dca() computes dynamic correlation signals from the data matrix. It can use PCA, SPCA, and kmeans clustering to find dominant signals. The subroutine find.xy() subsequently finds variable pairs that are associated with each latent signal. The subroutine plot.la() plots the dynamic correlation of two variables X and Y given the Z vector. 
}
\author{
Tianwei Yu <tianwei.yu@emory.edu>
}
\references{
Li, K.C. (2002) Genome-wide coexpression dynamics: theory and application, Proceedings of the National Academy of Sciences of the United States of America, 99, 16875-16880.
}
\keyword{ package }
