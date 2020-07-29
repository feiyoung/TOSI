\name{MultiMeanMin}
\alias{MultiMeanMin}
\title{Multi-split-based Two-Stage Minimum Mean Test method for the mean vector.}
\description{
  Conduct the simultaneous inference for a set of mean components as for the Null hypothesises H02, where H02 assume the set of mean components exist zero.
}
\usage{
  MultiMeanMin(X, G2, Nsplit = 5, alpha=0.05, frac.size=0.5, seed=1, standardized=F)
}
\arguments{
  \item{X}{a \code{n}-by-\code{p} matrix, the observed data}
  \item{G2}{a positive vector with values between 1 and p, the set of H02.}
  \item{Nsplit}{a positive integer, the random split times used, default as 5.}
  \item{alpha}{a positive number, the significant level.}
  \item{frac.size}{a positive number between 0 and 1, the proportion of the sample used in stage I.}
  \item{seed}{a non-negative integer, the random seed.}
  \item{standardized}{a logical value, whether standerdize in stage I.}
}

\value{
  return a matrix with class "Min-test", row names 'chiq_test' , and column names 'CriticalValue', 'TestStatistic', 'reject_status', and 'p-value', including all the information about testing.
}
\references{Wei Liu, Huazhen Lin, Jin Liu (2020). Two-directional simultaneous inference for high-dimensional models.
}
\author{
  Liu Wei
}
\note{
  nothing
}
\seealso{
  \code{\link{gendata_Mean}}, \code{\link{MeanMin}}, \code{\link{MultiMeanMax}}
}
\examples{
  ### Example
  n <- 100; p <- 100; i <- 2
  s0 <- 5 # First five components are nonzeros
  rho <- 5; tau <- 1;
  dat1 <- gendata_Mean(n, p, s0, seed=i, rho, tau)
  # ex1: H01 is false
  MultiMeanMin(dat1$X, 1:s0)
  # ex1: H01 is true
  MultiMeanMin(dat1$X, 1:p)
}
\keyword{Simultanous inference}
