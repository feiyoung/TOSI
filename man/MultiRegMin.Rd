\name{MultiRegMin}
\alias{MultiRegMin}
\title{Multi-split-based Two-Stage Minimum Reg Test method for the Reg vector.}
\description{
  Conduct the simultaneous inference for a set of regression coefficients components as for the NUll hypothesises H02, where H02 assume the set of regression coefficients components exist zero.
}
\usage{
  MultiRegMin(X, Y,  G2, Nsplit = 5, alpha=0.05, seed=1, sub.frac=0.5, standardized=F)
}
\arguments{
  \item{X}{a \code{n}-by-\code{p} matrix, the observed covariates matrix.}
  \item{Y}{a \code{n}-dimensional vector, the observed response vector.}
  \item{G2}{a positive vector with values between 1 and p, the set of H02.}
  \item{Nsplit}{a positive integer, the random split times used, default as 5.}
  \item{alpha}{a positive number, the significant level.}
  \item{sub.frac}{a positive number between 0 and 1, the proportion of the sample used in stage I.}
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
  \code{\link{gendata_Reg}}
}
\examples{
  ### Example
  n <- 100; p <- 20;i <- 1
  s0 <- 5 # First five components are nonzeros
  rho <- 5;
  dat1 <- gendata_Reg(n, p, s0, seed=i, rho)
  # ex1: H01 is false
  MultiRegMin(dat1$X, dat1$Y, 1:s0)
  # ex1: H01 is true
  MultiRegMin(dat1$X, dat1$Y, p)
}
\keyword{Simultanous inference}
