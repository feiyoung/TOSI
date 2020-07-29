\name{TSrowMinST}
\alias{TSrowMinST}
\title{Two-Stage Minimum Row Test method for the loading vectors.}
\description{
  Conduct the simultaneous inference for a set of loading vectors as for the NUll hypothesises H02, where H02 assume there is zero loading vector in the set of loading vectors.
}
\usage{
  TSrowMinST(X, G2,  alpha=0.05, seed=1, sub.frac=0.5)
}
\arguments{
  \item{X}{a \code{n}-by-\code{p} matrix, the observed data}
  \item{G2}{a positive vector with values between 1 and p, the set of H02.}
  \item{alpha}{a positive number, the significant level.}
  \item{seed}{a non-negative integer, the random seed.}
  \item{sub.frac}{a positive number between 0 and 1, the proportion of the sample used in stage I.}
}

\value{
   return a matrix with class "Min-test", row names 'chiq_test' , and column names 'CriticalValue', 'TestStatistic', 'reject_status', and 'p-value', including all the information about testing.
}
\references{Wei Liu, Huazhen Lin, Jin Liu (2020). Estimation and inference on high-dimensional sparse factor models.}
\author{
  Liu Wei
}
\note{
  nothing
}
\seealso{
  \code{\link{factor}}, \code{\link{Factorm}},\code{\link{simultestentryX}}
}
\examples{
  ### Example
  dat <- gendata_Fac(n = 300, p = 500)
  res <- Factorm(dat$X)
  X <- dat$X
  # ex1: H01 is false
  G2 <- 1:200; # all are nonzero loading vectors
  TSrowMinST(X, G2=G2, alpha=0.05, sub.frac=0.5)
  # ex2: H01 is true
  G2 <- 1:500 # all are zero loading vectors
  TSrowMinST(X, G2=G2, alpha=0.05, sub.frac=0.5)
}
\keyword{Simultanous inference}
\keyword{Feature}
