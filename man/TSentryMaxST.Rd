\name{TSentryMaxST}
\alias{TSentryMaxST}
\title{Two-Stage Maximum Entry Test method for the loading vectors.}
\description{
  Conduct the simultaneous inference for a set of loading vectors as for the NUll hypothesises H01, where H01 assume the set of loading vectors are all zeroes.
}
\usage{
  TSentryMaxST(X, S1,  alpha=0.05, seed=1, sub.frac=0.5)
}
\arguments{
  \item{X}{a \code{n}-by-\code{p} matrix, the observed data}
  \item{S1}{a matrix with two column whose rows represent bi-index, where the first index is between 1 and p, and the second index is between 1 and q, the set of H01.}
  \item{alpha}{a positive number, the significant level.}
  \item{seed}{a non-negative integer, the random seed.}
  \item{sub.frac}{a positive number between 0 and 1, the proportion of the sample used in stage I.}
}

\value{
  return a matrix with class "Max-test", row names 'chiq_test' , and column names 'CriticalValue', 'TestStatistic', 'reject_status', and 'p-value', including all the information about testing.
}
\references{Wei Liu, Huazhen Lin, Jin Liu (2020). Estimation and inference on high-dimensional sparse factor models.
}
\author{
  Liu Wei
}
\note{
  nothing
}
\seealso{
  \code{\link{factor}}, \code{\link{Factorm}},\code{\link{TSrowMaxST}}
}
\examples{
  ### Example
  p <- 500; q <- 2
  dat <- gendata_Fac(n = 300, p = p, q=q)
  nzindMat <- indvec2matFun(which(dat$B0!=0), nrow=p)
  zindMat <- indvec2matFun(which(dat$B0==0), p)
  # ex1: H01 is false
  S1 <- rbind(c(1,2), c(1,3)); # all are nonzero loading entries
  TSentryMaxST(dat$X, S1=S1, alpha=0.05, sub.frac=0.4)
  # ex2: H01 is true: all are zero entries
  TSentryMaxST(dat$X, S1=zindMat, alpha=0.05)
  # ex3: H01 is false: no zero entry.
  TSentryMaxST(dat$X, S1=nzindMat, alpha=0.05)
}
\keyword{Simultanous inference}
\keyword{Feature}
