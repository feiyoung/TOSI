\name{linReg_TOSI}
\alias{linReg_TOSI}
\title{TOSI-assited regression.}
\description{
 Variable selection based on dichotomy.The TOSI also give the direction of choosing penalty parameter $lambda$, e.g., if
$H0S$ is rejected, we should choose smaller value for $lambda$ until $ HOSc$ is accepted $S$ is the estimates of $S0$ from penalized regression such as Lasso); If
$ tilde H0S$ is accepted, we should choose larger value for $lambda$ until $ H0Sc$ is rejected. If $H0Sc$  is accepted and $tilde H0S$ is rejected, then the  optimal  $lambda$  in some sense is obtained.
}
\usage{
  linReg_TOSI(X, Y, alpha=0.05, seed=1, sub.frac=0.3, standardized=F)
}
\arguments{
  \item{X}{a \code{n}-by-\code{p} matrix, the observed covariates matrix.}
  \item{Y}{a \code{n}-dimensional vector, the observed response vector.}
  \item{G1}{a positive vector with values between 1 and p, the set of H01.}
  \item{alpha}{a positive number, the significant level.}
  \item{sub.frac}{a positive number between 0 and 1, the proportion of the sample used in stage I.}
  \item{seed}{a non-negative integer, the random seed.}
  \item{standardized}{a logical value, whether standerdize in stage I.}
}

\value{
  return a list including the following components:
  \item{lmreg}{a list including the  results of linear regresion from
  \code{\link{lm}} based on selected variables.}
  \item{glmnet}{a list from \code{\link{glmnet}}, TOSI finely tune penalty parameters based on the results of glmnet.}
  \item{Infer}{a list including the TOSI inference results for the selected variables.}
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
  \code{\link{gendata_Reg}}, \code{\link{RegMax}}, \code{\link{RegMin}}
}
\examples{
  ### Example
  n <- 100; p <- 50
  s0 <-  3#
  rho <- 3# 3
  i <- 3
  dat1 <- gendata_Reg(n, p, s0, rho, seed=i)
  res <- linReg_TOSI(dat1$X, dat1$Y, alpha=0.05, seed=1, sub.frac=0.3, standardized=F)
}
\keyword{Simultanous inference}
