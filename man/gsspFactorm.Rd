\name{gsspFactorm}
\alias{gsspFactorm}
\title{High Dimensional Sparse  Factor Analysis Model}
\description{
  sparse factor analysis to extract latent linear factor and estimate row-sparse and entry-wise-sparse loading matrix.
}
\usage{
  gsspFactorm(X, q=NULL, lambda1=nrow(X)^(1/4), lambda2=nrow(X)^(1/4))
}
\arguments{
  \item{X}{a \code{n}-by-\code{p} matrix, the observed data}
  \item{q}{an integer between 1 and \code{p} or \code{NULL}, default as \code{NULL} and automatically choose q by the eigenvalue ratio method.}
  \item{lambda1}{a non-negative number, the row-sparse penalty parameter, default as \code{n}^(1/4).}
  \item{lambda2}{a non-negative number, the entry-sparse penalty parameter, default as \code{n}^(1/4).}
}

\value{
  return a list with class named \code{fac}, including following components:

  \item{hH}{a \code{n}-by-\code{q} matrix, the extracted lantent factor matrix.}
  \item{sphB}{a \code{p}-by-\code{q} matrix, the estimated row-sparseloading matrix.}
  \item{hB}{a \code{p}-by-\code{q} matrix, the estimated loading matrix without penalty.}
  \item{q}{an integer between 1 and \code{p}, the number of factor extracted.}
  \item{propvar}{a positive number between 0 and 1, the explained propotion of cummulative variance by the \code{q} factors.}
  \item{egvalues}{a n-dimensional(n<=p) or p-dimensional(p<n) vector, the eigenvalues of sample covariance matrix.}
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
  \code{\link{factor}}, \code{\link{Factorm}}
}
\examples{
  dat <- gendata_Fac(n = 300, p = 500)
  res <- gsspFactorm(dat$X)
  ccorFun(res$hH, dat$H0) # the smallest canonical correlation
  ## comparison of l2 norm
  par(mar=c(5,5,2,2), mfrow=c(1,2))
  plot(rowSums(dat$B0^2), type='o', ylab='l2B', main='True')
  l2B <- rowSums(res$sphB^2)
  plot(l2B, type='o', main='Est.')

  Bind <- ifelse(dat$B0==0, 0, 1)
  hBind <- ifelse(res$sphB==0, 0, 1)
  ## structure of loading matrix
  image(x=1:nrow(Bind), y= 1:ncol(Bind), z=Bind , col = 1:3,#hcl.colors(12, "YlOrRd",rev=F),
        ylab='columns', xlab='rows',cex.lab=1.8, cex.axis=1.6 )
  image(x=1:nrow(Bind), y= 1:ncol(Bind), z=hBind , col = 1:3,#hcl.colors(12, "YlOrRd",rev=F),
        ylab='columns', xlab='rows',cex.lab=1.8, cex.axis=1.6 )

  ## Select good penalty parameters
  dat <- gendata_Fac(n = 300, p = 200)
  res <- gsspFactorm(dat$X, lambda1=0.04*nrow(dat$X)^(1/4) ,lambda2=1*nrow(dat$X)^(1/4))
  ccorFun(res$hH, dat$H0) # the smallest canonical correlation

  ## comparison of l2 norm
  plot(rowSums(dat$B0^2), type='o', ylab='l2B', main='True')
  l2B <- rowSums(res$sphB^2)
  plot(l2B, type='o', main='Est.')

  ## comparison of structure of loading matrix
  Bind <- ifelse(dat$B0==0, 0, 1)
  hBind <- ifelse(res$sphB==0, 0, 1)
  image(x=1:nrow(Bind), y= 1:ncol(Bind), z=Bind , col = 1:3,#hcl.colors(12, "YlOrRd",rev=F),
        ylab='columns', xlab='rows',cex.lab=1.8, cex.axis=1.6, main='True')
  image(x=1:nrow(Bind), y= 1:ncol(Bind), z=hBind , col = 1:3,#hcl.colors(12, "YlOrRd",rev=F),
        ylab='columns', xlab='rows',cex.lab=1.8, cex.axis=1.6, main='Est.' )


}
\keyword{Factor}
\keyword{Feature}
