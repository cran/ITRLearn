\name{maximin}
\alias{maximin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Maixmin projection learning for optimal individualized treatment regime
}
\description{
  Derives a meaningful and reliable individualized treatment regime for future patients based on
  estimated groupwise contrast function.
}
\usage{
maximin(B, c0)	
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{An \eqn{p_1}*\eqn{G} matrix containing parameters in the groupwise contast function. Here \eqn{p_1} is the dimension 
  of \code{x.tau} and \eqn{G} is the number of subgroups. It does not contain the intercept term. It can be computed by \code{\link{MPL}}.
 }
  \item{c0}{The common marginal treatment effect shared by all subgroups. It can be computed by \code{\link{MPL}}.
  \code{\link{maximin}} to compute the maximin effects.
 }
}
\details{Denoted by \eqn{\beta_g} the \eqn{g}-th column of \code{B}. This function computes
\deqn{\arg\max_{\|(\beta^T,c)^T\|=1} \min_{g\in\{1,\dots,G\}} (\beta_g^T \beta+c_0 c).}
The above optimaization problem can be efficiently computed based on quadratic programming.
}
\value{
  A vector of maximin effects.
}
\references{
Shi, C., Song, R., Lu, W., and Fu, B. (2018). Maximin Projection Learning for Optimal Treatment 
Decision with Heterogeneous Individualized Treatment Effects. \emph{Journal of the Royal Statistical
Society, Series B,} \bold{ 80:} 681-702.
}
\author{
Chengchun Shi
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{MPL}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(12345)
X <- matrix(rnorm(1600), 800, 2)
A <- rbinom(800, 1, 0.5)
h <- 1+sin(0.5*pi*X[,1]+0.5*pi*X[,2])
tau <- rep(0, 800)
B <- matrix(0, 2, 4)
B[,1] <- c(2,0)
B[,2] <- 2*c(cos(15*pi/180), sin(15*pi/180))
B[,3] <- 2*c(cos(70*pi/180), sin(70*pi/180))
B[,4] <- c(0,2)
for (g in 1:4){
    tau[((g-1)*200+1):(g*200)] <- X[((g-1)*200+1):(g*200),]\%*\%B[,g]
}
## mean and scale of the subgroup covariates are allowed to be different
X[1:200,1] <- X[1:200,1]+1
X[201:400,2] <- 2*X[201:400,2]-1
X[601:800,] <- X[601:800,]/2
Y <- h+A*tau+0.5*rnorm(800)
G <- c(rep(1,200), rep(2,200), rep(3,200), rep(4,200))
result <- MPL(Y~X|A|G)
maximin(result$B, result$c0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Maximin projection learning}
