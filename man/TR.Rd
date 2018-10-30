\name{TR}
\alias{TR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Individualized treatment regime based on MPL.
}
\description{
Recommend individualized treatment regime for future patients, based on the 
maximin projection learning method.
}
\usage{
TR(object, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Fitted object of class "MPL".
}
  \item{x}{A matrix consisting of future patients baseline covariates. If there's only one group
  of patients in the observed data fitted by "MPL", then future patients are assuming to coming from
  this group and there is no restrictions on the number of observations in \code{x}. Otherwise,
  we require the number of observations in \code{x} to be larger than or equal to the dimension of \code{x}.
}
}
\value{
A vector of individualized treatments tailored for future patients.
}
\author{
Chengchun Shi
}
\seealso{\code{\link{MPL}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
ITR <- TR(result, matrix(rnorm(200), 100, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Individualized treatment regime}
