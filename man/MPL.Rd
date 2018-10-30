\name{MPL}
\alias{MPL}
\alias{MPL.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Maixmin projection learning for optimal individualized treatment regime
}
\description{
  Derives a meaningful and reliable individualized treatment regime based on the observed
  dataset from different subgroups with heterogeneity in optimal individualized treatment
  decision making. When patients are coming from the same group, it implements the classical
  Q learning and A learning algorithm.
}
\usage{
MPL(formula, data, subset, na.action, method = c("Q", "A"), bootstrap = FALSE,
	 control = MPL.control(...), model = TRUE, y = TRUE, a = TRUE, 
	 g = TRUE, x.tau = TRUE, x.h = TRUE, x.pi = TRUE, random = FALSE, ...)
	
MPL.fit(y, x.tau, a, g=NULL, x.h=NULL, x.pi=NULL, method=c("Q", "A"),
			bootstrap=FALSE, random=FALSE, control=MPL.control())	
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A symbolic description of the model to be fitted(of type y ~ x.tau | a, or
  y ~ x.tau | a | g, or y ~ x.tau | a | g | x.h, or y ~ x.tau | a | g | x.h | x.pi, or
  y ~ x.tau | a | g | | x.pi. Details 
  are given in 'Details').
}
  \item{data}{An optional list or environment containing variables in \code{formula}.
}
  \item{subset, na.action}{Arguments controlling formula processing via \code{\link{model.frame}}.
}
  \item{method}{Method used for estimating the parameter in the groupwise contrast function. See 'Details'.
}
  \item{bootstrap}{A logical value indicating whether bootstrap will be used. Default is \code{FALSE}. 
  See 'Details'.
}
  \item{control}{A list of control argument via \code{\link{MPL.control}}.
}
  \item{model}{A logical value indicating whether \emph{model frame} should be included
  as a component of the return value.
}
  \item{y, a, g, x.tau, x.h, x.pi}{For \code{MPL}: logical values indicating whether the response, 
  the treatment, the subgroup indicator, covariates used to fit the contrast function, covariates
  used to fit the baseline function and covariates used to fit the propensity score function.
  
  For \code{MPL.fit}: y is the response vector (the larger the better), a is the treatment vector  
  denoting the treatment patients receive, g is the group indicator indicating which group each patient 
  belongs to, x.tau, x.h, x.pi are the design matrices used to fit the contrast, the baseline and
  the propensity score function.
}
  \item{random}{A logical value indicating whether using a constant to fit the propensity score function
  or not. In randomized studies, the propensity score is usually a constant function independent of baseline
  covariates. When \code{random=TRUE}, \code{MPL} uses a constant to fit the propensity score. Otherwise,
  it uses a logistic regression function based on covariates in \code{x.pi}.
}
  \item{\dots}{Argument passed to \code{\link{MPL.control}}.
}
}
\details{A saline feature of data from clinical trials and medical studies is inhomogeneity. Patients not 
only differ in baseline characteristics, but also the way they respond to the treatment. Individualized
treatment regimes are developed to select effective treatments based on patient's heterogeneity. Formally
speaking, an individualized treatment regime (ITR) is a function that maps patients' baseline covariates to 
the space of available treatment options. The goal in precision medicine is to identify the optimal ITR to
reach the best clinical outcomes. 

However, the optimal ITR might also vary for patients across different subgroups. This function implements 
the maximin projection learning method that derives a meanful and reliable ITR for future patients based on 
the observed data from different populations with heterogeneity in optimal individualized decision making.

The means and covariance matrices of patients baseline covariates are allowed to vary across different subgroups.
\code{MPL} will first standardize the groupwise baseline covariates to have zero mean and indentity covariance
matrix (based on Gram-Schmidt Orthonormalization) and then recommends an ITR for future groups of patients. Notice
that the resulting ITR cannot be directly applied to future patients. We need to standardize future patients 
baseline covariates (based on the same procedure) first and apply the transformed covariates to the ITR. This is
implemented by the \code{TR} function.

When the group indicator \code{g} is omitted (or it is a constant vector) in the formula, \code{MPL} assumes  
all the patients are coming from the same group and implements the classical Q-learning and A-learning algorithm.
Otherwise, \code{g} should be a numeric vector that has the same length of \code{y}, indicating which group
each patient belongs to. 

When \code{x.h} is omitted and the baseline \code{h.est} in \code{MPL.control} is not specified, \code{MPL} 
sets \code{x.h=x.tau}. When \code{x.pi} is omitted, the propensity score \code{pi.est} in \code{MPL.control} is
not specified, and \code{random=FALSE}, \code{MPL} sets \code{x.pi=x.tau}. 

Q-learning fits the entire Q function (the conditional mean of response given baseline covariates and treatment)
to derive the optimal ITR. A-learning is a more robust method that focuses directly on the contrast function (the 
difference between two Q functions). It requires to specify both the baseline and the propensity score function 
and the resulting estimator for the contrast function is consistent when either of the function is correctly 
specified. This is referred to as the doubly robustness property of A-learning. \code{MPL} uses Q-learning or 
A-learning to estimate the groupwise contrast function that shares the same marginal treatment effects across 
different subgroups. These estimators are further used to derive a ITR for future groups of patients. By default, 
\code{method="A"} and A-learning is implemented.

Inference for the maximin effects and the parameters in the groupwise contrast functions are conducted based on 
bootstrap. By default, \code{bootstrap=FALSE} and Bootstrap will not be conducted. 
}
\value{
  \item{Theta.tau.est}{An \eqn{(p_1+1)\times G} matrix containing estimated parameters in the groupwise contrast function. Here
  \eqn{p_1} is the dimension of \code{x.tau} and \eqn{G} is the number of subgroups. The first row contains the intercept term.}
  \item{Theta.h.est}{An \eqn{(p_2+1)\times G} matrix containing estimated parameters in the groupwise baseline function. Here
  \eqn{p_2} is the dimension of \code{x.h} and \eqn{G} is the number of subgroups. The first row contains the intercept term. It
  equals \code{NULL} when \code{h.est} in \code{\link{MPL.control}} is prespecified.}
  \item{Theta.pi.est}{An \eqn{(p_3+1)\times G} matrix containing estimated parameters in the groupwise propensity score function. 
  Here \eqn{p_3} is the dimension of \code{x.pi} and \eqn{G} is the number of subgroups. The first row contains the intercept term. 
  It equals \code{NULL} when \code{pi.est} in \code{\link{MPL.control}} is prespecified.}
  \item{h.est}{Estimated baseline function.}
  \item{pi.est}{Estimated propensity score function.}
  \item{B}{An \eqn{p_1\times G} matrix containing estimated parameters in the groupwise contast function. Here \eqn{p_1} is the dimension 
  of \code{x.tau} and \eqn{G} is the number of subgroups. It does not contain the intercept term. These parameters are the 
  corresponding coefficients of the transformed covariates and are thus different from \code{Theta.tau.est}. It can be
  used as the input of the function \code{\link{maximin}} to compute the maximin effects.}
  \item{c0}{The common marginal treatment effect shared by all subgroups. It can be used as the input of the function 
  \code{\link{maximin}} to compute the maximin effects.}
  \item{beta.est}{The estimated maximin effects used to contruct ITR for future patients.}
  \item{Theta.tau.boot}{An \eqn{(p_1+1)\times G\times B_0} array containing bootstrap samples for the estimated parameters in the groupwise 
  contrast function. Here \eqn{p_1} is the dimension of \code{x.tau}, \eqn{G} is the number of subgroups and \eqn{B_0} is the number of
  bootstrap samples. It equals \code{NULL} when \code{bootstrap=FALSE}.}
  \item{Theta.h.boot}{An \eqn{(p_2+1)\times G\times B_0} array containing bootstrap samples for the estimated parameters in the groupwise 
  baseline function. Here \eqn{p_2} is the dimension of \code{x.h}, \eqn{G} is the number of subgroups and \eqn{B_0} is the number of
  bootstrap samples. It equals \code{NULL} when \code{bootstrap=FALSE} or \code{h.est} in \code{\link{MPL.control}} is 
  prespecified.}
  \item{Theta.pi.boot}{An \eqn{(p_3+1)\times G\times B_0} array containing bootstrap samples for the estimated parameters in the groupwise 
  propensity score function. Here \eqn{p_3} is the dimension of \code{x.pi}, \eqn{G} is the number of subgroups and \eqn{B_0} is the number of
  bootstrap samples. It equals \code{NULL} when \code{bootstrap=FALSE} or \code{pi.est} in \code{\link{MPL.control}} is 
  prespecified.}
  \item{beta.boot}{An \eqn{p_1\times B_0} matrix containing bootstrap sample for the estimated maximin effects. Here \eqn{p_1} is the dimension 
  of \code{x.tau} and \eqn{B_0} is the number of bootstrap samples. It equals \code{NULL} when \code{bootstrap=FALSE}
  }
  \item{standardize}{A logical value indcating whether future patients covariates should be standardized first to be applied 
  to the ITR constructed by the maximin effects. \code{TRUE} if there are multiple subgroups and \code{FALSE} otherwise.}
  \item{model}{The full model frame (if \code{model = TRUE}).}
  \item{y}{Response vector (if \code{y = TRUE}).}
  \item{x.tau}{Covariates used to model the contrast function (if \code{x.tau = TRUE}).}
  \item{a}{Treatment vector (if \code{a = TRUE}).}
  \item{g}{Group Indicator (if \code{g = TRUE}).}
  \item{x.h}{Covariates used to model the baseline function (if \code{x.h = TRUE}).}
  \item{x.pi}{Covariates used to model the propensity score function (if \code{x.pi = TRUE}).}
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

\seealso{\code{\link{MPL.control}}, \code{\link{TR}}, \code{\link{maximin}}
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

## Q-learning
result <- MPL(Y~X|A|G, method="Q")

## A-learning
result <- MPL(Y~X|A|G)

\donttest{
result <- MPL(Y~X|A|G, bootstrap=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Maximin projection learning}
