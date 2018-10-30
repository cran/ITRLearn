\name{MPL.control}
\alias{MPL.control}
\title{Control parameters for maximin projection learning}
\description{
  Parameters that control fitting of maximin projection learning.
}
\usage{
MPL.control(pi.est = NULL, h.est = NULL, boot.sample = 600)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pi.est}{Estimated propentisy score for each patient. If not specified, a 
  logistic regression model is fitted to estimate the propensity score.
}
  \item{h.est}{Estimated baseline function for each patient. If not specified, a linear
  regression model is fitted to estimate the baseline function.
}
  \item{boot.sample}{Number of bootstrap samples used for inference of the maximin effects
  and the subgroup parameter. Default is 600.  
}
}
\value{
  A list with the arguments specified. 
}

\seealso{\code{\link{MPL}}, \code{\link{MPL.fit}}
}

