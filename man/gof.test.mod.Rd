\name{gof.test.mod}
\alias{gof.test.mod}
\title{Modified Signature for \code{\link{ks.test}} or other Goodness-of-fit tests for two samples}
\usage{
gof.test.mod(x, y = "pexp", what = "ks.test", par = list(), ...)
}
\arguments{
  \item{x}{sample}

  \item{y}{a vector (second sample) or a string name of
  cumulative probabiity function, e. g. \code{"exp"} for
  Expontantial Distribution (\code{\link{pexp}}).}

  \item{what}{name of the test function, e. g.
  \code{\link{ks.test}}.}

  \item{par}{list of distribution parameters}

  \item{...}{further argumets}
}
\value{
The value of the function given through \code{what}
argument.
}
\description{
Modified Signature for \code{\link{ks.test}} or other
Goodness-of-fit tests for two samples
}
\seealso{
\code{\link{ks.test}},\code{\link{do.call}}
}

