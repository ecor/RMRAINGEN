\name{generate.blockmatrix}
\alias{generate}
\alias{generate.blockmatrix}
\title{generate}
\usage{
\method{generate}{blockmatrix}(x, xprev = NULL, noise = noise, n = 10,
  x.noise.gen = NULL, ...)

\method{generate}{blockmatrix}(x, xprev = NULL, noise = noise, n = 10,
  x.noise.gen = NULL, ...)
}
\arguments{
  \item{x}{\code{blockmatrix} S3 object. See
  \code{\link{CoeffYWeq}}}

  \item{...}{further arguments}
}
\description{
generate
}
\details{
Implementation of \code{generate} method for
\code{YuleWalkerCoefficientBlockmatrices} or
\code{CCGammaObject} S3 object.  It generates a multivarite
random series according using a VAR model with coefficient
obtained by \code{\link{CoeffYWeq}}
(\code{YuleWalkerCoefficientBlockmatrices} S3 object) .
Alternatively it generates by applying a first-order Markov
Cain from \code{CCGammaObject} S3 Object.
}
\examples{
library(RMRAINGEN)

set.seed(125)
data(trentino)

year_min <- 1961
year_max <- 1990

period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
station <- names(PRECIPITATION)[!(names(PRECIPITATION) \%in\% c("day","month","year"))]
prec_mes <- PRECIPITATION[period,station]

## removing nonworking stations (e.g. time series with NA)
accepted <- array(TRUE,length(names(prec_mes)))
names(accepted) <- names(prec_mes)
for (it in names(prec_mes)) {
		 accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
}

prec_mes <- prec_mes[,accepted]
## the dateset is reduced!!!
prec_mes <- prec_mes[,1:2]

## Not Run in the examples, uncomment to run the following lines
# coeff <- CoeffYWeq(data=prec_mes,p=1,tolerance=0.001)

# generation <- generate(coeff,n=10,names=names(prec_mes))
}
\seealso{
\code{\link{generate}},\code{\link{CCGammaToBlockmatrix}},\code{\link{blockmatrix}}
}

