\name{generate.YuleWalkerCoefficientBlockmatrices}
\alias{generate}
\alias{generate.CCGammaObject}
\alias{generate.CCGammaObjectListPerEachMonth}
\alias{generate.YuleWalkerCoefficientBlockmatrices}
\alias{generate.YuleWalkerCoefficientBlockmatricesPerEachMonth}
\title{generate}
\usage{
\method{generate}{YuleWalkerCoefficientBlockmatrices}(x, FUN = rnorm,
  n = 100, names = NULL, xprev = NULL, names_x = c("A", "Sigma_u",
  "CCGammaInfo"), nearPD = TRUE, precipitation.indicator = FALSE, ...)

\method{generate}{YuleWalkerCoefficientBlockmatricesPerEachMonth}(x,
  FUN = rnorm, year_min = 1961, year_max = 1990, names = NULL,
  xprev = NULL, names_x = c("A", "Sigma_u", "CCGammaInfo"), nearPD = TRUE,
  precipitation.indicator = FALSE, ...)

\method{generate}{CCGammaObject}(x, n = 100, names = NULL, xprev = NULL,
  precipitation.indicator = TRUE, ...)

\method{generate}{CCGammaObjectListPerEachMonth}(x, year_min = 1961,
  year_max = 1990, names = NULL, xprev = NULL,
  precipitation.indicator = TRUE, ...)
}
\arguments{
  \item{x}{\code{YuleWalkerCoefficientBlockmatrices} S3
  object. See \code{\link{CoeffYWeq}}}

  \item{FUN}{random function of the probability
  distribution used for noise random generation. Default is
  \code{\link{rnorm}}. See
  \url{http://cran.r-project.org/web/views/Distributions.html}}

  \item{n}{number of generations requested}

  \item{names}{null object or string vectors or names of
  the variables to be generated simultaneously. Default is
  \code{NULL}.}

  \item{xprev}{null object or initial condition of the
  multivariate random process to be generated. Default is
  \code{NULL}.}

  \item{names_x}{names of the elements of a
  \code{YuleWalkerCoefficientBlockmatrices} S3 object. See
  examples.}

  \item{nearPD}{logical. If \code{TRUE} (Default) the
  function \code{\link{nearPD}} is applied to the
  covariance matrix of residuals. If \code{FALSE}, the
  estimated covariance matrix is not verified to be
  positive definite and error may occur during the function
  execution.}

  \item{precipitation.indicator}{logical value. Default is
  \code{FALSE}. If it is \code{TRUE}, the output is
  transformed to a number between 0 (no precipitation) and
  1 (high precipitation) taking into account the
  probabability of no occurence and the cumulate
  probability function referred to \code{FUN}.}

  \item{year_min,year_max}{first and last years of
  generation period}

  \item{...}{additional arguments for \code{FUN}}
}
\value{
a matrix or a data frame object
}
\description{
generate

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


## Not Run in the examples, uncomment to run the following lines
# origin <- paste(year_min,1,1,sep="-")


# coeff_monthly <- CoeffYWeq(data=prec_mes,p=1,tolerance=0.001,sample="monthly",origin=origin)


# generation_monthly <- generate(coeff_monthly,year_min=year_min,year_max=year_max,
#					names=names(prec_mes))


###  generation with CCGammaObject

# CCGamma <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE)

# generation_CCGamma <- generate(x=CCGamma,n=100,names=names(prec_mes))

# CCGamma_monthly <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE,
#                            sample="monthly",origin=origin)
## generation_CCGamma <- generate(x=CCGamma_monthly,year_min=year_min,year_max=year_max,
##                                names=names(prec_mes))
}
\seealso{
\code{\link{CoeffYWeq}},\code{\link{CCGammaToBlockmatrix}}
}

