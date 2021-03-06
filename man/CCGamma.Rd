\name{CCGamma}
\alias{CCGamma}
\title{This function extends \code{\link{continuity_ratio}} and adds the corresponding gaussian correlation matrix for no-precipitation occurence.}
\usage{
CCGamma(data, lag = 0, p0_v1 = NULL, p = NA, valmin = 0.5,
  nearPD = (lag >= 0), interval = c(-1, 1),
  tolerance = .Machine$double.eps, only.matrix = FALSE,
  return.value = NULL, null.gcorrelation = 1e-05, sample = NULL,
  origin = "1961-1-1", ...)
}
\arguments{
  \item{data}{data frame or 'zoo' R object containing daily
  precipitation time series for several gauges (one gauge
  time series per column). See
  \code{\link{continuity_ratio}}.}

  \item{lag}{numeric lag (expressed as number of days) used
  for computation for "cross" continuity ratio and joint
  probability of prercipitation (no)occurence. See
  \code{\link{continuity_ratio}}.}

  \item{p}{positive integer parameter. Default is
  \code{NA}, otherwise, \code{lag} is calculated as the
  vector \code{0:p}.}

  \item{valmin}{threshold precipitation value [mm] for
  wet/dry day indicator. If precipitation is lower than
  \code{valmin}, day is considered dry. Default is 0.5 mm.
  See \code{\link{continuity_ratio}}.}

  \item{p0_v1}{vector for marginal probablities, see
  \code{\link{omega}} and \code{\link{omega_inv}}.}

  \item{nearPD}{see \code{\link{omega_inv}}. Default is
  \code{(lag==0)}.}

  \item{interval,tolerance}{see \code{\link{omega_inv}}}

  \item{only.matrix}{logical value. If \code{TRUE} the
  function returns only the gaussian correlaton matrix.
  Deafaul is \code{FALSE}.}

  \item{return.value}{string. If it is not either
  \code{NULL} (Default) and \code{NA}, function returns
  only the argument indicated by this argument.}

  \item{null.gcorrelation}{numerical value
  \code{nooccurence_gcorrelation} under which is considered
  to be 0.}

  \item{sample}{character string indicated if function must
  be calculated differently for subset of the year, e.g.
  monthly. Admitted values are \code{NULL} (Default),
  \code{"all"} or \code{"monthly"}.}

  \item{origin}{character string (yyyy-dd-mm) indicated the
  date of the first row of \code{"data"}. It is used if
  \code{data} and \code{sample} are not \code{NULL}.}

  \item{...}{additional agruments of
  \code{\link{omega_inv}} or \code{\link{CCGamma}}}
}
\value{
An object which is a list containing the following fields:

\code{continuity_ratio} : \code{lag}-day lagged continuity
ratio, as returned by \code{\link{continuity_ratio}};

\code{occurence} : joint probability of \code{lag}-day
lagged precipitation occurence, as returned by
\code{\link{continuity_ratio}};

\code{nooccurence} : joint probability of \code{lag}-day
lagged no precipitation occurence, as returned by
\code{\link{continuity_ratio}};

\code{lag} : number of days lagged between the two compared
events (see argument \code{lag});

\code{p0_v1} : vector of marginal probability of no
precipitation occurence. If \code{lag} is 0, it corresponds
to the diagonal of \code{nooccurence} matrix (see argument
\code{p0_v1});

\code{nooccurence_gcorrelation} corresponding gaussian
correlation for no precipitation occurence obtained by
applying \code{\link{omega_inv}} to \code{nooccurence},

\code{TransintionMatrixMCFirstOrder} List of transition
matrices (with condinitioned probability values) of the
first-order Markov Chain.  See
\code{\link{TransitionMatrixMCFirstOrder}}.

If the argument \code{only.matrix} is \code{TRUE}, only
\code{nooccurence_gcorrelation} is returned as a matrix. In
case the argument \code{lag} is a vector wirh length more
than one, the function returns a list of the above-cited
return object for each value of the vector \code{lag}.
}
\description{
This function extends \code{\link{continuity_ratio}} and
adds the corresponding gaussian correlation matrix for
no-precipitation occurence.
}
\note{
This functon is useful to generate the serial
cross-correlation matrices for no precipitation occurrence
for Yule-Walker Equations. In case \code{lag} is a vactor,
\code{nearPD} must be a vector of the same size, default is
\code{(lag==0)}.

See the R code for major details
}
\examples{
library(RMAWGEN)
data(trentino)

year_min <- 1961
year_max <- 1990
origin <- paste(year_min,1,1,sep="-")

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

CCGamma <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE)
## Not Run in the examples, uncomment to run the following line
## CCGamma <- CCGamma(data=prec_mes,lag=0:2,tolerance=0.001,only.matrix=FALSE)

## Not Run in the examples, uncomment to run the following line
## CCGamma_monthly <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE,
#                     sample="monthly",origin=origin)
}
\author{
Emanuele Cordano
}
\seealso{
\code{\link{continuity_ratio}},\code{\link{omega_inv}},\code{\link{omega}},\code{\link{CCGammaToBlockmatrix}}
}

