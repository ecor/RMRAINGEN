\name{normalizeGaussian_severalstationsWilks}
\alias{normalizeGaussian_severalstationsWilks}
\title{Wilsks Correctio for Precipitation Marginal Gausianization}
\usage{
normalizeGaussian_severalstationsWilks(x, data = x, valmin = 0.5,
  tolerance = 0.001, ...)
}
\arguments{
  \item{x}{see
  \code{\link{normalizeGaussian_severalstations}}}

  \item{valmin,tolerance}{see \code{\link{CCGamma}}}

  \item{...}{further arguments for
  \code{\link{normalizeGaussian_severalstations}}}
}
\description{
Wrapper Gaussianization of
\code{\link{normalizeGaussian_severalstations}}
}
\examples{
library(RMRAINGEN)


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
valmin <- 0.5
prec_mes_gaussWilks <- normalizeGaussian_severalstationsWilks(x=prec_mes, data=prec_mes,valmin=valmin,sample="monthly",extremes=TRUE,origin_x = origin, origin_data = origin)
prec_mes_gauss <- normalizeGaussian_severalstations(x=prec_mes, data=prec_mes,step=0,sample="monthly",extremes=TRUE,origin_x = origin, origin_data = origin)

CCGamma <- CCGamma(data=prec_mes, lag = 0,valmin=valmin,only.matrix=TRUE,tolerance=0.001)
plot(prec_mes_gaussWilks[,1],prec_mes_gauss[,1])
plot(prec_mes_gaussWilks[,2],prec_mes_gauss[,2])
plot(prec_mes_gaussWilks[,19],prec_mes_gauss[,19])
plot(prec_mes_gaussWilks[,10],prec_mes_gauss[,10])
plot(cor(prec_mes_gaussWilks),cor(prec_mes_gauss))
prec_mes_ginv <- normalizeGaussian_severalstations(x=prec_mes_gaussWilks, data=prec_mes,step=0,sample="monthly",extremes=TRUE,origin_x = origin, origin_data = origin,inverse=TRUE)
 str(prec_mes_ginv)

plot(prec_mes[,1],prec_mes_ginv[,1])
 plot(prec_mes[,19],prec_mes_ginv[,19])
 plot(prec_mes[,5],prec_mes_ginv[,5])

CCGamma <- CCGamma(data=prec_mes, lag = 0,valmin=valmin,only.matrix=TRUE,tolerance=0.001)
}
\seealso{
\code{\link{normalizeGaussian_severalstations}},\code{\link{CCGamma}}
}

