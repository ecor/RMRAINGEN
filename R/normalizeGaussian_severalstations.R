# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#' Wilsks Correctio for Precipitation Marginal Gausianization
#' 
#' Wrapper Gaussianization of \code{\link{normalizeGaussian_severalstations}} 
#' 
#' 
#' 
#' @param x see  \code{\link{normalizeGaussian_severalstations}}
#' @param valmin,tolerance see \code{\link{CCGamma}}
#' @param ... further arguments for \code{\link{normalizeGaussian_severalstations}}
#' 
#' 
#' @export 
#' 
#' 
#' 
#' 
#' 
#' 
#' @seealso \code{\link{normalizeGaussian_severalstations}},\code{\link{CCGamma}}
#' 
#' @examples
#' 
#'
#'library(RMRAINGEN)
#'
#'
#'data(trentino)
#'
#'year_min <- 1961
#'year_max <- 1990
#'origin <- paste(year_min,1,1,sep="-")
#'
#'period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
#'station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
#'prec_mes <- PRECIPITATION[period,station]
#'
#'
#'
#' ## removing nonworking stations (e.g. time series with NA)
#'accepted <- array(TRUE,length(names(prec_mes)))
#' names(accepted) <- names(prec_mes)
#'for (it in names(prec_mes)) {
#'	accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
#'}
#'
#'prec_mes <- prec_mes[,accepted]
#'valmin <- 0.5
#'prec_mes_gauss <- normalizeGaussian_severalstations(x=prec_mes, data=prec_mes,valmin=valmin) 
### ,step = valmin, prec = 10^-4, type = 3,
#		extremes = TRUE, sample = "monthly", origin_x = origin, origin_data = origin)
#
#
#
#CCGamma <- CCGamma(data=prec_mes, lag = 0,valmin=valmin,only.matrix=TRUE,tolerance=0.001)
#
#
normalizeGaussian_severalstationsWilks <- function (x,valmin=0.5,tolerance=0.001,...) {
	
	lag=0
	out <- normalizeGaussian_severalstations(x,step=valmin,...)
	cor <- nearPD(cor(prec_mes_gauss),cor=TRUE)$mat
	
	ccgamma <- CCGamma(data=prec_mes, lag = lag,valmin=valmin,only.matrix=TRUE,tolerance=tolerance)
	
	chol_cor <- t(as.matrix(chol(cor))) 
	inv_chol_cor <- solve(chol_cor)
	chol_ccg <- t(as.matrix(chol(ccgamma)))
	
	for (r in 1:nrow(out)) {
		
		xv <- inv_chol_cor %*% as.vector(out[,r])
		
		out[r,] <- chol_ccg %*% as.matrix(xv)
		
	}
	
	return(out)
	
	
}

