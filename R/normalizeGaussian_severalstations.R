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
#'prec_mes_gaussWilks <- WilksGaussianization(x=prec_mes, data=prec_mes,valmin=valmin,sample="monthly",extremes=TRUE,origin_x = origin, origin_data = origin) 
#'prec_mes_gauss <- normalizeGaussian_severalstations(x=prec_mes, data=prec_mes,step=0,sample="monthly",extremes=TRUE,origin_x = origin, origin_data = origin)  
#' 
#' CCGamma <- CCGamma(data=prec_mes, lag = 0,valmin=valmin,only.matrix=TRUE,tolerance=0.001)
#' prec_mes_ginv <- normalizeGaussian_severalstations(x=prec_mes_gaussWilks, data=prec_mes,step=0,sample="monthly",extremes=TRUE,origin_x = origin, origin_data = origin,inverse=TRUE)
#' 
#' 
#' str(prec_mes_ginv)
#' 
#'  plot(prec_mes[,1],prec_mes_ginv[,1])
#'  plot(prec_mes[,19],prec_mes_ginv[,19])
#'  plot(prec_mes[,5],prec_mes_ginv[,5])

#' plot(prec_mes_gaussWilks[,1],prec_mes_gauss[,1])
#' plot(prec_mes_gaussWilks[,2],prec_mes_gauss[,2])
#' plot(prec_mes_gaussWilks[,19],prec_mes_gauss[,19])
#' plot(prec_mes_gaussWilks[,10],prec_mes_gauss[,10])
#' plot(cor(prec_mes_gaussWilks),cor(prec_mes_gauss))
 
#
#' 
# ### ,step = valmin, prec = 10^-4, type = 3,
#		extremes = TRUE, sample = "monthly", origin_x = origin, origin_data = origin)
#
#
#
#' CCGamma <- CCGamma(data=prec_mes, lag = 0,valmin=valmin,only.matrix=TRUE,tolerance=0.001)
#
#
WilksGaussianization <- function (x,data=x,valmin=0.5,tolerance=0.001,prec_tolerance=0.05,iterations=15,force.precipitation.value=TRUE,seed=1234,...) {
	
	set.see(seed)
	lag=0
	x[x<valmin] <- 0 
	data[data<valmin] <- 0 
	gauss <- normalizeGaussian_severalstations(x=x,data=data,step=0,...)
	
	
	ccgamma <- CCGamma(data=x, lag = lag,valmin=valmin,only.matrix=TRUE,tolerance=tolerance)
	
	chol_ccg <- t(as.matrix(chol(ccgamma)))
	
	
	
	out <- gauss
	for (iter in 1:iterations) {
		print(iter)
		cor <- nearPD(cor(gauss),cor=TRUE)$mat
		chol_cor <- t(as.matrix(chol(cor))) 
		inv_chol_cor <- solve(chol_cor)
		
		for (r in 1:nrow(out)) {
			
			
			temp <- as.vector(gauss[r,])
			temp <- t(as.matrix(temp))
			
			#	str(temp)
			#	print(temp)
			#	str(inv_chol_cor)
			#	print(inv_chol_cor)
			xv <- inv_chol_cor %*% temp
			
			out[r,] <- chol_ccg %*% as.matrix(xv)
			
		}
	##	if (iterations>0) {
			
			message <- sprintf("Iteration: %d on %d",iter,iterations)
			
			x_rec <- normalizeGaussian_severalstations(x=out, data=data,step=0,inverse=TRUE,...)
			
			cond <- abs(x-x_rec)<prec_tolerance
			condhigh <- abs(x-x_rec)<prec_tolerance*50
			points <- length(cond)
			adjusted <- length(which(!cond))
			adjustedhigh <- length(which(!condhigh))
			message2 <- sprintf("Points adjusted: %d (%d) on %d",adjusted,adjustedhigh,points)
			print(message)  
			print(message2)
			gauss[cond] <- out[cond]
			
			
			
			
		
		
		
###		}
		
		
	}
	
	
	
	
	
	
	
	
	if (force.precipitation.value) {
		
	   cond <- x>valmin
	   out[cond] <- gass[cond]
		
	}
		
###	}  
	return(out)
	
	
}

