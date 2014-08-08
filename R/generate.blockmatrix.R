NULL
#' generate
#' 
#' Implementation of \code{generate} method for \code{YuleWalkerCoefficientBlockmatrices} or \code{CCGammaObject} S3 object.  It generates a multivarite random series according using a VAR model with coefficient obtained by \code{\link{CoeffYWeq}} (\code{YuleWalkerCoefficientBlockmatrices} S3 object) . Alternatively it generates by applying a first-order Markov Cain from \code{CCGammaObject} S3 Object.
#'
#' @param x \code{blockmatrix} S3 object. See \code{\link{CoeffYWeq}}
#' @param ... further arguments 
#' 
#' 
#' @title generate
# @name generate
# @rdname generate
#' @rdname generate.blockmatrix
#' @method generate blockmatrix
#' @S3method generate blockmatrix
#' @aliases generate generate.blockmatrix 
#' @importFrom RGENERATE generate
#' @export
# @import methods
#' @seealso \code{\link{generate}},\code{\link{CCGammaToBlockmatrix}},\code{\link{blockmatrix}}
#' 
#' @examples 
#' 
#' library(RMRAINGEN)
#' 
#' set.seed(125)
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
#' 
#' period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
#' station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
#' prec_mes <- PRECIPITATION[period,station]  
#' 
#' ## removing nonworking stations (e.g. time series with NA)
#' accepted <- array(TRUE,length(names(prec_mes)))
#' names(accepted) <- names(prec_mes)
#' for (it in names(prec_mes)) {
#' 		 accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
#' }
#'
#' prec_mes <- prec_mes[,accepted]
#' ## the dateset is reduced!!! 
#' prec_mes <- prec_mes[,1:2]
#' 
#' ## Not Run in the examples, uncomment to run the following lines 
#'  coeff <- CoeffYWeq(data=prec_mes,p=2,tolerance=0.001)
#' 
#' generation <- generate(coeff$A)
#' 
#'



generate.blockmatrix <- function(x,xprev=NULL,noise=NULL,n=10,x.noise.gen=NULL,is.VAR=TRUE,...) {
	
	
	out <- NULL
	p <- nrow(x)
	K <- nrow(x[1,1])
	print(p)
	
	if (is.null(xprev)) xprev <- generate(x.noise.gen,n=p,K=K,...)
	if (!is.null(xprev)) { 
		
		
		if (nrow(xprev)>1) {
			xprevv <- as.vector(xprev[1,])
			names(xprevv) <- names(xprev)
			for (r in 2:nrow(xprev)) {
				temp <- as.vector(xprev[r,])
				print(temp)
				print(xprevv)
				names(temp) <- paste(names(temp),(r-1),sep=".l")
				xprevv <- cbind(xprevv,temp)
				
				
				
			}
			xprev <- as.vector(xprevv)
			xprev <- t(xprev)
		}	else if (nrow(prev)==1) {
			
			xprev <- t(xprev)
		}
		
	
	} 
	
	if (length(xprev)!=p*K) { stop("Error in generate.blockamatix: xprev not consistent!!!")}
	
	
	
	if (p>1) count <- 1:p 
	
	
	if (!is.null(noise)) n <- nrow(noise)
	
	if (is.null(noise))  noise <- generate(x.noise.gen,n=n,K=K,...)
		
	if (p>1) {
		
		temp <- noise	
		for (l in 1:(p-1)) {
				
			###count <- (p+1):n-l
			temp  <- noise
			if (l>0) names(temp) <- paste(names(temp),l,sep=".l")
			if (l>0) temp <- temp*0
			noisep <- cbind(noisep,temp)
				
		}
			
			####print(noisep)
			###### TO DOOO
			### GENERARE DA P+1 IN POI ......
			
	}
	
	out <- generate(as.matrix(x),noise=noisep,names=names(noisep),xprev=xprev)
	
	
	###out <- rbind(out[1:p,]*NA,out)
	rownames(out) <- rownames(noise)
	if (is.VAR) out <- out[,1:K]
		
	return(out)
} 

