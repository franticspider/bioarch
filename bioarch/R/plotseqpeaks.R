




#' plot the peaks of a mass spec
#' 
#' @param data a data strucyture in MALDI format
#' @param myxlim c(min,max) for the range of the data you want to plot
#' @keywords bruker
#' @export
#' @examples
#' plotseqpeaks(cd1,myxlim)
ba_plotseqpeaks <- function(data,myxlim){

	x <- data$mass
	y <- data$prob
	
	plot(x, y, xlim = myxlim, ylim = c(0,1))# ylim=c(0, max(y)), t="n", axes=FALSE, ann=FALSE)
	#axis(1)
	#axis(2, at=pretty(c(0, max(y))))
	#mtext("# Mutations", side=2, line=2.5)
	#mtext("P53_Human", side=3, line=0.5, adj=0, font=2)

	#This is the function that draws the lines down: 
	segments(x0=x, y0=y, y1=0, col=8)
	points(x, y, pch=21, col=1, bg=2)

}
