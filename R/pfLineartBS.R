
pfLineartBS<- function(dat, particles=1000, plot=FALSE, onlinePlot) {

    # if no data supplied, use default
    if (missing(dat)){
		data(pfdata)
		dat <- pfdata
	}

    if (missing(onlinePlot)) {
        useOnline <- FALSE
        onlinePlot <- function() { NULL }
    } else {
        useOnline <- TRUE
        # set up graphics window
        dev.new(width=3,height=3)
        par(mar=c(3,3,1,1),cex=0.8, pch=19, ask=FALSE)
    }

    # more eloquent tests can be added
    stopifnot(nrow(dat) > 0,
              ncol(dat) == 2,
              colnames(dat) == c("x", "y"),
              class(onlinePlot) == "function")
				 
    res <- pfLineartBS_cpp(as.matrix(dat), particles, useOnline, onlinePlot)

    if (plot) {
        ## plot 5.1 from vignette / paper
        with(dat, plot(x, y, col="red"))
        with(res, lines(Xm, Ym, lty="dashed"))
    }

    invisible(res)
}

pfLineartRange <- function(rrng)
{
   min <- rrng[1]
   max <- rrng[2]

   if(min > 0) {
      rmin = exp(floor(log(min)))
   } else if (min < 0) {
      rmin = -exp(ceiling(log(-min)))
   } else {
      rmin = 0
   }

   if(max > 0) {
      rmax = exp(ceiling(log(max)))
   } else if (max < 0){
      rmax = exp(floor(log(-max)))
   } else {
      rmax = 0;
   }

   invisible(c(rmin,rmax))
}

pfLineartBSOnlinePlot <- function(xm, ym) {
    plot(xm, ym, xlim = pfLineartRange(range(xm)), ylim=pfLineartRange(range(ym)))
    # FIXME sleep time should also be a variable
    Sys.sleep(0.05)
}

