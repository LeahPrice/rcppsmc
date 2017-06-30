
LinRegLABS<- function(dat, temperatures, particles=1000) {

    # if no data supplied, use default
    if (missing(dat)){
		data(radiata)
		dat <- radiata
	}
    if (missing(temperatures)) temperatures <- seq(0,1,0.05)^5

    # more eloquent tests can be added
    stopifnot(nrow(dat) > 0,
              ncol(dat) == 2,
              colnames(dat) == c("x", "y"),
			  nrow(temperatures) > 0)

    res <- LinRegLABS_cpp(as.matrix(dat), as.matrix(temperatures), particles)

    invisible(res)
}