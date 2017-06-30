
LinRegBS<- function(dat, particles=1000) {

    # if no data supplied, use default
    if (missing(dat)){
		data(radiata)
		dat <- radiata
	}

    # more eloquent tests can be added
    stopifnot(nrow(dat) > 0,
              ncol(dat) == 2,
              colnames(dat) == c("x", "y"))
	
	res <- LinRegBS_cpp(as.matrix(dat),particles)			 
	
    invisible(res)
}
