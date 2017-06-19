
radiataBS<- function(data, particles=1000) {

    # if no data supplied, use default
    if (missing(data)) data <- getradiataBSData()

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) == 2,
              colnames(data) == c("x", "y"))
	
	res <- radiataBS_cpp(as.matrix(data),particles)			 
				 

    invisible(res)
}

# simple convenience function, should probably make the data a
# data component of the package...
getradiataBSData <- function() {
    file <- system.file("sampleData", "radiata-data.csv", package="RcppSMC")
    dat <- read.table(file, skip=1, header=FALSE, col.names=c("x","y"))
    invisible(dat)
}