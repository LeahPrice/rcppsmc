
LinRegLABS<- function(data, temperatures, particles=1000) {

    # if no data supplied, use default
    if (missing(data)) data <- getRadiataData()
    if (missing(temperatures)) temperatures <- getRadiataTemps()

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) == 2,
              colnames(data) == c("x", "y"),
			  nrow(temperatures) > 0)

    res <- LinRegLABS_cpp(as.matrix(data), as.matrix(temperatures), particles)

    invisible(res)
}

# simple convenience function, should probably make the data a
# data component of the package...
getRadiataData <- function() {
    file <- system.file("sampleData", "radiata-data.csv", package="RcppSMC")
    dat <- read.table(file, skip=1, header=FALSE, col.names=c("x","y"))
    invisible(dat)
}
# simple convenience function, should probably make the data a
# data component of the package...
getRadiataTemps <- function() {
    file <- system.file("sampleData", "radiata-temps.csv", package="RcppSMC")
    temperatures <- read.table(file, skip=1, header=FALSE, col.names=c("temperatures"))
    invisible(temperatures)
}