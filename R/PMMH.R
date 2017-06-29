
PMMH<- function(data, particles=5000, iterations=10000) {

    # if no data supplied, use default
    if (missing(data)) data <- getPMMHData()

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) == 1,
              colnames(data) == c("x"))

    res <- PMMH_cpp(as.matrix(data), particles, iterations)

    invisible(res)
}

# simple convenience function, should probably make the data a
# data component of the package...
getPMMHData <- function() {
    file <- system.file("sampleData", "PMMH-data.csv", package="RcppSMC")
    dat <- read.table(file, skip=1, header=FALSE, col.names=c("x"))
    invisible(dat)
}