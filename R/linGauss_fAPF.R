linGauss_fAPF<- function(data, initial = c(0.9,0.3,0.7,0.1,0.2,0.6,0.4,0.1,0.1,0.3,0.1,0.2,0.5,0.2,0), particles=5000, iterations=100000) {

    if (missing(data)) {
        A <- diag(5)
        A[upper.tri(A,diag=T)] <- c(0.9,0.3,0.7,0.1,0.2,0.6,0.4,0.1,0.1,0.3,0.1,0.2,0.5,0.2,0)
        A <- t(A)
        warning("data argument contained no data, using data simulated from the default model.")
        data <- simLinGauss(len=100,A=A)$data
    }

    res <- linGauss_fAPF_impl(as.matrix(data), as.matrix(initial), particles, iterations)

    invisible(res)
}
