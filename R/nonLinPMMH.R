
nonLinPMMH<- function(data, particles=5000, iterations=10000) {

    if (missing(data)) {
         warning("data argument contained no data, using data simulated from the model.")
         data <- simNonlinPMMH(len=500)$data
    }

    res <- nonLinPMMH_cpp(as.matrix(data), particles, iterations)
	
	# Need to add some plots here - probably bivariate density plots, histogram and trace plots
	# like in Figure 4(c) of Andrieu et al. (2010). Also probably some autocorrelation plots
	# of the parameters, like in Figures 5(b) and 5(d).

    invisible(res)
}
