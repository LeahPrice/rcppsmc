\name{NEWS}
\title{News for Package \pkg{RcppSMC}}

\section{Changes in RcppSMC version 0.1.8 (2017-06-16)}{
  \itemize{
	\item Switching to armadillo vectors. The particle class is now
		  on the whole population of particles level.
    \item Added an example \code{radiataPPBS} for a likelihood annealing
		  version of the radiata pine linear regression problem.
  }
}

\section{Changes in RcppSMC version 0.1.7 (2017-06-06)}{
  \itemize{
    \item Added an example \code{radiataBS} for a linear regression problem
	\item Added code to estimate the ratio of normalising constants between the initial and final targets
  }
}

\section{Changes in RcppSMC version 0.1.6 (2017-xx-yy)}{
  \itemize{
    \item Also use \code{.registration=TRUE} in \code{useDynLib} in
    \code{NAMESPACE} 
  }
}
  
\section{Changes in RcppSMC version 0.1.5 (2017-03-03)}{
  \itemize{
    \item Correct \code{.Call} to use \code{PACKAGE=} argument
    \item \code{DESCRIPTION}, \code{NAMESPACE}, \code{README.md} changes
    to comply with current \code{R CMD check} levels
    \item Added file \code{init.c} with calls to \code{R_registerRoutines()}
    \code{and R_useDynamicSymbols()}
    \item Updated \code{.travis.yml} file for continuous integration
  }
}

\section{Changes in RcppSMC version 0.1.4 (2014-05-02)}{
  \itemize{
    \item Added missing delete operator to destructor in \code{sampler} 
  }
}

\section{Changes in RcppSMC version 0.1.3 (2014-05-01)}{
  \itemize{
    \item Bugfix in \code{Sampler.iterate()} for memory overrun detected by valgrind
  }
}

\section{Changes in RcppSMC version 0.1.2 (2014-02-06)}{
  \itemize{
    \item Updated for \CRANpkg{Rcpp} 0.11.0 with explicit \code{importFrom}
    in \code{NAMESPACE} and corresponding versioned \code{Imports:} in \code{DESCRIPTION};
    also removed linking instruction from \code{src/Makevars} as no it
    is longer needed with this new \CRANpkg{Rcpp} version
    \item Addded GitHub / Travis CI support
    \item Use more portable \code{dev.new()} rather than \code{x11()} in \code{pfLinearBS.R}
    \item Applied some corrections to \code{pfNonlinBS.R} example
    \item Converted \code{NEWS} to \code{NEWS.Rd}
  }
}

\section{Changes in RcppSMC version 0.1.1 (2012-05-14)}{
  \itemize{
    \item Minor \code{g++-4.7} build fix of using \code{std::max()} explicitly
  }
}

\section{Changes in RcppSMC version 0.1.0 (2012-03-21)}{
  \itemize{
    \item Version 0.1.0 and initial CRAN upload
  }
}
