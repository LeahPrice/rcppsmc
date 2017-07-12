// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// pflineart.cpp: Rcpp wrapper for SMC library -- first example of Johansen (2009)
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2012         Dirk Eddelbuettel and Adam Johansen
// Copyright (C) 2017         Dirk Eddelbuettel, Adam Johansen and Leah South
//
// This file is part of RcppSMC.
//
// RcppSMC is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppSMC is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppSMC.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>

#include "smctc.h"
#include "pflineart.h"
#include <typeinfo>
#include <cxxabi.h> // for demangling the type name

#include <cstdio> 
#include <cstdlib>
#include <cstring>

namespace pflineart {
	const double var_s0 = 4;
	const double var_u0 = 1;
	const double var_s  = 0.02;
	const double var_u  = 0.001;

	const double scale_y = 0.1;
	const double nu_y = 10.0;
	const double Delta = 0.1;
}

using namespace std;
using namespace arma;
using namespace pflineart;


// pf() function callable from R via Rcpp:: essentially the same as main() from pf.cc 
// minor interface change to pass data down as matrix, rather than a filename
// [[Rcpp::export]]
Rcpp::List pfLineartBS_cpp(arma::mat data, unsigned long inlNumber, bool useF, Rcpp::Function f) { 	

	long lIterates;

	try {
		
		lNumber = inlNumber;
		
		lIterates = data.n_rows;
		y.x_pos = data.col(0);
		y.y_pos = data.col(1);
		
		///Initialise and run the sampler 
		smc::sampler<cv_state> Sampler(lNumber, HistoryType::NONE);  
		smc::moveset<cv_state> Moveset(fInitialise, fMove, NULL);
		
		
		//Sampler.SetSMCParams(ResampleType::RESIDUAL, 0.5);
		Sampler.SetResampParams(ResampleType::RESIDUAL, 0.999);
		Sampler.SetMoveSet(Moveset);
		Sampler.Initialise();
		
		Rcpp::NumericVector Xm(lIterates), Xv(lIterates), Ym(lIterates), Yv(lIterates), ESS(lIterates);
		
		Xm(0) = Sampler.Integrate(integrand_mean_x, NULL);
		Xv(0) = Sampler.Integrate(integrand_var_x, (void*)&Xm(0));
		Ym(0) = Sampler.Integrate(integrand_mean_y, NULL);
		Yv(0) = Sampler.Integrate(integrand_var_y, (void*)&Ym(0));
		ESS(0) = Sampler.GetESS();
		
		for(int n=1; n < lIterates; ++n) {
			Sampler.Iterate();
			
			Xm(n) = Sampler.Integrate(integrand_mean_x, NULL);
			Xv(n) = Sampler.Integrate(integrand_var_x, (void*)&Xm(n));
			Ym(n) = Sampler.Integrate(integrand_mean_y, NULL);
			Yv(n) = Sampler.Integrate(integrand_var_y, (void*)&Ym(n));
			ESS(n) = Sampler.GetESS();
			
			if (useF) f(Xm, Ym);
		}
		
		//Rcpp::Rcout << Sampler << std::endl;
		//Sampler.StreamParticle(Rcpp::Rcout,0);
		
		double logNC = Sampler.GetLogNCPath();
		
		return Rcpp::List::create(Rcpp::Named("Xm") = Xm,
		Rcpp::Named("Xv") = Xv,
		Rcpp::Named("Ym") = Ym,
		Rcpp::Named("Yv") = Yv,
		Rcpp::Named("ESS") = ESS,
		Rcpp::Named("logNC") = logNC);
	}
	catch(smc::exception  e) {
		Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o
	}
	return R_NilValue;          	// to provide a return 
}

namespace pflineart {
	double integrand_mean_x(const cv_state& s, void *){ return s.x_pos;}
	double integrand_mean_y(const cv_state& s, void *){ return s.y_pos;}

	double integrand_var_x(const cv_state& s, void* vmx){
		double* dmx = (double*)vmx;
		double d = (s.x_pos - (*dmx));
		return d*d;
	}

	double integrand_var_y(const cv_state& s, void* vmy)
	{
		double* dmy = (double*)vmy;
		double d = (s.y_pos - (*dmy));
		return d*d;
	}

#include <iostream>
#include <cmath>
	//#include <gsl/gsl_randist.h>

	using namespace std;

	double logLikelihood(long lTime, const cv_state & value)
	{
		return - 0.5 * (nu_y + 1.0) * (log(1 + pow((value.x_pos - y.x_pos[lTime])/scale_y,2) / nu_y) + log(1 + pow((value.y_pos - y.y_pos[lTime])/scale_y,2) / nu_y));
	}

	void fInitialise(cv_state & value, double & logweight)
	{ 
		value.x_pos = R::rnorm(0,sqrt(var_s0));
		value.y_pos  = R::rnorm(0,sqrt(var_s0));
		value.x_vel  = R::rnorm(0,sqrt(var_u0));
		value.y_vel  = R::rnorm(0,sqrt(var_u0));
		logweight = - 0.5 * (nu_y + 1.0) * (log(1 + pow((value.x_pos - y.x_pos[0])/scale_y,2) / nu_y) + log(1 + pow((value.y_pos - y.y_pos[0])/scale_y,2) / nu_y));
	}

	void fMove(long lTime, cv_state & value, double & logweight)
	{
		value.x_pos += value.x_vel * Delta + R::rnorm(0,sqrt(var_s));
		value.x_vel += R::rnorm(0,sqrt(var_u));
		value.y_pos += value.y_vel * Delta + R::rnorm(0,sqrt(var_s));
		value.y_vel += R::rnorm(0,sqrt(var_u));

		logweight += logLikelihood(lTime, value);
	}
}



namespace std {
	/// Produce a human readable display of an the particle values using the standard stream operators

	/// \param os The output stream to which the display should be made.
	/// \param p  The particle which is to be displayed.
	std::ostream & operator << (std::ostream & os, cv_state & value)
	{
		double x_pos = value.x_pos;
		double y_pos = value.y_pos;
		double x_vel = value.x_vel;
		double y_vel = value.y_vel;
		os << "Position: (" << x_pos << ", " << y_pos << ") Velocity: (" <<  x_vel << ", " << y_vel << ")";
		return os;
	}
}