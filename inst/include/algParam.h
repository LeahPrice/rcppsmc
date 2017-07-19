// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// algParam: A class that holds all of the algorithm parameters that can be adapted
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
// along with RInside.  If not, see <http://www.gnu.org/licenses/>.


//! \file 
//! \brief Algorithm parameter class
//!

#ifndef __SMC_ALGPARAM_H
#define __SMC_ALGPARAM_H 1.0

#include <RcppArmadillo.h>
#include <population.h>


namespace smc {

	///A class which contains the algorithm parameters.
	template <class Space, class Params> class algParam {
		
	protected:
		Params param;
		
	public:
		
		///Free the workspace allocated for the algorithm parameters
		virtual ~algParam() {
		}

		///Holder function for updates to be done before the MCMC step
		virtual void updateForMCMC(const population<Space> & pop, int nAccepted, int nResampled) {}
		
		///Holder function for updates to be done before the move step
		virtual void updateForMove(const population<Space> & pop) {}

		///Holder function for updates to be done at the end of each iteration
		virtual void updateEnd(const population<Space> & pop) {}
		
		///Function to get the values of the parameters
		Params GetParams(void){return param;}
	};
}


#endif
