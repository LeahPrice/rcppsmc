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
//! \brief Random number generation.
//!
//! This file contains the definitions for the smc::rng and smc::rnginfo class.
//! It wraps the random number generation facilities provided by the GSL and provides a convenient interfaces to access several of its more commonly-used features.

#ifndef __SMC_ALGPARAM_H
#define __SMC_ALGPARAM_H 1.0

#include <RcppArmadillo.h>
#include <population.h>

namespace smc {

	///A class which contains the algorithm parameters.
	template <class Space> class algParam {
	private:
		///A base update function (updates that are always done)
		void updateForMCMCBase(const population<Space> & pop) {
		}
		
	public:
		///Initialise the random number generator using default settings
		algParam(void) {
		}
		
		///Free the workspace allocated for the algorithm parameters
		~algParam() {
			
		}

		///Holder function for additional updates that can be done (to be changed by the user)
		virtual void updateForMCMCExtra(const population<Space> & pop) {
			
		}
		
		
		///The function called from within the sampler object which combines the two
		void updateForMCMC(const population<Space> & pop) {
			updateForMCMCBase(pop);
			updateForMCMCExtra(pop);
		}


	};
}


#endif
