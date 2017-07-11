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

///Specifiers for various resampling algorithms:
namespace ResampleType
{
	enum Enum {MULTINOMIAL = 0, 
		RESIDUAL, 
		STRATIFIED, 
		SYSTEMATIC };
};

namespace smc {

	///A class which contains the algorithm parameters.
	template <class Space> class algParam {
	private:
		long myN;
		ResampleType::Enum myResample;
		double myResampleThreshold;
		
		///A base update function (updates that are always done)
		void updateForMCMCBase(const population<Space> & pop) {}
		
		///A base update function (updates that are always done)
		void updateForMoveBase(const population<Space> & pop) {}
		
		///A base update function (updates that are always done)
		void updateEndBase(const population<Space> & pop) {}
		
		
		
	public:
		///Initialise the algorithm parameters
		algParam(void){};
		
		///Initialise the algorithm parameters
		algParam(long n) : myN(n), myResample(ResampleType::STRATIFIED), myResampleThreshold(0.5*n){};
		
		///Initialise the algorithm parameters
		algParam(long n, ResampleType::Enum restype, double resthresh) : myN(n), myResample(restype), myResampleThreshold(resthresh){};
		
		void SetN(long n){ myN = n;}
		
		long GetN(void){return myN;}
		
		void SetResThresh(double dThreshold){
		if(dThreshold < 1)
		 myResampleThreshold = dThreshold * myN;
		else
		myResampleThreshold = dThreshold;
		}
		
		double GetResThresh(void){return myResampleThreshold;}
		
		void SetResample(ResampleType::Enum inRes){ myResample = inRes;}
		
		ResampleType::Enum GetResample(void){return myResample;}
		
		///Free the workspace allocated for the algorithm parameters
		virtual ~algParam() {
			
		}

		///Holder function for additional updates that can be done (to be changed by the user)
		virtual void updateForMCMCExtra(const population<Space> & pop) {}
		
		///Holder function for additional updates that can be done (to be changed by the user)
		virtual void updateForMoveExtra(const population<Space> & pop) {}

		///Holder function for additional updates that can be done (to be changed by the user)
		virtual void updateEndExtra(const population<Space> & pop) {}
		
		///The function called from within the sampler object which combines the two
		void updateForMCMC(const population<Space> & pop) {
			updateForMCMCBase(pop);
			updateForMCMCExtra(pop);
		}

		///The function called from within the sampler object which combines the two
		void updateForMove(const population<Space> & pop) {
			updateForMoveBase(pop);
			updateForMoveExtra(pop);
		}
		///The function called from within the sampler object which combines the two
		void updateEnd(const population<Space> & pop) {
			updateEndBase(pop);
			updateEndExtra(pop);
		}
	};
}


#endif
