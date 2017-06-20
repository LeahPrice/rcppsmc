// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// particle.h: Rcpp integration of SMC library -- storing and manipulating particles 
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2017		  Adam Johansen, Dirk Eddelbuettel and Leah South
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
//

//! \file
//! \brief Class used to store and manipulate a population of particles.
//!
//! This file contains the smc::population class which is used internally and passed to move functions.

#ifndef __SMC_POPULATION_HH
#define __SMC_POPULATION_HH 1.0

#include <RcppArmadillo.h>

#include <float.h>
#include <limits>
#include <cmath>
namespace smc {	
  /// A template class for the particles of an SMC algorithm
  template <class Space> class population
    {
    private:
      /// Value of the particles
      std::vector<Space>    value;
      /// Natural logarithm of this particles' weights.
      arma::vec   logweight;

    public:
      population();
      /// Constructor which initialises the particle values and weights.
      population(const std::vector<Space> &sInit,const arma::vec & dLogWeight);
      /// The copy constructor performs a shallow copy.
      population(const population<Space> & pFrom);
      /// The assignment operator performs a shallow copy.
      population<Space> & operator= (const population<Space> & pFrom);

      ~population();

      /// Returns the population values 
      std::vector<Space> const & GetValue(void) const {return value;}
      /// Returns the value of the nth particle in the population
      Space const & GetValueN(int n) {return value[n];}
      /// Returns the value of the nth particle in the population
      Space & GetValueRefN(int n) {return value[n];}
      /// Returns a pointer to the values to allow for more efficient changes
      std::vector<Space>* GetValuePointer(void) {return &value;}
      /// Returns a pointer to the value to allow for more efficient changes
      Space* GetValuePointerN(int n) {return &value[n];}
      /// Returns the particles' log weights.
      arma::vec GetLogWeight(void) const {return logweight;}
      /// Returns the nth particle's log weight.
      double GetLogWeightN(int n) const {return logweight(n);}
      /// Returns the nth particle's log weight.
      double & GetLogWeightRefN(int n) {return logweight(n);}
      /// Returns the particles' unnormalised weights.
      arma::vec GetWeight(void) const {return exp(logweight);}
      /// Returns the nth particle's unnormalised weight.
      double GetWeightN(int n) const {return exp(logweight(n));}
      
      /// \brief Sets the particle values and weight explicitly
      ///
      /// \param sValue The particle values to use 
      /// \param dLogWeight The natural logarithm of the new particle weights
      void Set(const std::vector<Space> &sValue,const arma::vec & dLogWeight){value = sValue; logweight = dLogWeight;}
      /// \brief Sets the particle's value explicitly
      ///
      /// \param sValue The particle values to use
      void SetValue(const std::vector<Space> & sValue){value = sValue;}
      /// \brief Sets the particles' values explicitly
      ///
      /// \param sValue The particle values to use
      void SetValueN(const Space & sValue, int n){value[n] = sValue;}
      /// \brief Sets the particle log weights explicitly
      ///
      /// \param dLogWeight The natural logarithm of the new particle weights
      void SetLogWeight(const arma::vec & dLogWeight) {logweight = dLogWeight;}
      /// \brief Sets the particle weights explicitly
      ///
      /// \param dWeight The new particle weights
      void SetWeight(const arma::vec & dWeight) {logweight = log(dWeight);}

      /// \brief Increase the log weights by a specified amount
      ///
      /// \param dIncrement The amount to add to the log weights.
      void AddToLogWeight(const arma::vec & dIncrement) { logweight += dIncrement;}
      /// \brief Multiply the weights by a specified factor
      ///
      /// \param dMultiplier The factor to multiply the weights by.
      void MultiplyWeightBy(const arma::vec & dMultiplier) { logweight += log(dMultiplier);}
  };


/// Create a particle with undefined values and weights NAN
  template <class Space>
    population<Space>::population()
    {	
		arma::vec logweight;
    }
  

  ///Copy constructor
  template <class Space>
  population<Space>::population(const population<Space> & pFrom)
  {
    value = pFrom.value;
    logweight = pFrom.logweight;
  }
  
  /// Create particles with values sInit and log weights dLogWeight 
  /// \param sInit The initial values of the particles
  /// \param dLogWeight The initial values of the natural logarithm of the particle weights
  template <class Space>
    population<Space>::population(const std::vector<Space> & sInit, const arma::vec & dLogWeight)
    {
      value = sInit;
      logweight =dLogWeight;
    }

  /// Dispose of particles which are no longer required
  template <class Space>
  population<Space>::~population()
  {
  }

  /// Copy the values of pFrom to the values of this to set this particle identical to pFrom in a deep
  /// copy sense.
  template <class Space>
  population<Space> & population<Space>::operator= (const population<Space> & pFrom)
  {	  
    this->value = pFrom.value;
    this->logweight = pFrom.logweight;
    return *this;
  }
}

#endif
