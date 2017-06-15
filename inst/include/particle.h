//   SMCTC: particle.hh
//
//   Copyright Adam Johansen, 2008-2009.
// 
//   This file is part of SMCTC.
//
//   SMCTC is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   SMCTC is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with SMCTC.  If not, see <http://www.gnu.org/licenses/>.
//

//! \file
//! \brief Class used to store and manipulate a single particle.
//!
//! This file contains the smc::particle class which is used internally and passed to move functions.

#ifndef __SMC_PARTICLE_HH
#define __SMC_PARTICLE_HH 1.0

#include <RcppArmadillo.h>

#include <float.h>
#include <limits>
#include <cmath>
namespace smc {
  /// A template class for the particles of an SMC algorithm
  template <class Space> class particle
    {
    private:
      /// Value of this particle
      std::vector<Space>    value;
      /// Natural logarithm of this particle's weight.
      arma::vec   logweight;

    public:
      particle();
      /// Constructor which initialises the particles value and weight.
      particle(std::vector<Space> sInit,arma::vec dLogWeight);
      /// The copy constructor performs a shallow copy.
      particle(const particle<Space> & pFrom);
      /// The assignment operator performs a shallow copy.
      particle<Space> & operator= (const particle<Space> & pFrom);

      ~particle();

      /// Returns the particle's value 
      std::vector<Space> const & GetValue(void) const {return value;}
      /// Returns the particle's value 
      Space GetValueN(int n) {return value[n];}
      /// Returns a pointer to the value to allow for more efficient changes
      std::vector<Space>* GetValuePointer(void) {return &value;}
      /// Returns a pointer to the value to allow for more efficient changes
      Space* GetValuePointerN(int n) {return &value[n];}
      /// Returns the particle's log weight.
      arma::vec GetLogWeight(void) const {return logweight;}
      /// Returns the particle's log weight.
      double GetLogWeightN(int n) const {return logweight(n);}
      /// Returns the particle's unnormalised weight.
      arma::vec GetWeight(void) const {return exp(logweight);}
      /// Returns the particle's unnormalised weight.
      double GetWeightN(int n) const {return exp(logweight(n));}
      
      /// \brief Sets the particle's value and weight explicitly
      ///
      /// \param sValue The particle value to use 
      /// \param dLogWeight The natural logarithm of the new particle weight
      void Set(std::vector<Space> sValue,arma::vec dLogWeight){value = sValue; logweight = dLogWeight;}
      /// \brief Sets the particle's value explicitly
      ///
      /// \param sValue The particle value to use
      void SetValue(const std::vector<Space> & sValue){value = sValue;}
      /// \brief Sets the particle's value explicitly
      ///
      /// \param sValue The particle value to use
      void SetValueN(const Space & sValue, int n){value[n] = sValue;}
      /// \brief Sets the particle's log weight explicitly
      ///
      /// \param dLogWeight The natural logarithm of the new particle weight
      void SetLogWeight(const arma::vec & dLogWeight) {logweight = dLogWeight;}
      /// \brief Sets the particles weight explicitly
      ///
      /// \param dWeight The new particle weight
      void SetWeight(const arma::vec & dWeight) {logweight = log(dWeight);}

      /// \brief Increase the log weight by a specified amount
      ///
      /// \param dIncrement The amount to add to the log weight.
      void AddToLogWeight(arma::vec dIncrement) { logweight += dIncrement;}
      /// \brief Multiply the weight by a specified factor
      ///
      /// \param dMultiplier The factor to multiply the weight by.
      void MultiplyWeightBy(arma::vec dMultiplier) { logweight += log(dMultiplier);}
  };


/// Create a particle with undefined value and weight NAN
  template <class Space>
    particle<Space>::particle()
    {	
		arma::vec logweight;
    }
  

  ///Copy constructor
  template <class Space>
  particle<Space>::particle(const particle<Space> & pFrom)
  {
    value = pFrom.value;
    logweight = pFrom.logweight;
  }
  
  /// Create a particle with value sInit and log weight dLogWeight 
  /// \param sInit The initial value of the particle
  /// \param dLogWeight The initial value of the natural logarithm of the particle weight
  template <class Space>
    particle<Space>::particle(std::vector<Space> sInit, arma::vec dLogWeight)
    {
      value = sInit;
      logweight =dLogWeight;
    }

  /// Dispose of a particle which is no longer required
  template <class Space>
  particle<Space>::~particle()
  {
  }

  /// Copy the values of pFrom to the values of this to set this particle identical to pFrom in a deep
  /// copy sense.
  template <class Space>
  particle<Space> & particle<Space>::operator= (const particle<Space> & pFrom)
  {	  
    this->value = pFrom.value;
    this->logweight = pFrom.logweight;
    return *this;
  }
}

#endif
