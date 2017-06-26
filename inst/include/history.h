// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// history.h: Rcpp integration of SMC library -- sampler history 
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
//! \brief Classes and function related to the history of the sampler.
//!
//! This file contains template definitions for the classes used to store the history of an SMCTC sampler.
//! It defines smc::history, smc::historyelement and smc::history.

#ifndef __SMC_HISTORY_HH
#define __SMC_HISTORY_HH 1.0

#include <RcppArmadillo.h>

namespace smc {
  /// The historyflags class holds a set of flags which describe various properties of the population at a given time.
  class historyflags
  {
  private:
    /// true if the population was resampled during the described iteration.
    unsigned int Resampled : 1;
  public:
    ///Create a new set of history flags corresponding to the specified properties
    historyflags(int wasResampled);

    ///This function returns true if the flag set indicates that the ensemble was resampled during the described iteration.
    int WasResampled(void) const {return Resampled;}
  };
  
  
  /// A template class for the elements of a linked list to be used for the history of the sampler.
  template <class Population>class historyelement
  {
  private:
    long number; //!< The number of particles (presently redundant as this is not a function of iteration)
    int nAccepted; //!< Number of MCMC moves accepted during this iteration.
    Population value; //!< The particles themselves (values and weights)
    historyflags flags; //!< Flags associated with this iteration.
    
  public:
    /// The null constructor creates an empty history element.
    historyelement();
    /// A constructor with four arguments initialises the particle set.
    historyelement(long lNumber, Population New, int nAccepts, historyflags hf);

    /// The destructor tidies up.
    ~historyelement();

    /// Returns the effective sample size of this particle generation.
    double GetESS(void) const;
    /// Returns the flags
    historyflags GetFlags(void) const {return flags;}
    /// Returns the number of particles present.
    long GetNumber(void) const {return number;} 
    /// Returns a pointer to the current particle set.
    Population * GetPointers(void) const { return &value; }
	/// Returns a pointer to the current particle set.
    Population & GetRefs(void) { return value; }
    /// Integrate the supplied function according to the empirical measure of the particle ensemble.
    long double Integrate(long lTime, double (*pIntegrand)(long,const Population&,long,void*), void* pAuxiliary) const;
    /// Sets the particle set to the specified values.  
    void Set(long lNumber, const Population &New, int inAccepted, const historyflags &histflags){number = lNumber; value = New; nAccepted = inAccepted; flags = histflags;};
    /// Returns the number of MCMC moves accepted during this iteration.
    int AcceptCount(void) const {return nAccepted; }
    /// Returns true if the particle set 
    int WasResampled(void) const {return flags.WasResampled(); }

  };

  template <class Population>
  historyelement<Population>::historyelement(): flags(0)
  {
    number = 0;
    nAccepted = 0;
  }


  /// \param lNumber The number of particles present in the particle generation
  /// \param New    The array of particles which are present in the particle generation
  /// \param nAccepts The number of MCMC moves that were accepted during this particle generation
  /// \param hf      The historyflags associated with the particle generation

  template <class Population>
  historyelement<Population>::historyelement(long lNumber, Population New, int nAccepts, historyflags hf) :
    flags(hf)
  {
    number = lNumber;
    Population value = New;
    nAccepted = nAccepts;
	flags = hf;
  }

  template <class Population>
  historyelement<Population>::~historyelement(void)
  {
  }

  template <class Population>
  double historyelement<Population>::GetESS(void) const
  {
	double sum = arma::sum(exp(value.GetLogWeight()));
	double sumsq = arma::sum(exp(2.0*value.GetLogWeight()));
  return expl(-log(sumsq) + 2*log(sum));
  }

  /// \param lTime The timestep at which the integration is to be carried out
  /// \param pIntegrand The function which is to be integrated
  /// \param pAuxiliary A pointer to additional information which is passed to the integrand function

  template <class Population>
  long double historyelement<Population>::Integrate(long lTime, double (*pIntegrand)(long,const Population&,long,void*), void* pAuxiliary) const
  {
    long double rValue = 0;
	long double wSum = 0;
	for(long i =0; i < number; i++)
	{
		rValue += expl(value.GetLogWeightN(i)) * (long double)pIntegrand(lTime, value,i, pAuxiliary); // NEEDS FIXING TO HAVE RIGHT INPUT
		wSum  += expl(value.GetLogWeightN(i));
	}
  
	rValue /= wSum;
	return rValue;
  }
  
  
  
  
  
  
  /// A template class for the history associated with a population evolving in SMC.

  ///  The history is a template class which should have an associated class type corresponding to
  ///    a _population_ of the desired type, not the type itself.
  ///
  ///    Essentially, this is implemented as a doubly linked list. 


  template <class Population> class history
    {
    private:
      ///The first time step
	  std::list<historyelement<Population> > hist;

    public:
      ///The argument free constructor creates an empty list.
      history();
	  
	  ///The standard list operations
	  const std::list<historyelement<Population> > & GetHistory(void) const {return hist;}
	  unsigned long size(void) const {return hist.size();}
	  void clear(void) {hist.clear();}
	  void push_back(historyelement<Population> & histel) {hist.push_back(histel);}
	  void pop_back(void) {hist.pop_back();}
	  const historyelement<Population> & back(void) const {return hist.back();}

      /// Returns the effective sample size of the specified particle generation.
      double GetESS(long lGeneration) const;
      ///Integrate the supplied function over the path of the particle ensemble.
      double IntegratePathSampling(double (*pIntegrand)(long,const Population&,long,void*), double (*pWidth)(long,void*), void* pAuxiliary) const;
      double IntegratePathSamplingFinalStep(double (*pIntegrand)(long,const Population&,long,void*), void* pAuxiliary) const;

      ///Output a vector indicating the number of accepted MCMC moves at each time instance
      void OstreamMCMCRecordToStream(std::ostream &os) const;
      ///Output a 0-1 value vector indicating the times at which resampling occured to an output stream
      void OstreamResamplingRecordToStream(std::ostream &os) const;

      ///Display the list of particles in a human readable form.
      //  void StreamPopulations(std::ostream & os);
    };

  /// This constructor simply sets the root and leaf pointers to NULL and the list length to zero.
  template <class Population>
  history<Population>::history()
  {
	  
  }
   /// Returns the effective sample size of the specified particle generation.
  template <class Population>
  double  history<Population>::GetESS(long lGeneration) const
  {
    typename std::list<historyelement<Population> >::const_iterator it = hist.begin();
	std::advance(it,lGeneration);
    return it->GetESS(); 
  }



  /// This function records the MCMC acceptance history to the specified output stream as a list of
  /// the number of moves accepted at each time instant.
  ///
  /// \param os The output stream to send the data to.
  template <class Population>
  void history<Population>:: OstreamMCMCRecordToStream(std::ostream &os) const
  {
	os << "Accepted MCMC proposals history:" << std::endl;
    os << "======================" << std::endl;
	for(typename std::list<historyelement<Population> >::const_iterator it = hist.begin(); it!=hist.end(); it++){
		os << it->AcceptCount() << std::endl;
    }
  }
  /// This function records the resampling history to the specified output stream as a 0-1 valued list which takes
  /// the value 1 for those time instances when resampling occured and 0 otherwise.
  ///
  /// \param os The output stream to send the data to.
  template <class Population>
  void history<Population>:: OstreamResamplingRecordToStream(std::ostream &os) const
  {
	os << "Resampling history:" << std::endl;
    os << "======================" << std::endl;
	os << "Flag\t" << "ESS\t" << std::endl;
	for(typename std::list<historyelement<Population> >::const_iterator it = hist.begin(); it!=hist.end(); it++){ 
	  if(it->WasResampled())
			os << "1\t";
      else
			os << "0\t";

		os << it->GetESS() << std::endl;
    }
  }

  /// This function performs a trapezoidal integration of the type which is useful when using path sampling to estimate the
  /// normalising constant of a potential function in those cases where a sequence of distributions is produced by deforming
  /// the initial distribution by a sequence of progressively more influential potential functions which are proportional
  /// to the density of some other distribution with respect to the starting distribution.
  ///
  /// The integrand is integrated at every time point in the population history. The results of this integration are
  /// taken to be point-evaluations of the path sampling integrand which are spaced on a grid of intervals given by the
  /// width function. The path sampling integral is then calculated by performing a suitable trapezoidal integration and
  /// the results of this integration is returned.
  ///
  /// pAuxiliary is passed to both of the user specified functions to allow the user to pass additional data to either or
  /// both of these functions in a convenient manner. It is safe to use NULL if no such data is used by either function.
  ///
  /// \param pIntegrand  The function to integrated. The first argument is evolution time, the second the population at which the function is to be evaluated, the third is the particle index and the final argument is always pAuxiliary.
  /// \param pWidth      The function which returns the width of the path sampling grid at the specified evolution time. The final argument is always pAuxiliary
  /// \param pAuxiliary  A pointer to auxiliary data to pass to both of the above functions

  template <class Population>
  double history<Population>::IntegratePathSampling(double (*pIntegrand)(long,const Population&,long,void*), double (*pWidth)(long,void*), void* pAuxiliary) const
  {
	long lTime = 1;
	long double rValue = 0.0;
    //for(typename std::list<historyelement<Population> >::const_iterator it = ++hist.begin(); it!=--hist.end(); it++){
    for(typename std::list<historyelement<Population> >::const_iterator it = ++hist.begin(); it!=hist.end(); it++){
		rValue += it->Integrate(lTime, pIntegrand, pAuxiliary) * (long double)pWidth(lTime,pAuxiliary);
		lTime++;
    }	
	return ((double)rValue);
  }

  template <class Population>
  double history<Population>::IntegratePathSamplingFinalStep(double (*pIntegrand)(long,const Population&,long,void*), void* pAuxiliary) const
  {
    return hist.back().Integrate(hist.size()-1,pIntegrand,pAuxiliary);
  }

}
  
#endif
