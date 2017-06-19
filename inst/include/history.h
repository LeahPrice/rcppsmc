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
    int WasResampled(void) {return Resampled;}
  };

  
  
  
  
  /// A template class for the elements of a linked list to be used for the history of the sampler.
  template <class Population>class historyelement
  {
  private:
    long number; //!< The number of particles (presently redundant as this is not a function of iteration)
    int nAccepted; //!< Number of MCMC moves accepted during this iteration.
    Population value; //!< The particles themselves (values and weights)
    historyflags flags; //!< Flags associated with this iteration.
    historyelement<Population> *pLast; //!< The parent of this node.
    historyelement<Population> *pNext; //!< The child of this node.
    
  public:
    /// The null constructor creates an empty history element.
    historyelement();
    /// A constructor with four arguments initialises the particle set.
    historyelement(long lNumber, const Population & pNew, int nAccepts, historyflags hf);
    /// A constructor with six arguments initialises the whole structure.
    historyelement(long lNumber, const Population & pNew, historyelement<Population>* pParent, historyelement<Population>* pChild, int nAccepts, historyflags  hf);

    /// The destructor tidies up.
    ~historyelement();

    /// Returns the effective sample size of this particle generation.
    double GetESS(void) const;
    /// Returns the flags
    historyflags GetFlags(void) const{return flags;}
    /// Returns the parent of this element.
    historyelement<Population> * GetLast(void) const { return pLast; }
    /// Returns the child of this element.
    historyelement<Population> * GetNext(void) const { return pNext; }
    /// Returns the number of particles present.
    long GetNumber(void) const {return number;} 
    /// Returns a pointer to the current particle set.
    Population const & GetValues(void) const { return value; }
    /// Returns a pointer to the current particle set.
    Population * GetPointers(void) { return &value; }
    /// Add a new history element with the specified values after the current one and maintain the list structure.
    void InsertAfter(long lNumber, const Population & pNew);
    /// Integrate the supplied function according to the empirical measure of the particle ensemble.
    long double Integrate(long lTime, double (*pIntegrand)(long,const Population &,void*), void* pAuxiliary);
    /// Sets the elements parent.
    void SetLast(historyelement<Population>* pParent) {pLast = pParent; }
    /// Sets the elements child.
    void SetNext(historyelement<Population>* pChild) {pNext = pChild; }
    /// Sets the particle set to the specified values.    
    void SetValue(long lNumber, const Population & pNew);

    /// Returns the number of MCMC moves accepted during this iteration.
    int AcceptCount(void) {return nAccepted; }
    /// Returns true if the particle set 
    int WasResampled(void) {return flags.WasResampled(); }
    /// \brief The left shift operator returns the element a number of positions prior to this one.
    ///
    /// \param ulCount The number of positions by which to shift.
    /// \return The element a number of positions before this one.

    historyelement<Population> operator<<(unsigned long ulCount)
    {
      if(ulCount)
	return this->pLast << (--ulCount);
	else
	  return *this;
    }

    ///\brief The right shift operator returns the element a number of positions after this one.
    ///
    /// \param ulCount The number of positions by which to shift.
    /// \return The right shift operator returns the element a number of positions after this one.
    historyelement<Population> operator>>(unsigned long ulCount)
    {
      if(ulCount)
	return this->pNext << (--ulCount);
      else
	return *this;
    }

  };

  template <class Population>
  historyelement<Population>::historyelement()
  {
    number = 0;
    value = NULL;
    nAccepted = 0;
    pLast = NULL;
    pNext = NULL;
  }


  /// \param lNumber The number of particles present in the particle generation
  /// \param pNew    The array of particles which are present in the particle generation
  /// \param nAccepts The number of MCMC moves that were accepted during this particle generation
  /// \param hf      The historyflags associated with the particle generation

  template <class Population>
  historyelement<Population>::historyelement(long lNumber, const Population & pNew, int nAccepts, historyflags hf) :
    flags(hf)
  {
    number = lNumber;
    // value = new Population[number];
    // for(int i = 0; i < number; i++)
    //   value[i] = pNew[i];
  
	Population value;
	value = pNew;

    nAccepted = nAccepts;
    pLast = NULL;
    pNext = NULL;
  }

  /// \param lNumber The number of particles present in the particle generation
  /// \param pNew    The array of particles which are present in the particle generation
  /// \param pParent A pointer to the previous particle generation
  /// \param pChild  A pointer to the next particle generation
  /// \param nAccepts The number of MCMC moves that were accepted during this particle generation
  /// \param hf      The historyflags associated with the particle generation
  template <class Population>
  historyelement<Population>::historyelement(long lNumber, const Population & pNew,
					   historyelement<Population>* pParent, historyelement<Population>* pChild,
					   int nAccepts, historyflags hf) :
    flags(hf)
  {
    number = lNumber;
  
	Population value;
	value = pNew;

    nAccepted = nAccepts;
    pLast = pParent;
    pNext = pChild;
  }

  template <class Population>
  historyelement<Population>::~historyelement(void)
  {
    // if(value)
    //   delete [] value;
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
  long double historyelement<Population>::Integrate(long lTime, double (*pIntegrand)(long,const Population &,void*), void* pAuxiliary)
  {
	long double rValue = 0;
	long double wSum = 0;
	for(int i =0; i < number; i++)
	{
		rValue += expl(value.GetLogWeightN(i)) * (long double)pIntegrand(lTime, value.GetValueN(i), pAuxiliary); // NEEDS FIXING TO HAVE RIGHT INPUT
		// rValue += expl(value[i].GetLogWeight()) * (long double)pIntegrand(lTime, value[i], pAuxiliary);
		wSum  += expl(value.GetLogWeightN(i));
	}
  
	rValue /= wSum;
	return rValue;
  }

  /// \param lNumber The number of particles in the generation to be inserted
  /// \param pNew The value of the particle generation to be inserted

  template <class Population>
  void historyelement<Population>::InsertAfter(long lNumber, const Population & pNew)
  {
    pNext = new historyelement<Population>(lNumber, pNew, this, pNext);    
  }

  /// A template class for the history associated with a population evolving in SMC.

  ///  The history is a template class which should have an associated class type corresponding to
  ///    a _population_ of the desired type, not the type itself.
  ///
  ///    Essentially, this is implemented as a doubly linked list. 


  template <class Population> class history
    {
    private:
      ///The length of the history in time steps
      long  lLength;
      ///The first time step
      historyelement<Population> *pRoot;
      ///The most recent time step
      historyelement<Population> *pLeaf;

    public:
      ///The argument free constructor creates an empty list.
      history();

      ///This function returns a pointer to the first element of the history.
      const historyelement<Population > * GetElement(void) const {return pRoot; }

      /// Returns the effective sample size of the specified particle generation.
      double GetESS(long lGeneration) const;
      ///Returns the number of generations recorded within the history.
      long GetLength (void) const { return lLength; }
      ///Integrate the supplied function over the path of the particle ensemble.
      double IntegratePathSampling(double (*pIntegrand)(long,const Population &,void*), double (*pWidth)(long,void*), void* pAuxiliary);
      double IntegratePathSamplingFinalStep(double (*pIntegrand)(long,const Population&,void*), void* pAuxiliary) const;

      ///Output a vector indicating the number of accepted MCMC moves at each time instance
      void OstreamMCMCRecordToStream(std::ostream &os) const;
      ///Output a 0-1 value vector indicating the times at which resampling occured to an output stream
      void OstreamResamplingRecordToStream(std::ostream &os) const;

      ///Remove the terminal particle generation from the list and return that particle.
      Population* Pop(void);
      ///Remove the terminal particle generation from the list and stick it in the supplied data structures
      void Pop(long* plNumber, Population* ppNew, int* pnAccept, historyflags * phf);
      ///Append the supplied particle generation to the end of the list.
      void Push(long lNumber, const Population & pNew, int nAccept, historyflags hf);


      ///Display the list of particles in a human readable form.
      //  void StreamPopulations(std::ostream & os);
    };

  /// This constructor simply sets the root and leaf pointers to NULL and the list length to zero.
  template <class Population>
  history<Population>::history()
  {
    pRoot = NULL;
    pLeaf = NULL;
    lLength = 0;
  }
   /// Returns the effective sample size of the specified particle generation.
  template <class Population>
  double  history<Population>::GetESS(long lGeneration) const
  {
    historyelement<Population> * pCurrent = pRoot;
    for(long l = 0; l < lGeneration; l++)
      pCurrent = pCurrent->GetNext();
    return pRoot->GetESS();
  }



  /// This function records the MCMC acceptance history to the specified output stream as a list of
  /// the number of moves accepted at each time instant.
  ///
  /// \param os The output stream to send the data to.
  template <class Population>
  void history<Population>:: OstreamMCMCRecordToStream(std::ostream &os) const
  {
    historyelement<Population> * pCurrent = pRoot;

    while(pCurrent) {
      os << pCurrent->AcceptCount() << std::endl;
      pCurrent = pCurrent->GetNext();
    }
  }
  /// This function records the resampling history to the specified output stream as a 0-1 valued list which takes
  /// the value 1 for those time instances when resampling occured and 0 otherwise.
  ///
  /// \param os The output stream to send the data to.
  template <class Population>
  void history<Population>:: OstreamResamplingRecordToStream(std::ostream &os) const
  {
    historyelement<Population> * pCurrent = pRoot;

    while(pCurrent) {
      if(pCurrent->WasResampled())
	os << "1\t";
      else
	os << "0\t";

      os << pCurrent->GetESS() << std::endl;

      pCurrent = pCurrent->GetNext();
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
  /// \param pIntegrand  The function to integrated. The first argument is evolution time, the second a set of particles at which the function is to be evaluated and the final argument is always pAuxiliary.
  /// \param pWidth      The function which returns the width of the path sampling grid at the specified evolution time. The final argument is always pAuxiliary
  /// \param pAuxiliary  A pointer to auxiliary data to pass to both of the above functions

  template <class Population>
  double history<Population>::IntegratePathSampling(double (*pIntegrand)(long,const Population &,void*), double (*pWidth)(long,void*), void* pAuxiliary)
  {
    long lTime = 0;
    long double rValue = 0.0;
    
    historyelement<Population> * pCurrent = pRoot;

    lTime = 1;
    pCurrent=pCurrent->GetNext();
    while(pCurrent) 
      {
	rValue += pCurrent->Integrate(lTime, pIntegrand, pAuxiliary) * (long double)pWidth(lTime,pAuxiliary);
	pCurrent = pCurrent->GetNext();
	lTime++; 
      }
    return ((double)rValue);
  }

  template <class Population>
  double history<Population>::IntegratePathSamplingFinalStep(double (*pIntegrand)(long,const Population&,void*), void* pAuxiliary) const
  {
    return pLeaf->Integrate(lLength-1,pIntegrand,pAuxiliary);
  }


  /// Pop() operates just as the standard stack operation does. It removes the particle generation currently occupying
  /// the terminal position in the chain, decrements the length counter and returns the particle set as an array.
  template <class Population>
  Population* history<Population>::Pop(void)
  {
    if(lLength == 0)
      return NULL;

    Population * rValue = pLeaf->GetPointers();

    lLength--;

    if(lLength == 0)
      pRoot = pLeaf = 0;
    else {
      pLeaf = pLeaf->GetLast();
      delete pLeaf->GetNext();
      pLeaf->SetNext(NULL);

    }
    return rValue;
  }

  /// Pop operates as the usual stack operation
  ///
  /// If called with four pointers of this sort, it removes the last particle generation from the history stack and
  /// places them in the reference objects.
  template <class Population>
  void history<Population>::Pop(long* plNumber, Population* ppNew, int* pnAccept, historyflags * phf)
  {
    if(plNumber)
      (*plNumber) = pLeaf->GetNumber();
    if(ppNew) {
      //for(long l = 0; l < *plNumber; l++)
      (*ppNew)    = pLeaf->GetValues();
    }
    if(pnAccept)
      (*pnAccept) = pLeaf->AcceptCount();
    if(phf)
      (*phf)      = pLeaf->GetFlags();

    if(lLength <= 1) {
      pRoot = NULL;
      pLeaf = NULL;
    }
    else {
      pLeaf = pLeaf->GetLast();
      //      delete pLeaf->GetNext();
      pLeaf->SetNext(NULL);

    }

    lLength--;

    return;

  }

  /// Push operates just like the standard stack operation: it adds the specified particle set generation to the history
  /// of the particle set and increments the stack counter.
  ///
  /// \param lNumber The number of particles present in this generation of the system.
  /// \param pNew    An array containing the particles present in this generation of the system.
  /// \param nAccepts The number of accepted MCMC moves during this iteration of the system
  /// \param hf      The historyflags associated with this generation of the system.

  template <class Population>
  void history<Population>::Push(long lNumber, const Population & pNew, int nAccepts, historyflags hf)
  {
    if(lLength == 0) {
      pRoot = new historyelement<Population>(lNumber, pNew, nAccepts, hf);
      pLeaf = pRoot;
    }
    else {
      pLeaf = new historyelement<Population>(lNumber, pNew, pLeaf, NULL, nAccepts, hf);
      pLeaf->GetLast()->SetNext(pLeaf);
    } 
    lLength++;
  }
}



namespace std {
  /// This function will ultimately allow the standard stream operators to be used to display a particle history in a human readable form.

  /// It isn't yet fully implemented.
  template <class Population>
  ostream & operator<<(ostream & os, smc::history<Population> h)
  {
    h.StreamPopulations(os);
    return os;
  }
  
  
  
  
  
  
  
  
  
}
#endif
