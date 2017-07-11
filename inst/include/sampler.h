// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// sampler.h: Rcpp integration of SMC library -- sampler object 
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

//! \file
//! \brief Defines the overall sampler object.
//!
//! This file defines the smc::sampler class which is used to implement entire particle systems.

#ifndef __SMC_SAMPLER_HH

#define __SMC_SAMPLER_HH 1.0

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>

#include "history.h"
#include "moveset.h"
#include "algParam.h"
#include "population.h"
#include "smc-exception.h"

///Specifiers for various path sampling methods:
namespace PathSamplingType
{
	enum Enum { TRAPEZOID2 = 0, 
		TRAPEZOID1, 
		RECTANGLE};
};


///Storage types for the history of the particle system. 
namespace HistoryType
{
	enum Enum {NONE = 0, 
		RAM};
};

namespace smc {

	/// A template class for an interacting particle system suitable for SMC sampling
	template <class Space> 
	class sampler
	{
	protected:
		///Number of particles in the system.
		long N;
		///The current evolution time of the system.
		long T;
		
		///Structure used internally for resampling.
		arma::vec dRSWeights;
		///Structure used internally for resampling.
		arma::Col<unsigned int> uRSCount;
		///Structure used internally for resampling.
		arma::Col<unsigned int> uRSIndices;

		///The particles within the system.
		population<Space> pPopulation;
		///The set of moves available.
		moveset<Space> Moves;
		///The algorithm parameters.
		algParam<Space>* pAlgParams;
		///The historical process associated with the particle system.
		std::vector<historyelement< population<Space> > > History;

		///The number of MCMC moves which have been accepted during this iteration
		int nAccepted;
		///A flag which tracks whether the ensemble was resampled during this iteration
		int nResampled;

		///A mode flag which indicates whether historical information is stored
		HistoryType::Enum htHistoryMode;
		///An estimate of the overall ratio of normalising constants
		double dlogNCPath;
		///An estimate of the latest iteration's ratio of normalising constants
		double dlogNCIt;

	public:
		///Create a particle system containing lSize uninitialised particles with the specified mode.
		sampler(long lSize, HistoryType::Enum htHistoryMode);
		///Dispose of a sampler.
		~sampler();
		///Calculates and Returns the Effective Sample Size.
		double GetESS(void) const;
		/// Returns the effective sample size of the specified particle generation.
		double GetESS(long lGeneration) const;
		///Returns a pointer to the History of the particle system
		const std::vector<historyelement< population<Space> > > & GetHistory(void) const { return History; }
		///Returns a pointer to the algorithm parameters
		algParam<Space> * GetAlgs(void) const { return pAlgParams; }																											   
		///Returns the current estimate of the log normalising constant ratio over the entire path
		double GetLogNCPath(void) const { return dlogNCPath; }
		///Returns the current estimate of the log normalising constant ratio over the last step
		double GetLogNCStep(void) const { return dlogNCIt; }
		///Returns the current estimate of the normalising constant ratio over the entire path
		double GetNCPath(void) const { return exp(dlogNCPath); }
		///Returns the current estimate of the normalising constant ratio over the last step
		double GetNCStep(void) const { return exp(dlogNCIt); }
		///Returns the number of particles within the system.
		long GetNumber(void) const {return N;}
		///Returns the number of particles within the system.
		long GetHistoryLength(void) const {return History.size();}
		///Return the values of particles
		const std::vector<Space> &  GetPopulationValue(void) const { return pPopulation.GetValue(); }
		///Return the values of particles
		const population<Space> &  GetPopulation(void) const { return pPopulation; }
		///Return the values of particles
		const Space &  GetPopulationValueN(long n) const { return pPopulation.GetValueN(n); }
		///Return the logarithmic unnormalized weights of particles
		const arma::vec & GetPopulationLogWeight(void) const { return pPopulation.GetLogWeight(); }
		///Return the unnormalized weights of particls
		arma::vec GetPopulationWeight(void) const { return pPopulation.GetWeight(); }  
		///Return the unnormalized weights of particls
		double GetPopulationWeightN(int n) const { return pPopulation.GetWeightN(n); }  
		///Returns the current evolution time of the system.
		long GetTime(void) const {return T;}
		///Initialise the sampler and its constituent particles.
		void Initialise(void);
		///Integrate the supplied function with respect to the current particle set.
		double Integrate(double(*pIntegrand)(const Space &,void*), void* pAuxiliary) const;
		///Integrate the supplied function over the path path using the supplied width function and integration method.
		double IntegratePathSampling(PathSamplingType::Enum, double (*pIntegrand)(long,const population<Space>&,long,void*), double (*pWidth)(long,void*), void* pAuxiliary);
		///Integrate the supplied function over the path path using the supplied width function and the default integration method (the corrected trapezoid rule).
		double IntegratePathSampling(double (*pIntegrand)(long,const population<Space>&,long,void*), double (*pWidth)(long,void*), void* pAuxiliary) {return IntegratePathSampling(PathSamplingType::TRAPEZOID2, pIntegrand, pWidth, pAuxiliary);}
		///Perform one iteration of the simulation algorithm.
		void Iterate(void);
		///Cancel one iteration of the simulation algorithm.
		void IterateBack(void);
		///Perform one iteration of the simulation algorithm and return the resulting ess
		double IterateEss(void);
		///Perform iterations until the specified evolution time is reached
		void IterateUntil(long lTerminate);
		///Move the particle set by proposing an applying an appropriate move to each particle.
		void MovePopulations(void);
		///Resample the particle set using the specified resmpling scheme.
		void Resample(ResampleType::Enum lMode);
		///Sets the entire moveset to the one which is supplied
		void SetMoveSet(moveset<Space>& pNewMoveset) {Moves = pNewMoveset;}
		///Set Resampling Parameters
		void SetSMCParams(ResampleType::Enum rtMode, double dThreshold);
		///Set Resampling Parameters
		void SetSMCParams(algParam<Space>* inParams);
		///Dump a specified particle to the specified output stream in a human readable form
		const std::ostream & StreamParticle(std::ostream & os, long n) const;
		///Dump the entire particle set to the specified output stream in a human readable form
		const std::ostream & StreamParticles(std::ostream & os) const;
		///Output a vector indicating the number of accepted MCMC moves at each time instance
		void OstreamMCMCRecordToStream(std::ostream &os) const;
		///Output a 0-1 value vector indicating the times at which resampling occured to an output stream
		void OstreamResamplingRecordToStream(std::ostream &os) const;
		

	private:
		///Duplication of smc::sampler is not currently permitted.
		sampler(const sampler<Space> & sFrom);
		///Duplication of smc::sampler is not currently permitted.
		sampler<Space> & operator=(const sampler<Space> & sFrom);
		///Generate a multinomial random vector with parameters (n,w[1:k]) and store it in X
		void Multinomial(unsigned n, unsigned k, arma::vec w, unsigned int * X);		
		

protected:
	///Returns the crude normalising constant ratio estimate implied by the weights.
	double CalcLogNC(void) const;
};

template <class Space>
void sampler<Space>::Multinomial(unsigned n, unsigned k, arma::vec w, unsigned int * X) {		
	Rcpp::IntegerVector v(k);
	w = w/arma::sum(w);
	
	double * w_mem = w.memptr();
	
	// R sources:  rmultinom(int n, double* prob, int K, int* rN);
	rmultinom(static_cast<int>(n), const_cast<double*>(w_mem), static_cast<int>(k), &(v[0]));
	
	for (unsigned int i=0; i<k; i++) {
		X[i] = static_cast<unsigned int>(v[i]);
	}	
	}


	/// The constructor prepares a sampler for use but does not assign any moves to the moveset, initialise the particles
	/// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
	/// store the particle set and to initialise a random number generator.
	///
	/// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
	/// \param htHM The history mode to use: set this to HistoryType::RAM to store the whole history of the system and SMC_HISTORY_NONE to avoid doing so.
	/// \tparam Space The class used to represent a point in the sample space.
	template <class Space>
	sampler<Space>::sampler(long lSize, HistoryType::Enum htHM)
	{N = lSize;
	pAlgParams = new algParam<Space>(lSize);
		uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
		htHistoryMode = htHM;
	}

	template <class Space>
	sampler<Space>::~sampler()
	{
		delete pAlgParams;
	}


	template <class Space>
	double sampler<Space>::GetESS(void) const
	{
		double sum = arma::sum(exp(pPopulation.GetLogWeight()));
		double sumsq = arma::sum(exp(2.0*pPopulation.GetLogWeight()));

		return expl(-log(sumsq) + 2*log(sum));
	}



	/// Returns the effective sample size of the specified particle generation.
	template <class Space>
	double  sampler<Space>::GetESS(long lGeneration) const
	{
		typename std::vector<historyelement<population<Space> > >::const_iterator it = History.begin();
		std::advance(it,lGeneration);
		return it->GetESS(); 
	}


	template <class Space>
	double sampler<Space>::CalcLogNC(void) const
	{
		double dMaxWeight = arma::max(pPopulation.GetLogWeight());
		double sum = arma::sum(exp(pPopulation.GetLogWeight() - dMaxWeight*arma::ones(N)));

		return (dMaxWeight + log(sum));
	}


	/// The initialise function:
	///          -# resets the system evolution time to 0 and calls the moveset initialisor to assign each particle in the ensemble.
	///         -# checks the effective sample size and resamples if necessary
	///         -# performs a mcmc step if required
	///         -# appends the particle set to the history if desired
	///
	/// Note that the initialisation function must be specified before calling this function.

	template <class Space>
	void sampler<Space>::Initialise(void)
	{
		T = 0;
		dlogNCIt = 0;
		dlogNCPath = 0;
		nAccepted = 0;

		std::vector<Space> InitVal(N);
		arma::vec InitWeights(N);
		pPopulation = population<Space>(InitVal,InitWeights);
		Moves.DoInit(pPopulation,N);

		//Scaling weights by 1/N (mostly for evidence computation)
		pPopulation.SetLogWeight(pPopulation.GetLogWeight() - log(static_cast<double>(N))*arma::ones(N));

		//Estimate the normalising constant
		dlogNCIt = CalcLogNC();
		dlogNCPath += dlogNCIt;

		//Normalise the weights to sensible values....
		pPopulation.SetLogWeight(pPopulation.GetLogWeight() - dlogNCIt*arma::ones(N));
		
		pAlgParams->updateForMCMC(pPopulation);
		//Check if the ESS is below some reasonable threshold and resample if necessary.
		double ESS = GetESS();
		if(ESS < pAlgParams->GetResThresh()) {
			nResampled = 1;
			Resample(pAlgParams->GetResample());
		}
		else {
			nResampled = 0;
		}
		//A possible MCMC step could be included here.
		nAccepted += Moves.DoMCMC(0,pPopulation, N); 

		if(htHistoryMode != HistoryType::NONE) {
			History.clear();
			historyelement<population<Space> > histel;
			histel.Set(N, pPopulation, nAccepted, historyflags(nResampled));
			History.push_back(histel);
			//History.emplace_back(historyelement<population<Space> >(N, pPopulation, nAccepted, historyflags(nResampled)));
			
		}  
		
		pAlgParams->updateEnd(pPopulation);
		return;
	}

	/// This function returns the result of integrating the supplied function under the empirical measure associated with the
	/// particle set at the present time. The final argument of the integrand function is a pointer which will be supplied
	/// with pAuxiliary to allow for arbitrary additional information to be passed to the function being integrated.
	///
	/// \param pIntegrand The function to integrate with respect to the particle set 
	/// \param pAuxiliary A pointer to any auxiliary data which should be passed to the function

	template <class Space>
	double sampler<Space>::Integrate(double(*pIntegrand)(const Space&,void*), void * pAuxiliary) const
	{
		long double rValue = 0;
		for(int i =0; i < N; i++)
		{
			rValue += expl(pPopulation.GetLogWeightN(i)) * pIntegrand(pPopulation.GetValueN(i), pAuxiliary);
		}

		rValue /= expl(CalcLogNC());
		return (double)rValue;
	}


	/// This function is intended to be used to estimate integrals of the sort which must be evaluated to determine the
	/// normalising constant of a distribution obtain using a sequence of potential functions proportional to densities with respect
	/// to the initial distribution to define a sequence of distributions leading up to the terminal, interesting distribution.
	///
	/// In this context, the particle set at each time is used to make an estimate of the path sampling integrand, and a
	/// trapezoidal integration is then performed to obtain an estimate of the path sampling integral which is the natural logarithm
	/// of the ratio of normalising densities.

	/// The function performs a trapezoidal integration of the type which is useful when using path sampling to estimate the
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

	template <class Space>
	double sampler<Space>::IntegratePathSampling(PathSamplingType::Enum PStype, double (*pIntegrand)(long,const population<Space> &,long,void*), double (*pWidth)(long,void*), void* pAuxiliary)
	{
		if(htHistoryMode == HistoryType::NONE)
		throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "The path sampling integral cannot be computed as the history of the system was not stored.");

		// historyelement<population<Space> > histel;
		// histel.Set(N, pPopulation, nAccepted, historyflags(nResampled));
		// History.push_back(histel);
		
		long lTime = 1;
		long double rValue = 0.0;
		typename std::vector<historyelement<population<Space> > >::const_iterator it;	
		
		switch(PStype) {
		case PathSamplingType::RECTANGLE:
			{
				for(it = ++History.begin(); it!=History.end(); it++){
					rValue += it->Integrate(lTime, pIntegrand, pAuxiliary) * (long double)pWidth(lTime,pAuxiliary);
					lTime++;
				}
				break;
			}
			
			
		case PathSamplingType::TRAPEZOID1:
			{	

				long double previous_expt = History.begin()->Integrate(0,pIntegrand,pAuxiliary);
				long double current_expt;
				for(it = ++History.begin(); it!=History.end(); it++){
					current_expt = it->Integrate(lTime, pIntegrand, pAuxiliary);
					rValue += (long double)pWidth(lTime,pAuxiliary)/2.0 * (previous_expt + current_expt) ;
					lTime++;
					previous_expt = current_expt;
				}
				
				break;
			}

		case PathSamplingType::TRAPEZOID2:
		default:
			{
				long double previous_expt = History.begin()->Integrate(0,pIntegrand,pAuxiliary);
				long double previous_var = History.begin()->Integrate_Var(0,pIntegrand,previous_expt,pAuxiliary);
				long double current_expt;
				long double current_var;
				long double width = 0.0;
				for(it = ++History.begin(); it!=History.end(); it++){
					current_expt = it->Integrate(lTime, pIntegrand, pAuxiliary);
					current_var = it->Integrate_Var(lTime, pIntegrand, current_expt, pAuxiliary);
					width = (long double)pWidth(lTime,pAuxiliary);
					rValue += width/2.0 * (previous_expt + current_expt) - pow(width,2)/12.0*(current_var - previous_var);
					lTime++;
					previous_expt = current_expt;
					previous_var = current_var;
				}	
				
				break;
			}

		}
		
		// History.pop_back();
		
		return ((double)rValue);
	}

	/// The iterate function:
	///          -# moves the current particle set
	///         -# checks the effective sample size and resamples if necessary
	///         -# performs a mcmc step if required
	///         -# appends the current particle set to the history if desired
	///         -# increments the current evolution time
	template <class Space>
	void sampler<Space>::Iterate(void)
	{
		IterateEss();
		return;
	}

	template <class Space>
	void sampler<Space>::IterateBack(void)
	{
		if(htHistoryMode == HistoryType::NONE)
		throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "An attempt to undo an iteration was made; unforunately, the system history has not been stored.");

		History.pop_back();
		historyelement<population<Space> > recent = History.back();
		pPopulation = recent.GetRefs();
		N =recent.GetNumber();
		nAccepted = recent.AcceptCount();
		nResampled = recent.WasResampled();
		T--;
		return;
	}

	template <class Space>
	double sampler<Space>::IterateEss(void)
	{
		nAccepted = 0;
		pAlgParams->updateForMove(pPopulation);
		//Move the particle set.
		MovePopulations();

		//Estimate the normalising constant
		dlogNCIt = CalcLogNC();
		dlogNCPath += dlogNCIt;
		
		//Normalise the weights
		pPopulation.SetLogWeight(pPopulation.GetLogWeight()  - dlogNCIt*arma::ones(N));

		pAlgParams->updateForMCMC(pPopulation);
		//Check if the ESS is below some reasonable threshold and resample if necessary.
		//A mechanism for setting this threshold is required.
		double ESS = GetESS();
		if(ESS < pAlgParams->GetResThresh()) {
			nResampled = 1;
			Resample(pAlgParams->GetResample());
		}
		else
		nResampled = 0;
		//A possible MCMC step could be included here.
		nAccepted += Moves.DoMCMC(T+1,pPopulation,N);
		// Increment the evolution time.
		T++;

		//Finally, the current particle set should be appended to the historical process.
		if(htHistoryMode != HistoryType::NONE){
			historyelement<population<Space> > histel;
			histel.Set(N, pPopulation, nAccepted, historyflags(nResampled));
			History.push_back(histel);
			//History.emplace_back(historyelement<population<Space> >(N, pPopulation, nAccepted, historyflags(nResampled)));
		}
		
		pAlgParams->updateEnd(pPopulation);

		return ESS;
	}

	template <class Space>
	void sampler<Space>::IterateUntil(long lTerminate)
	{
		while(T < lTerminate)
		Iterate();
	}

	template <class Space>
	void sampler<Space>::MovePopulations(void)
	{
		Moves.DoMove(T+1,pPopulation, N);
	}

	template <class Space>
	void sampler<Space>::Resample(ResampleType::Enum lMode)
	{
		//Resampling is done in place.
		unsigned uMultinomialCount;

		
		//First obtain a count of the number of children each particle has.
		switch(lMode) {
		case ResampleType::MULTINOMIAL:
			//Sample from a suitable multinomial vector
			dRSWeights = pPopulation.GetWeight();
			//uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
			Multinomial(N,N,dRSWeights,uRSCount.memptr());
			//uRSCount = pRng->Multinomial(N,N,dRSWeights);
			break;
			
			
		case ResampleType::RESIDUAL:
			dRSWeights = exp(log((double)N)*arma::ones(N) + pPopulation.GetLogWeight() - CalcLogNC()*arma::ones(N));
			
			uRSIndices = arma::zeros<arma::Col<unsigned int> >((int)N);
			
			
			for(int i = 0; i < N; ++i)
			uRSIndices(i) = unsigned(floor(dRSWeights(i)));
			dRSWeights = dRSWeights - uRSIndices;
			uMultinomialCount = N - arma::sum(uRSIndices);
			
			//uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
			Multinomial(uMultinomialCount,N,dRSWeights,uRSCount.memptr());
			uRSCount += uRSIndices;
			
			break;
			
		case ResampleType::STRATIFIED:
		default:
			{
				// Procedure for stratified sampling
				//Generate a random number between 0 and 1/N times the sum of the weights
				double dRand = R::runif(0,1/(double)N);
				
				arma::vec dWeightCumulative = arma::cumsum(exp(pPopulation.GetLogWeight() - CalcLogNC()));
				
				int k = 0;
				int j = 0;
				uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
				while(j < N) {
					while((dWeightCumulative(k) - dRand) > ((double)j)/((double)N) && j < N) {
						uRSCount(k)++;
						j++;
						dRand = R::runif(0,1/(double)N);
					}
					k++;
				}
				break;
			}

		case ResampleType::SYSTEMATIC:
			{
				// Procedure for stratified sampling but with a common RV for each stratum
				//Generate a random number between 0 and 1/N times the sum of the weights
				double dRand = R::runif(0,1/(double)N);
				
				int j = 0, k = 0;
				uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
				arma::vec dWeightCumulative = arma::cumsum(exp(pPopulation.GetLogWeight() - CalcLogNC()*arma::ones(N)));
				while(k < N) {
					while((dWeightCumulative(k) - dRand) > ((double)j)/((double)N) && j < N) {
						uRSCount(k)++;
						j++;
					}
					k++;
				}
				break;
			}

		}

		uRSIndices = arma::zeros<arma::Col<unsigned int> >((int)N);
		//Map count to indices to allow in-place resampling
		for (int i=0, j=0; i<N; ++i) {
			if (uRSCount(i)>0) {
				uRSIndices(i) = i;
				while (uRSCount(i)>1) {
					while (uRSCount(j)>0) ++j; // find next free spot
					uRSIndices(j++) = i; // assign index
					--uRSCount(i); // decrement number of remaining offsprings
				}
			}
		}
		Space Current;
		//Perform the replication of the chosen.
		for(int i = 0; i < N ; ++i) {
			if(uRSIndices(i) != static_cast<unsigned int>(i)){
				Current = pPopulation.GetValueN((int)uRSIndices(i));
				pPopulation.SetValueN(Current,i);
			}
		}
		//Set equal normalised weights
		pPopulation.SetLogWeight(- log(static_cast<double>(N))*arma::ones(N));
	} 				
	
	/// This function configures the resampling parameters, allowing the specification of both the resampling
	/// mode and the threshold at which resampling is used.
	///
	/// \param rtMode The resampling mode to be used.
	/// \param dThreshold The threshold at which resampling is deemed necesary.
	///
	/// The rtMode parameter should be set to one of the following:
	/// -# ResampleType::MULTINOMIAL to use multinomial resampling  
	/// -# ResampleType::RESIDUAL to use residual resampling
	/// -# ResampleType::STRATIFIED to use stratified resampling
	/// -# ResampleType::SYSTEMATIC to use systematic resampling
	///
	/// The dThreshold parameter can be set to a value in the range [0,1) corresponding to a fraction of the size of
	/// the particle set or it may be set to an integer corresponding to an actual effective sample size.

	template <class Space>
	void sampler<Space>::SetSMCParams(algParam<Space>* inParams)
	{
		pAlgParams = inParams; //later want to add lSize, history type to this
	}
	
	template <class Space>
	void sampler<Space>::SetSMCParams(ResampleType::Enum rtMode, double dThreshold)
	{
		pAlgParams->SetResample(rtMode);
		pAlgParams->SetResThresh(dThreshold);
	}

	template <class Space>
	const std::ostream & sampler<Space>::StreamParticle(std::ostream & os, long n) const
	{
		Space val = pPopulation.GetValueN(n);
		double weight = pPopulation.GetWeightN(n);
		os << val << "," << weight;
		return os;
	}

	template <class Space>
	const std::ostream & sampler<Space>::StreamParticles(std::ostream & os) const
	{
		Space value;
		double weight;
		for(int i = 0; i < N - 1; i++){
			value = pPopulation.GetValueN(i);
			weight = pPopulation.GetWeightN(i);
			os << value << "," << weight << std::endl;
		}
		
		return os;
	}

	/// This function records the MCMC acceptance history to the specified output stream as a list of
	/// the number of moves accepted at each time instant.
	///
	/// \param os The output stream to send the data to.
	template <class Space>
	void sampler<Space>:: OstreamMCMCRecordToStream(std::ostream &os) const
	{
		os << "Accepted MCMC proposals history:" << std::endl;
		os << "======================" << std::endl;
		for(typename std::vector<historyelement<population<Space> > >::const_iterator it = History.begin(); it!=History.end(); it++){
			os << it->AcceptCount() << std::endl;
		}
	}
	/// This function records the resampling history to the specified output stream as a 0-1 valued list which takes
	/// the value 1 for those time instances when resampling occured and 0 otherwise.
	///
	/// \param os The output stream to send the data to.
	template <class Space>
	void sampler<Space>:: OstreamResamplingRecordToStream(std::ostream &os) const
	{
		os << "Resampling history:" << std::endl;
		os << "======================" << std::endl;
		os << "Flag\t" << "ESS\t" << std::endl;
		for(typename std::vector<historyelement<population<Space> > >::const_iterator it = History.begin(); it!=History.end(); it++){ 
			if(it->WasResampled())
			os << "1\t";
			else
			os << "0\t";

			os << it->GetESS() << std::endl;
		}
	}


}


namespace std {
	/// Produce a human-readable display of the state of an smc::sampler class using the stream operator.

	/// \param os The output stream to which the display should be made.
	/// \param s  The sampler which is to be displayed.
	template <class Space>
	std::ostream & operator<< (std::ostream & os, smc::sampler<Space> & s)
	{
		os << "Sampler Configuration:" << std::endl;
		os << "======================" << std::endl;
		os << "Evolution Time:   " << s.GetTime() << std::endl;
		os << "Particle Set Size: " << s.GetNumber() << std::endl;
		os << "Effective Sample Size: " << s.GetESS() << std::endl;
		os << std::endl;
		os << "Particle Set: " << std::endl;
		s.StreamParticles(os);
		os << std::endl;
		return os;
	}
}


#endif
