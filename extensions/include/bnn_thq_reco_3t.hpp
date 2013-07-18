/*******************************************************************************
 * This file was automatically generated by bnn-hep package
 * https://github.com/andrey-popov/bnn-hep
 * 
 * Thu Jul 18 02:19:29 2013
 * 
 * 
 * The list of input variables:
 *  var #0: abs(EtaHiggs)
 *  var #1: abs(EtaLJet)
 *  var #2: DeltaEtaLepTop
 *  var #3: DeltaRBJetsHiggs
 *  var #4: log(MassHiggs)
 *  var #5: log(MassTop)
 *  var #6: log(MinPtBJet)
 *  var #7: NPassBTagHiggs
 *  var #8: PassBTagTop
 *  var #9: RelHt
 * 
 * From the code below, the user is expected to only use the class BNN.
 ******************************************************************************/

#ifndef THQ_RECO_3T_HPP
#define THQ_RECO_3T_HPP

#include <Rtypes.h>
#include <TMath.h>

#include <algorithm>
#include <list>
#include <string>


#ifdef USE_INTERNAL_BINARY_DISCRIMINATOR_ABSTRACT_BASE
#ifndef BINARY_DISCRIMINATOR_ABSTRACT_BASE
#define BINARY_DISCRIMINATOR_ABSTRACT_BASE

class BinaryDiscriminator
{
	public:
		virtual Double_t operator()(Double_t const *) const = 0;
};
#endif

#else
#include "BinaryDiscriminator.hpp"
#endif


namespace thq_reco_3t{

extern std::string const taskName;
extern std::list<std::string> inputVarNames;
extern std::list<std::string> sgnFileNames;
extern std::list<std::string> bkgFileNames;

void Initialize();

class Transform0
{
	public:
		Transform0();
		void operator()(Double_t *vars) const;

	private:
		Double_t mean[10], sigma[10];
};


class NN
{
	public:
		NN();
	
	public:
		void SetWeightsL1(Double_t const [15][10]);
		void SetBiasesL1(Double_t const [15]);
		void SetWeightsL2(Double_t const [1][15]);
		void SetBiasesL2(Double_t const [1]);
		Double_t const * Apply(Double_t const *) const;
	
	private:
		Double_t weightsL1[15][10];
		Double_t biasesL1[15];
		Double_t weightsL2[1][15];
		Double_t biasesL2[1];
		mutable Double_t bufferIn[15];
		mutable Double_t bufferOut[15];
};


class BNN: public BinaryDiscriminator
{	 public:
		BNN(UInt_t netBegin_ = 0, UInt_t netEnd_ = 50);
	
	public:
		void SetNetRange(UInt_t netBegin_, UInt_t netEnd_);
		Double_t operator()(Double_t const *vars) const;
		Double_t operator()(Double_t var0, Double_t var1, Double_t var2, Double_t var3, Double_t var4, Double_t var5, Double_t var6, Double_t var7, Double_t var8, Double_t var9) const;
	
	private:
		Double_t Apply(Double_t const *vars) const;
	
	private:
		NN nets[50];
		UInt_t netBegin, netEnd;
		Transform0 trans0;
};

}

#endif
