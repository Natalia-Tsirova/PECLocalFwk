/**
 * \file SingleTopTChanPlugin_wjets.hpp
 * \author Andrey Popov, Natalia Tsirova
 * 
 * The module defines a class to calculate a set of variables that are used in single-top t-channel
 * analysis.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <WjetsHFPlugin.hpp>
#include <BTagger.hpp>
//#include <SystDefinition.hpp>

#include <TFile.h>
#include <TTree.h>

#include <string>
#include <memory>


/**
 * \class SingleTopTChanPlugin_wjets
 * \brief Calculates and stores variables that are used in single-top t-channel analysis
 * 
 * The class is expected to serve as an illustration rather than to be used in a real-life analysis.
 */
class SingleTopTChanPlugin_wjets: public Plugin
{
    public:
        /// Constructor
        SingleTopTChanPlugin_wjets(std::string const &outDirectory, std::shared_ptr<BTagger const> &bTagger, bool const isWeightSyst, std::string const hftype_str);
    
    public:
        /**
         * \brief Creates a newly-initialized copy
         * 
         * Consult documentation of the overriden method for details.
         */
        Plugin *Clone() const;
        
        /**
         * \brief Notifies this that a dataset has been opened
         * 
         * Consult documentation of the overriden method for details.
         */
        void BeginRun(Dataset const &dataset);
        
        /**
         * \brief Notifies this that a dataset has been closed
         * 
         * Consult documentation of the overriden method for details.
         */
        void EndRun();
        
        /**
         * \brief Processes the current event
         * 
         * Consult documentation of the overriden method for details.
         */
        bool ProcessEvent();
    
    private:
        /// Pointer to PECReaderPlugin
        PECReaderPlugin const *reader;
        
        /// Pointer to a plugin for ttbar classification
        WjetsHFPlugin const *WjetsHFClassifier;
        
        /// An object to perform b-tagging
        std::shared_ptr<BTagger const> bTagger;
        
        /// Directory to store output files
        std::string outDirectory;



	/// Systematics
	//SystTypeAlgo const &syst;
	bool const isWeightSyst;
	
	/// heavy flavour type
        WjetsHFPlugin::Type hftype;
        std::string const hftype_str;
        
        /// Current output file
        TFile *file;
        
        /// Current output tree
        TTree *tree;
        
        // Output buffers
        ULong64_t eventNumber, runNumber, lumiSection;
        
        Float_t Pt_Lep, Eta_Lep, RelIso_Lep;
        Float_t MET, MtW;
        Float_t Phi_MET, DPhi_LepNu;
        
        Float_t Pt_J1, Eta_J1, Pt_J2, Eta_J2;
        Float_t Pt_LJ, Eta_LJ;
        Float_t Pt_BJ1, Pt_BJ2;
        
        Float_t M_J1J2, DR_J1J2, Pt_J1J2;
        
        Float_t DR_LepJ1, DR_LepJ2, DPhi_LepJ1;
        
        Int_t N_J, N_BJ, N_LJ, Charge_Lep;
        
        Float_t Ht, Ht_J, Ht_JNotBest, M_J, M_JNotBest, Pt_JNotBest;
        Float_t M_JW;
        
        Float_t Mtop_BJ1, Mtop_BestJ, Pttop_BJ1;
        Float_t Cos_LepLJ_BJ1, Cos_WLJ_BJ1;
        
        Float_t Ht_J1J2, Pt_W, Cos_LepJ1;
        
        Float_t Sphericity, Planarity, Aplanarity;
        
        Int_t nPV, WHFClass;
        Float_t weight, weight_PileUpUp, weight_PileUpDown;
	Float_t weight_TagRateUp, weight_TagRateDown, weight_MistagRateUp, weight_MistagRateDown;
};
