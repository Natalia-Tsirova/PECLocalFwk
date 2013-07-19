/**
 * \file THEvalBNNPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a class that applies a BNN of thq extraction and stores its decision in ROOT tuples.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <THRecoPlugin.hpp>
#include <TTbarRecoPlugin.hpp>
#include <BTagger.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <string>
#include <vector>


#include <bnn_ttbar_discr_3t.hpp>


/**
 * \class THEvalBNNPlugin
 * \brief Creates tuples to train a BNN for extraction of thq events
 * 
 * The class calculates and stores a set of input variables to perform tqh vs ttbar discrimination.
 * Each event is reconstructed under two alternative hypotheses in parallel: as a thq event and as
 * a semileptonic ttbar one. (The reconstruction is done with the help of dedicated plugins.) For
 * each hypothesis a set of observables is calculated and stored in a ROOT file. In addition to it,
 * several global variables that do not rely on event interpretation are calculated.
 */
class THEvalBNNPlugin: public Plugin
{
    public:
        /// Constructor
        THEvalBNNPlugin(std::string const &outDirectory, BTagger const &bTagger);

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
        
        /// Pointer to a plugin for thq reconstruction
        THRecoPlugin const *thqReconstructor;
        
        /// Pointer to a plugin for ttbar reconstrucion
        TTbarRecoPlugin const *ttbarReconstructor;
        
        /// An object to perform b-tagging
        BTagger const &bTagger;
        
        /// Directory to store output files
        std::string outDirectory;
        
        /// The BNN to discriminate thq from ttbar
        ttbar_discr_3t::BNN bnnDiscr;
        
        
        /// Vector off all jets in an event
        std::vector<Jet> allJets;
        
        /// Current output file
        TFile *file;
        
        /// Current output tree
        TTree *tree;
        
        
        // Output buffers
        ULong64_t eventNumber, runNumber, lumiSection;
        
        Float_t NJets30, NTags30;
        //^ Needed to allow the same tuples for 3t and 4t bins
        
        Float_t bnnDecision;
        
        Float_t weight;
};