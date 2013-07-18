/**
 * \file THInputVarsPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a class that creates tuples with observables for disrimination of thq from its
 * backgrounds.
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


/**
 * \class THInputVarsPlugin
 * \brief Creates tuples to train a BNN for thq MVA reconstruction
 * 
 * There are four energetic jets in the final state of thq process. The class considers all the ways
 * to choose four jets out of all jets in an event (combination of PECReader::GetJets() and
 * PECReader::GetAdditionalJets()). All the possible associations of these four jets to parent
 * objects (top-quark, Higgs boson, or recoil quark) are constructed. Such combinatorics defines a
 * set of interpretations of an event. For each interpretations a distance between reconstructed and
 * generator-level top-quark, Higgs boson, and recoil quark is calculated (metrics $\Delta R \oplus
 * \Delta p_T^{rel}$ is used). The interpretations are ordered in the distance (in increasing
 * order). Finally, a set of reconstruction-level observables is calculated for each interpretation
 * and stored in a ROOT tree (along with information about the distance). The user can also
 * configure the plugin to store only one interpretation per event (see description of
 * THInputVarsPlugin::pruned data member for details).
 * 
 * Details are available a talk by Andrey Popov in [1].
 * [1] https://indico.cern.ch/conferenceDisplay.py?confId=251808 (password is "tHmeeting")
 */
class THInputVarsPlugin: public Plugin
{
    public:
        /// Constructor
        THInputVarsPlugin(std::string const &outDirectory, BTagger const &bTagger);

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
        
        
        // Observables defined under thq hypothesis
        Float_t thq_MassHiggs, thq_PtHiggs, thq_EtaHiggs;
        Float_t thq_PtLJet, thq_EtaLJet;
        
        Float_t thq_DeltaRTopHiggs, thq_DeltaRBJetsHiggs;
        
        Float_t thq_CosLepLJetTH;
        //^ Angle between three-momenta of the lepton and the light-flavour jet in the rest frame of
        //t+h system
        
        Float_t thq_MassTopHiggs;
        
        Float_t tt_MassTopHad, tt_PtTopHad, tt_EtaTopHad;
        Float_t tt_MassWHad, tt_PtWHad, tt_EtaWHad;
        
        Float_t tt_RelHt;
        //^ Sum of pt of the two top-quarks devided by the total Ht of the event
        
        Float_t tt_DeltaRLightJets;
        
        Float_t tt_MaxMassBHadQ;
        //^ max m(b + q), where b is the b-jet from hadronically decaying top-quark and q is one of
        //light-flavour jets from the subsequent W-boson decay
        
        Float_t glb_PtJ1, glb_PtJ2;
        
        Float_t glb_SqrtSHat;
        //^ Calculated as mass of the sum of all objects
        
        Float_t glb_Sphericity;
        
        
        Float_t weight;
};