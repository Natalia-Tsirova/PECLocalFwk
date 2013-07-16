/**
 * \file TTbarRecoTrainPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a class to create tuples to train a BNN for ttbar MVA reconstruction.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <BTagger.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <string>
#include <vector>
#include <memory>


/**
 * \class Interpretation
 * \brief Auxiliary class to describe an interpretation of an event
 */
struct Interpretation
{
    /// Index of b-jet from the top-quark that decays semi-leptonically
    unsigned bTopLep;
    
    /// Index of b-jet from the top-quark that decays hadronically
    unsigned bTopHad;
    
    /// Indices of light-flavour jets from hadronically decaying top-quark
    unsigned q1TopHad, q2TopHad;
    
    /// Distanse to generator-level configuration
    double distance;
    
    /**
     * \brief Four-momenta of reconstructed top-quarks
     * 
     * It is not needed to define an interpretation, but it is saved for computational efficiency.
     */
    TLorentzVector p4RecoTopLep, p4RecoTopHad;
    
    /**
     * \brief Four-momentum of reconstructed W boson from hadronically decaying top-quark
     * 
     * It is not needed to define an interpretation, but it is saved for computational efficiency.
     */
    TLorentzVector p4RecoWHad;
};


/**
 * \class TTbarRecoTrainPlugin
 * \brief Creates tuples to train a BNN for ttbar MVA reconstruction
 * 
 * There are four energetic jets in the final state of semileptonic ttbar. The class considers all
 * the ways to choose four jets out of all jets in an event (combination of PECReader::GetJets() and
 * PECReader::GetAdditionalJets()). All the possible associations of these four jets to parent
 * objects (two top-quarks) are constructed. Such combinatorics defines a set of interpretations of
 * an event. For each interpretations a distance between reconstructed and generator-level
 * top-quarks and light-flavour jets from hadronically decaying top-quark is calculated (metrics
 * $\Delta R \oplus \Delta p_T^{rel}$ is used). The interpretations are ordered in the distance (in
 * increasing order). Finally, a set of reconstruction-level observables is calculated for each
 * interpretation and stored in a ROOT tree (along with information about the distance). The user
 * can also configure the plugin to store only one interpretation per event (see description of
 * TTbarRecoTrainPlugin::pruned data member for details).
 * 
 * Some details are available a talk by Andrey Popov in [1]; however, the procedure has evolved a
 * bit since that time.
 * [1] https://indico.cern.ch/conferenceDisplay.py?confId=250153 (password is "tHmeeting")
 */
class TTbarRecoTrainPlugin: public Plugin
{
    public:
        /// Constructor
        TTbarRecoTrainPlugin(std::string const &outDirectory, BTagger const &bTagger,
         bool pruned = false);

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
        /// Calculates reconstruction-level observables for a given event interpretation
        void CalculateRecoVars(Interpretation const &interpr, Lepton const &lepton,
         TLorentzVector const &p4RecoWLep, double Ht);
    
    private:
        /// Pointer to PECReaderPlugin
        PECReaderPlugin const *reader;
        
        /// An object to perform b-tagging
        BTagger const &bTagger;
        
        /// Directory to store output files
        std::string outDirectory;
        
        /**
         * \brief Indicates if the output tree must be pruned
         * 
         * In the default workflow all the interpretations of an event are written to the output
         * tree. But if this flag is set to true (by default it is false), only one interpretation
         * per event is stored. The interpretation is chosen randomly. The best interpretation is
         * picked up with a probability of 0.5; the rest of probability is shared uniformly among
         * the background combinations.
         */
        bool pruned;
        
        /// Random-number generator
        std::unique_ptr<TRandom3> rGen;
        
        /// Vector off all jets in an event
        std::vector<Jet> allJets;
        
        /// Vectory of indices of jets that are considered in the current interpretation (unordered)
        std::vector<unsigned> unmaskedJetIndices;
        
        /// Vector to store all interpretations of the current event
        std::vector<Interpretation> interpretations;
        
        /// Current output file
        TFile *file;
        
        /// Current output tree
        TTree *tree;
        
        // Output buffers
        ULong64_t eventNumber, runNumber, lumiSection;
        
        Float_t NJets30, NTags30;
        //^ Needed to allow the same tuples for 3t and 4t bins
        
        Int_t InterpretationRank;
        //^ 0 for the best interprothesis, 2 for the worst one, 1 for the rest of them. The branch
        //is used mostly to demarcate interpretations of different events
        
        Float_t Distance;
        
        Float_t MassTopLep, PtTopLep, EtaTopLep;
        Float_t MassTopHad, PtTopHad, EtaTopHad;
        Float_t MassWHad, PtWHad, EtaWHad;
        
        Float_t DeltaRTopTop;
        Float_t DeltaRTopLepWLep, DeltaRTopHadWHad, DeltaRBJetTopLepWLep, DeltaRBJetTopHadWHad;
        
        Float_t DeltaRLightJets;
        
        Float_t MinEtaTop, MaxEtaTop;
        Float_t DEtaTopTop;
        
        Float_t CosLepTopLepWLep;
        //^ Angle between 3-momenta of the lepton and the top-quark decaying leptonically in the rest
        //frame of the daughter W-boson
        
        Float_t RelHt;
        //^ Sum of pt(t) and pt(h) devided by the total Ht of the event
        
        Float_t MinPtBJet, MinPtLightJet;
        
        Float_t PassBTagTopLep, PassBTagTopHad;
        //^ Whether presumable b-jets from decays of each top-quark are b-tagged
        
        Float_t NLightPassBTagTopHad;
        //^ Number of presumably light-flavour jets from hadronic decay of a top-quark that pass b-tag
        
        Float_t CSVBJetTopLep, CSVBJetTopHad, MaxCSVLightJetsTopHad;
        
        Float_t weight;
};