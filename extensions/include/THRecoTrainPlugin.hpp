/**
 * \file THRecoTrainPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a class to create tuples to train a BNN for thq MVA reconstruction.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <BTagger.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <string>
#include <vector>


/**
 * \class Interpretation
 * \brief Auxiliary class to describe an interpretation of an event
 */
struct Interpretation
{
    /// Index of b-jet from the top-quark
    unsigned bTop;
    
    /// Index of recoil light-flavour jet
    unsigned qRecoil;
    
    /// Index of b-jets from the higgs decay
    unsigned b1Higgs, b2Higgs;
    
    /// Distanse to generator-level configuration
    double distance;
    
    /**
     * \brief Four-momentum of reconstructed top-quark
     * 
     * It is not needed to define an interpretation, but it is saved for computational efficiency.
     */
    TLorentzVector p4RecoTop;
    
    /**
     * \brief Four-momentum of reconstructed Higgs boson
     * 
     * It is not needed to define an interpretation, but it is saved for computational efficiency.
     */
    TLorentzVector p4RecoHiggs;
};


/**
 * \class THRecoTrainPlugin
 * \brief Creates tuples to train a BNN for thq MVA reconstruction
 * 
 * There are four energetic jets in the final state of thq process. The class considers all the ways
 * to choose four jets out of all jets in an event (combination of PECReader::GetJets() and
 * PECReader::GetAdditionalJets()). All the possible associations of these four jets to parent
 * objects (top-quark, Higgs boson, or recoil quark) are constructed. Such combinatorics defines a
 * set of interpretations of an event. For each interpretations a distance between reconstructed and
 * generator-level top-quark, Higgs boson, and recoil quark is calculated (metrics $\Delta R \oplus
 * \Delta p_T^{rel}$ is used). The interpretations are ordered in the distance (in increasing
 * order). Finally, a set of recontruction-level observables is calculated for each interpretation
 * and stored in a ROOT tree (along with information about the distance).
 * 
 * Details are available a talk by Andrey Popov in [1].
 * [1] https://indico.cern.ch/conferenceDisplay.py?confId=251808 (password is "tHmeeting")
 */
class THRecoTrainPlugin: public Plugin
{
    public:
        /// Constructor
        THRecoTrainPlugin(std::string const &outDirectory, BTagger const &bTagger);

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
        
        /// An object to perform b-tagging
        BTagger const &bTagger;
        
        /// Directory to store output files
        std::string outDirectory;
        
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
        
        Int_t InterpretationRank;
        //^ 0 for the best interprothesis, 2 for the worst one, 1 for the rest of them. The branch
        //is used mostly to demarcate interpretations of different events
        
        Float_t Distance;
        
        Float_t MassTop, PtTop, EtaTop;
        Float_t MassHiggs, PtHiggs, EtaHiggs;
        Float_t PtLJet, EtaLJet;
        
        Float_t DeltaRTopHiggs;
        Float_t DeltaRTopW, DeltaRBJetTopW;
        Float_t DeltaEtaLepTop;
        Float_t DeltaRBJetsHiggs;
        
        Float_t CosLepTopW;
        //^ Angle between 3-momenta of the lepton and the top-quark in the rest frame of the
        //daughter W-boson. Sign of the cosine is inverted w.r.t. to this description because this
        //is the common definition
        
        Float_t RelHt;
        //^ Sum of pt(t) and pt(h) devided by the total Ht of the event
        
        Float_t MinPtBJet;
        
        Float_t PassBTagTop, PassBTagLJet;
        //^ Whether the psesumable b-jet from top-quark decay and the light-flavoured recoil jets
        //are b-tagged
        
        Float_t NPassBTagHiggs;
        //^ Number of presumable b-jets from Higgs boson decay that pass b-tag
        
        Float_t CSVBJetTop, CSVLJet, MinCSVBJetsHiggs;
        
        Float_t weight;
};