/**
 * \file THRecoPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a class to perform MVA reconstruction of thq events.
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
#include <memory>


#include <bnn_thq_reco_3t.hpp>


/**
 * \class THRecoPlugin
 * \brief Performs an MVA reconstruction of a thq event
 * 
 * There are four energetic jets in the final state of thq process. The class considers all the ways
 * to choose four jets out of all jets in an event (combination of PECReader::GetJets() and
 * PECReader::GetAdditionalJets()). All the possible associations of these four jets to parent
 * objects (top-quark, Higgs boson, or recoil quark) are constructed. Such combinatorics defines a
 * set of interpretations of an event. A likelihood of each interpretation is calculated with the
 * help of a dedicated BNN. Finally, the interpretation with the highest BNN score is chosen. It
 * defines four-momenta of the reconstructed top-quark and Higgs boson, which are accessible by the
 * user via dedicated getters.
 * 
 * Details are available a talk by Andrey Popov in [1].
 * [1] https://indico.cern.ch/conferenceDisplay.py?confId=251808 (password is "tHmeeting")
 */
class THRecoPlugin: public Plugin
{
    public:
        /**
         * \struct Interpretation
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
    
    public:
        /// Constructor
        THRecoPlugin(BTagger const &bTagger);

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
        
        /// Returns reconstructed top-quark
        Candidate const &GetRecoTopQuark() const;
        
        /// Returns reconstruced Higgs boson
        Candidate const &GetRecoHiggsBoson() const;
        
        /// Returns reconstructed jet from the recoil light-flavour quark
        Candidate const &GetRecoRecoilQuark() const;
        
        /**
         * \brief Returns the best interpretation
         * 
         * Expected to be used on rare occasions when four-momenta of reconstruced objects are not
         * sufficient.
         */
        Interpretation const &GetInterpretation() const;
    
    private:
        /// Calculates reconstruction-level observables for a given event interpretation
        void CalculateRecoVars(Interpretation const &interpr, Lepton const &lepton,
         TLorentzVector const &p4RecoW, double Ht);
    
    private:
        /// Pointer to PECReaderPlugin
        PECReaderPlugin const *reader;
        
        /// An object to perform b-tagging
        BTagger const &bTagger;
        
        /// Vector off all jets in an event
        std::vector<Jet> allJets;
        
        /// Vectory of indices of jets that are considered in the current interpretation (unordered)
        std::vector<unsigned> unmaskedJetIndices;
        
        /// The neural network to perform reconstruction
        thq_reco_3t::BNN bnnReco;
        
        /// Best interpretation of the current event
        Interpretation bestInterpretation;
        
        /// Reconstructed top quark
        Candidate recoTopQuark;
        
        /// Reconstructed Higgs boson
        Candidate recoHiggsBoson;
        
        /// Reconstructed recoil quark
        Candidate recoRecoilQuark;
        
        
        // Input variables (not all of them are used in the BNN)
        Float_t MassTop, PtTop, EtaTop;
        Float_t MassHiggs, PtHiggs, EtaHiggs;
        Float_t PtLJet, EtaLJet;
        
        Float_t DeltaRTopHiggs;
        Float_t DeltaRTopW, DeltaRBJetTopW;
        Float_t DeltaEtaLepTop;
        Float_t DeltaRBJetsHiggs;
        
        Float_t RelHt;
        //^ Sum of pt(t) and pt(h) devided by the total Ht of the event
        
        Float_t MinPtBJet;
        
        Float_t PassBTagTop;
        //^ Whether the psesumable b-jet from top-quark decay is b-tagged
        
        Float_t NPassBTagHiggs;
        //^ Number of presumable b-jets from Higgs boson decay that pass b-tag
};