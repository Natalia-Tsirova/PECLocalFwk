/**
 * \file TTbarRecoPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a class to perform MVA reconstruction of ttbar events in semileptonic channel
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


#include <bnn_ttbar_reco_3t.hpp>
//#include <bnn_ttbar_reco_CSVM_3t.hpp>


/**
 * \class TTbarRecoPlugin
 * \brief Performs an MVA reconstruction of ttbar events in semileptonic channel
 * 
 * There are four energetic jets in the final state of semileptonic ttbar. The class considers all
 * the ways to choose four jets out of all jets in an event (combination of PECReader::GetJets() and
 * PECReader::GetAdditionalJets()). All the possible associations of these four jets to parent
 * objects (two top-quarks) are constructed. Such combinatorics defines a set of interpretations of
 * an event. A likelihood of each interpretation is calculated by a dedicated BNN, and the one with
 * highest BNN score is chosen. This interpretation defines four-momenta of reconstructed objects,
 * which are accessible by the user.
 * 
 * Some details are available a talk by Andrey Popov in [1]; however, the procedure has evolved a
 * bit since that time.
 * [1] https://indico.cern.ch/conferenceDisplay.py?confId=250153 (password is "tHmeeting")
 */
class TTbarRecoPlugin: public Plugin
{
    public:
        /**
         * \struct Interpretation
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
            
            /**
             * \brief Four-momenta of reconstructed top-quarks
             * 
             * It is not needed to define an interpretation, but it is saved for computational
             * efficiency.
             */
            TLorentzVector p4RecoTopLep, p4RecoTopHad;
            
            /**
             * \brief Four-momentum of reconstructed W boson from hadronically decaying top-quark
             * 
             * It is not needed to define an interpretation, but it is saved for computational
             * efficiency.
             */
            TLorentzVector p4RecoWHad;
        };
            
    public:
        /// Constructor
        TTbarRecoPlugin(BTagger const &bTagger);

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
        
        /// Returns reconstructed semileptonically decaying top-quark
        Candidate const &GetRecoTopQuarkLep() const;
        
        /// Returns reconstructed hadronically decaying top-quark
        Candidate const &GetRecoTopQuarkHad() const;
        
        /// Returns reconstructed W-boson from hadronic decay of a top-quark
        Candidate const &GetRecoWBosonHad() const;
        
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
         TLorentzVector const &p4RecoWLep, double Ht);
    
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
        ttbar_reco_3t::BNN bnnReco;
        //ttbar_reco_CSVM_3t::BNN bnnReco;
        
        /// Reconstructed top-quark that decays semileptonically
        Candidate recoTopQuarkLep;
        
        /// Reconstructed top-quark that decays hadronically
        Candidate recoTopQuarkHad;
        
        /// Reconstructed W-boson from hadronic decay of a top-quark
        Candidate recoWBosonHad;
        
        /// Best interpretation of the current event
        Interpretation bestInterpretation;
        
        
        // Input variables (not all of them are used in the BNN)
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
};