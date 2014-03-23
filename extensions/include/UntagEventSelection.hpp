/**
 * \file UntagEventSelection.hpp
 * \author Andrey Popov, Natalia Tsirova 
 * The module defines a class to implement a generic event selection with no bTag requirements.
 */

#pragma once

#include <EventSelectionInterface.hpp>


#include <list>
#include <array>
#include <map>
#include <memory>


/**
 * \class UntagEventSelection
 * \brief Allows to implement an event selection in a flexible way
 * 
 * Along with the constructor, the user must employ methods AddLeptonThreshold and AddJetTagBin
 * in order to specify the event selection.
 * 
 * Each instance of class PECReader must exploit its own copy of class UntagEventSelection, which
 * can be achieved with the help of method Clone.
 * 
 * Consult documentation for the base class for details of the interface.
 * 
 * The class is copyable.
 */
class UntagEventSelection: public EventSelectionInterface
{
    private:
        /**
         * \struct JetTagBin
         * \brief Defines a single bin for jet selection
         */
        struct JetBin
        {
            /// Default constructor
            JetBin() = default;
            
            /// Constructor with parameters
            JetBin(unsigned nJets_);
            
            unsigned nJets;  ///< Number of jets
            
        };
    
    public:
        /**
         * \brief Constructor
         * 
         * Parameter is the jet pt threshold.
         */
        UntagEventSelection(double jetPtThreshold_);
        
        /// Default copy constructor
        UntagEventSelection(UntagEventSelection const &) = default;
        
        /// Assignment operator is deleted
        UntagEventSelection &operator=(UntagEventSelection const &) = delete;
    
    public:
        /**
         * \brief Performs the event selection on leptons
         * 
         * Checks the number and the transverse momenta of the tight leptons, vetoes additional
         * loose leptons. For example, if the user has requested two muons with thresholds 25 and
         * 15 GeV/c with the help of method AddLeptonThreshold, then an event will be checked to
         * contain exactly two tight muons with pt > 15 GeV/c and, in addition, exactly one of them
         * will be required to have pt > 25 GeV/c.
         * 
         * \note See also the documentation in the base base class.
         */
        virtual bool PassLeptonStep(std::vector<Lepton> const &tightLeptons,
         std::vector<Lepton> const &looseLeptons) const;
        
        /**
         * \brief Performs the event selection on jets
         * 
         * An event is checked against the allowed jet-tag bins. The input collection must be
         * filtered with the help of IsAnalysisJet method. To find the number of b-tagged jets the
         * specified b-tagging object is used.
         */
        virtual bool PassJetStep(std::vector<Jet> const &jets) const;
        
        /**
         * \brief Checks if a jet is to be used in high-level analysis
         * 
         * See also documentation of the overridden method in the base class.
         */
        virtual bool IsAnalysisJet(Jet const &jet) const;
        
        /**
         * \brief Adds additional lepton to the selection
         * 
         * The method increases the number of required leptons of the given flavour by one and
         * sets the pt threshold for the new lepton.
         * 
         * \note The method is implemented in such a way that the underlying lists of thresholds
         * are ordered at every moment.
         */
        void AddLeptonThreshold(Lepton::Flavour flavour, double ptThreshold);
        
        /// Extends the selection on jets with a given jet bin
        void AddJetBin(unsigned nJets);
        
        /**
         * \brief Creates a newly-initialized copy of this
         * 
         * Consult documentation for the base class for details.
         */
        EventSelectionInterface *Clone() const;
    
    private:
        /// Map from lepton flavours to integers starting from zero
        std::map<Lepton::Flavour, unsigned> flavourMap;
        
        /// Iterators of the current positions in the lists of lepton pt thresholds
        mutable std::array<std::list<double>::const_iterator, 3> leptonThresholdIts;
        
        /// The three lists of lepton pt thresholds (one for each lepton flavour)
        std::array<std::list<double>, 3> leptonPtThresholds;
        
        double jetPtThreshold;  ///< Minimum pt for analysis-level jets
        std::list<JetBin> jetBins;  ///< List of allowed jet-tag bins
};