/**
 * \file SingleTopTChanBTagEffPlugin.hpp
 * \author Andrey Popov, Natalia Tsirova
 * 
 * The module defines a class to calculate a set of variables that are used in single-top t-channel
 * analysis.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>

#include <TFile.h>
#include <TTree.h>

#include <string>


/**
 * \class SingleTopTChanBTagEffPlugin
 * \brief Calculates and stores b-tag efficiency for event selection used in SingleTop t-chan analysis
 * 
 * The class is expected to serve as an illustration rather than to be used in a real-life analysis.
 */
class SingleTopTChanBTagEffPlugin: public Plugin
{
    public:
        /// Constructor
        SingleTopTChanBTagEffPlugin(std::string const &outDirectory);
    
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
        
        /// Directory to store output files
        std::string outDirectory;
        
        /// Current output file
        TFile *file;
        
        /// Output histograms
        TH2D *histB, *histTagB, *histC, *histTagC, *histUDS, *histTagUDS, *histG, *histTagG;
        
        // Output buffers
        ULong64_t eventNumber, runNumber, lumiSection;

        Float_t weight;
};
