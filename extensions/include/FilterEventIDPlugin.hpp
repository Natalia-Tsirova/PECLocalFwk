/**
 * \file FilterEventIDPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a class that filters event based on their ID.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <EventID.hpp>

#include <string>
#include <map>
#include <vector>


/**
 * \class FilterEventIDPlugin
 * \brief Filters events on their ID
 * 
 * The class filters events based on their ID (run, lumi, event numbers). The user specifies a list
 * of event IDs in a text file. [Describe the format.]
 * 
 * \warning The class can filter atomic datasets only (i.e. ones that contain a single file).
 */
class FilterEventIDPlugin: public Plugin
{
    public:
        /// Constructor
        FilterEventIDPlugin(std::string const &name, std::string const &eventIDsFileName,
         bool rejectKnownEvent = true);
    
    private:
        /// A private constructor to be used in method Clone
        FilterEventIDPlugin(std::string const &name, std::map<std::string,
         std::vector<EventID>> const &eventIDsAllFiles, bool rejectKnownEvent);

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
        
        /**
         * \brief Switch defining the logic of filtering
         * 
         * If set to true, the plugin rejects events whose IDs are mentioned in the provided file.
         * If set to false, it keeps only events with IDs written in the file.
         */
        bool rejectKnownEvent;
        
        /// Map to store event IDs. The key of the map is a short name of corresponding ROOT file
        std::map<std::string, std::vector<EventID>> eventIDsAllFiles;
        
        /// Pointer to vector of event IDs for the current ROOT file (note it might be null)
        std::vector<EventID> const *eventIDsCurFile;
};