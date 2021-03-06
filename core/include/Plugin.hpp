/**
 * \file Plugin.hpp
 * \author Andrey Popov
 * 
 * The module describes an interface for a plugin for Processor.
 */

#pragma once

#include <PluginForward.hpp>

#include <ProcessorForward.hpp>
#include <Dataset.hpp>

#include <string>


/**
 * \class Plugin
 * \brief An abstract class to define a plugin to be used in class Processor
 * 
 * The class specifies a simple interface that allows to perform a certain processing for each event
 * in a dataset and, additionally, notify the plugin when processing of the dataset starts of
 * finishes. It contains a pointer to the parent Processor instance that can be used to access
 * other plugins by their name.
 * 
 * A single instance of a derived class might be used to process several files, which requires a
 * certain initialization and termination actions to be performed in BeginRun/EndRun methods in
 * case of non-trivial plugins. The pointer to the parent Processor instance is initialized before
 * the first file is read and is expected to be valid during the lifetime of the plugin.
 * 
 * A derived class might or might not be copyable, but it must implement a method to clone an
 * instance. The clonning must address only configuration of the processing algorithm but not data
 * members specific for a dataset or an event (e.g. handlers of output files). Such a functionality
 * is requered to multiplicate the plugin structure for each thread represented by class Processor.
 * 
 * A derived class must be capable of working in a multi-thread mode. The user should pay attention
 * to the fact that ROOT is not thread-safe. For this reason all the critical blocks (which include,
 * for example, creation of any ROOT objects) must be guarded with the help of class ROOTLock.
 * 
 * A derived class must define a valid move constructor.
 */
class Plugin
{
    public:
        /// Constructor
        Plugin(std::string const &name);
        
        /// Default copy constructor
        Plugin(Plugin const &) = default;
        
        /// Default move constructor
        Plugin(Plugin &&) = default;
        
        /// Default assignment operator
        Plugin &operator=(Plugin const &) = default;
        
        /// Trivial destructor
        virtual ~Plugin();
    
    public:
        /**
         * \brief Provides a pointer to an instance of Processor class that owns the plugin
         * 
         * The pointer is guaranteed to be initialized before the first call to BeginRun. It stays
         * valid for the lifetime of the object.
         */
        void SetParent(Processor const *processor);
        
        /// Returns name of the plugin
        std::string const &GetName() const;
        
        /**
         * \brief Clones the object. Must be implemented by the user
         * 
         * The method must create a new instance of the (derived) class with the same constructor
         * parameters. The method must not address any parameters specific to a run or an event.
         * Thechnically it means the method must create a new instance of the class exactly in the
         * same way this has been created and initialized.
         * 
         * The method is used when unique copies of plugins are created for each instance of class
         * Processor. Clonning is performed before call to SetParent and before the first call to
         * BeginRun.
         */
        virtual Plugin *Clone() const = 0;
        
        /**
         * \brief Called before processing of a new dataset is started
         * 
         * The method is trivial in the default implementation
         */
        virtual void BeginRun(Dataset const &dataset);
        
        /**
         * \brief Called after processing of a dataset is finished
         * 
         * The method is trivial in the default implementation
         */
        virtual void EndRun();
        
        /**
         * \brief Called for each event. Must be implemented by the user
         * 
         * The boolean return value allows to implement event filtering. The plugin must return
         * true unless it suggests the event to be discarded.
         */
        virtual bool ProcessEvent() = 0;
    
    protected:
        /// Unique name to identify the plugin
        std::string const name;
        
        /// Parent Processor object
        Processor const *processor;
};