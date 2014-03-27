/**
 * \file WjetsHFPlugin.hpp
 * \author Andrey Popov, Natalia Tsirova
 * 
 * Defines a plugin to classify a wjets event depending on heavy flavours.
 */

#pragma once


#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <PhysicsObjects.hpp>

#include <list>
#include <vector>


/**
 * \class WjetsHFPlugin
 * \brief Classifies a wjets event depending on flavours.
 * 
 * Classification is based on partons after parton shower and reconstructed jets. Supported
 * categories are listed in the enumeration Type.
 */
class WjetsHFPlugin: public Plugin
{
public:
    /**
     * \brief Supported categories to classify an event
     * 
     * The classification is top-to-bottom: if an event matches some category, all the following
     * ones are not checked.
     */
    enum class Type
    {
        W_qq, ///< W+bbar or W+ccbar
        W_c,         ///< W+c
        W_other,   ///< W+heavy quark which is an immediate proton's daughter
        W_light,      ///< No heavy flavours
        
    };

public:
    /// Constructor
    WjetsHFPlugin(std::string const &name) noexcept;
    
public:
    /// Creates a newly-initialised copy of *this
    virtual Plugin *Clone() const;
    
    /// Saves a pointer to PECReaderPlugin for convenience
    virtual void BeginRun(Dataset const &dataset);
    
    /**
     * \brief Classifies the current event
     * 
     * Always returns true.
     */
    virtual bool ProcessEvent();
    
    /// Returns the classification decision
    Type GetDecision() const noexcept;
    
private:
    /// Pointer to the reader plugin
    PECReaderPlugin const *reader;
    
    /// Decision on classification of the current event
    Type type;
   
};
