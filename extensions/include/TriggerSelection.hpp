/**
 * \file TriggerSelection.hpp
 * \author Andrey Popov
 * 
 * The module implements classes to perform a trigger selection in case of real data and simulation.
 */

#pragma once

#include <TriggerSelectionInterface.hpp>
#include <TriggerRange.hpp>
#include <PECReader.hpp>

#include <Rtypes.h>

#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <stdexcept>


/**
 * \class TriggerSelectionData
 * \brief Implements a generic trigger selection for real data
 * 
 * The selection is described by a collection of TriggerRange objects. They must provide a valid
 * data trigger for each event in the dataset that will be processed; if the trigger corresponding
 * to a specific event is not found in the trigger tree, an exception is thrown.
 * 
 * Only necessary branches of the trigger tree are read. The GetWeight method only performs an
 * (optinal) event selection specified in the corresponding TriggerRange object; thus, the weight
 * returned is either 0. or 1.
 * 
 * The class is not thread-safe.
 */
class TriggerSelectionData: public TriggerSelectionInterface
{
    public:
        /**
         * \brief Constructs an instance from a collection of TriggerRange objects
         * 
         * The argument is a reference to a collection of TriggerRange objects. The collection must
         * be a valid container. Internally, only pointers to the TriggerRange objects are stored;
         * therefore, the objects must be available during the lifetime of this. They are not owned
         * by this.
         * 
         * The trigger ranges given to the constructor are expected not to overlap. Otherwise the
         * behaviour of trigger selection is undefined.
         */
        template<typename C>
        TriggerSelectionData(C const &ranges);
        
        /**
         * \brief Constructs an instance from a vector of pointers to TriggerRange objets
         * 
         * The objects must be available during the lifetime of this. They are not owned by this.
         */
        TriggerSelectionData(std::vector<TriggerRange const *> const &ranges);
    
    public:
        /**
         * \brief Resets the pointer to the trigger tree
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method. The value of second argument is ignored.
         * 
         * Implimentation in this class simply updates the pointer and corresponding counters.
         */
        virtual void UpdateTree(TTree *triggerTree, bool);
        
        /**
         * \brief Reads next event from the trigger tree
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * Before reading the next entry from the tree, the method checks if the trigger range
         * corresponding to the given event ID is different from the one used with previous event.
         * If it is the case, it changes the branch which is read from the tree.
         */
        virtual bool ReadNextEvent(EventID const &eventID);
        
        /**
         * \brief Checks if the current event is accepted by the corresponding trigger
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         */
        virtual bool PassTrigger() const;
        
        /**
         * \brief Checks if an event passes an additional selection
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * The method returns either 0. or 1.
         */
        virtual double GetWeight(PECReader const &reader) const;
        
        /**
         * \brief Clones this
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         */
        virtual TriggerSelectionInterface *Clone() const;
    
    private:
        /**
         * \brief Pointers to TriggerRange objects that define the trigger selection
         * 
         * The pointed-to objects must be available during the lifetime of this. They are not owned
         * by this.
         */
        std::vector<TriggerRange const *> ranges;
        
        /**
         * \brief Pointer to the current TriggerRange object
         * 
         * If the current event does not match any TriggerRange object or no attempt to find the
         * corresponding range has been made yet, the poiter is null.
         */
        TriggerRange const *currentRange;
        
        /// The buffer into which results of the selected trigger are read
        Bool_t eventAccepted;
};


template<typename C>
TriggerSelectionData::TriggerSelectionData(C const &ranges_):
    TriggerSelectionInterface(),
    currentRange(nullptr)
{
    // First a sanity check
    if (ranges_.size() == 0)
        throw std::logic_error("TriggerSelectionData::TriggerSelectionData: The provided "
         "collection of TriggerRange objects is empty.");
    
    
    // Save pointers to the TriggerRange objects
    ranges.reserve(ranges_.size());
    
    for (typename C::const_iterator it = ranges_.cbegin(); it != ranges_.cend(); ++it)
        ranges.push_back(&*it);
}


/**
 * \class TriggerSelectionMC
 * \brief Implements a generic trigger selection for simulation
 * 
 * The selection is described by a collection of TriggerRange objects. The MC triggers specified
 * there must be available in all datasets analysed with this class; otherwise an exception is
 * thrown.
 * 
 * An event passes the first step of the trigger selection if it is accepted by at least one of
 * the MC triggers. At the second step it is assigned a weight calculated as a sum of integrated
 * luminosities for all triggers that accept the event. It might further be corrected by the
 * trigger scale factors (not implemented yet). In addition, each TriggerRange object might define
 * an optional event selection, which is also evaluated when the weight is calculated.
 * 
 * The class is not thread-safe.
 */
class TriggerSelectionMC: public TriggerSelectionInterface
{
    public:
        /**
         * \brief Constructs an instance from a collection of TriggerRange objects
         * 
         * The argument is a reference to a collection of TriggerRange objects. The collection must
         * be iterable. Internally, only pointers to the TriggerRange objects are stored; therefore,
         * the objects must be available during the lifetime of this. They are not owned by this.
         */
        template<typename C>
        TriggerSelectionMC(C const &ranges);
        
        /**
         * \brief Constructs an instance from a vector of pointers to TriggerRange objets
         * 
         * The objects must be available during the lifetime of this. They are not owned by this.
         */
        TriggerSelectionMC(std::vector<TriggerRange const *> const &ranges);
        
        /// Destructor
        ~TriggerSelectionMC();
    
    public:
        /**
         * \brief Resets the pointer to the trigger tree
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method. The value of second argument is ignored.
         * 
         * Updates the pointer and corresponding counters. Sets addresses for branches of the
         * new trigger tree containining information about all the MC triggers specified in the
         * TriggerRange objects. The second argument is ignored.
         */
        virtual void UpdateTree(TTree *triggerTree, bool);
        
        /**
         * \brief Reads next event from the trigger tree
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * Reads all the active branches to the corresponding buffers. Throws an exception in case
         * the tree has not been specified. The argument is ignored.
         */
        virtual bool ReadNextEvent(EventID const &);
        
        /**
         * \brief Checks if the current event is accepted by the corresponding triggers
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * The method returns true if at least one of the specified MC triggers accepts the event.
         */
        virtual bool PassTrigger() const;
        
        /**
         * \brief Checks if an event passes an additional selection
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * All the TriggerRange objects are considered. The weight is calculated as a sum of
         * integrated luminosities multiplied by the scale factors; the sum is taken over those
         * ranges for which the event passes the additional selection and is accepted by the
         * corresponding MC trigger. No systematical variation is implemented.
         */
        virtual double GetWeight(PECReader const &reader) const;
        
        /**
         * \brief Clones this
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         */
        virtual TriggerSelectionInterface *Clone() const;
    
    protected:
        /**
         * \brief Calculates the trigger scale factor
         * 
         * This is a placeholder which always returns 1. The user can override this method. The
         * functionality of the class should also be extended to allow to specify arbitrary
         * callable expression, similar to how it is done in TriggerRange class. Note the
         * systematical variations are not implemented and are likely to require changes in the
         * trigger interface.
         */
        virtual double ScaleFactor(TriggerRange const *range, PECReader const &reader) const;
    
    private:
        /**
         * \brief Pointers to TriggerRange objects that define the trigger selection and buffers
         * that store trigger states
         * 
         * The pointed-to TriggerRange objects must be available during the lifetime of this. They
         * are not owned by this.
         * 
         * The buffers to which pair::second point store decisions of the corresponding MC triggers.
         * All the buffers are allocated in a block, to which the buffer variable points. Such
         * organisation is preferred to using Bool_t in pair::second because the same MC trigger
         * can be specified in several trigger ranges.
         */
        std::vector<std::pair<TriggerRange const *, Bool_t *>> ranges;
        
        /**
         * \brief Buffer into which the trigger decision are read
         * 
         * Points to a dinamycally allocated array. The buffers are normally accessed via pointers
         * in ranges::value_type::second. This variable is only used to initialise them and set
         * branches' addresses.
         */
        Bool_t *buffer;
        
        TFile *sf_file;
        TH1D *hist_id, *hist_iso, *hist_trig;
};


template<typename C>
TriggerSelectionMC::TriggerSelectionMC(C const &ranges_):
    TriggerSelectionInterface(),
    buffer(nullptr)
{
    // First a sanity check
    if (ranges_.size() == 0)
        throw std::logic_error("TriggerSelectionMC::TriggerSelectionMC: The provided "
         "collection of TriggerRange objects is empty.");
    
    
    // Save pointers to the TriggerRange objects
    ranges.reserve(ranges_.size());
    
    for (typename C::const_iterator it = ranges_.cbegin(); it != ranges_.cend(); ++it)
        ranges.emplace_back(&*it, nullptr);
}


/**
 * \class TriggerSelection
 * \brief A class to perform a generic trigger selection
 * 
 * When the method UpdateTree is called, the class creates an instance of either TriggerSelectionMC
 * or TriggerSelectionData, depending on the context. Later on it works as a simple proxy for this
 * selection object. The tree pointer and corresponding counters inherited from
 * TriggerSelectionInterface are never used.
 * 
 * The class is not thread-safe.
 */
class TriggerSelection: public TriggerSelectionInterface
{
    public:
        /**
         * \brief Constructs an instance from a collection of TriggerRange objects
         * 
         * The argument is a reference to a collection of TriggerRange objects. The collection must
         * be iterable. Internally, only pointers to the TriggerRange objects are stored; therefore,
         * the objects must be available during the lifetime of this. They are not owned by this.
         */
        template<typename C>
        TriggerSelection(C const &ranges);
        
        /**
         * \brief Constructs an instance from a vector of pointers to TriggerRange objets
         * 
         * The objects must be available during the lifetime of this. They are not owned by this.
         */
        TriggerSelection(std::vector<TriggerRange const *> const &ranges);
    
    public:
        /**
         * \brief Resets the pointer to the trigger tree
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method. The value of second argument is ignored.
         * 
         * Creates an object to perform the trigger selection on data or on simulation depending on
         * the provided flag.
         */
        virtual void UpdateTree(TTree *triggerTree, bool isData);
        
        /**
         * \brief Reads next event from the trigger tree
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * Simply calls ReadNextEvent of the selection object.
         */
        virtual bool ReadNextEvent(EventID const &eventID);
        
        /**
         * \brief Checks if the current event is accepted by the corresponding triggers
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * Simply calls PassTrigger of the selection object.
         */
        virtual bool PassTrigger() const;
        
        /**
         * \brief Checks if an event passes an additional selection
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         * 
         * All the TriggerRange objects are considered. The weight is calculated as a sum of
         * integrated luminosities multiplied by the scale factors; the sum is taken over those
         * ranges for which the event passes the additional selection and is accepted by the
         * corresponding MC trigger. No systematical variation is implemented.
         */
        virtual double GetWeight(PECReader const &reader) const;
        
        /**
         * \brief Clones this
         * 
         * Consult documentation of the overridden method in the base class for a description of
         * the purpose of this method.
         */
        virtual TriggerSelectionInterface *Clone() const;
    
    private:
        /**
         * \brief Pointers to TriggerRange objects that define the trigger selection
         * 
         * The pointed-to objects must be available during the lifetime of this. They are not owned
         * by this.
         */
        std::vector<TriggerRange const *> ranges;
        
        /**
         * \brief An object that performs the actual trigger selection
         * 
         * It is an instance of either TriggerSelectionMC or TriggerSelectionData.
         */
        std::unique_ptr<TriggerSelectionInterface> selection;
};


template<typename C>
TriggerSelection::TriggerSelection(C const &ranges_):
    TriggerSelectionInterface()
{
    // First a sanity check
    if (ranges_.size() == 0)
        throw std::logic_error("TriggerSelection::TriggerSelection: The provided "
         "collection of TriggerRange objects is empty.");
    
    
    // Save pointers to the TriggerRange objects
    ranges.reserve(ranges_.size());
    
    for (typename C::const_iterator it = ranges_.cbegin(); it != ranges_.cend(); ++it)
        ranges.push_back(&*it);
}