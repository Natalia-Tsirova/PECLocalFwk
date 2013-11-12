/**
 * \file EventCounter.hpp
 * \author Andrey Popov
 * 
 * Defines a class to build a cut-flow table and (optionally) store event IDs.
 */

#include <EventID.hpp>
 
#include <string>
#include <list>
#include <map>


/**
 * \struct SelectionStep
 * \brief An aggregate to accumulate information about a single step of event selection
 */
struct SelectionStep
{
    /**
     * \brief Default constructor
     * 
     * Needed to put an instance of SelectionStep into a map.
     */
    SelectionStep();
    
    /// Constructor takes a human-readable description of the selection step
    SelectionStep(std::string const &description);
    
    /// Move constructor
    SelectionStep(SelectionStep &&src);
    
    /// Copy constructor is deleted
    SelectionStep(SelectionStep const &) = delete;
    
    /// Assignment operator is deleted
    SelectionStep &operator=(SelectionStep const &) = delete;

    /// A human-readable description of the selection step
    std::string description;
    
    /// Number of events passing the step
    unsigned long nPassed;
    
    /// A list of IDs of events that pass the selection step
    std::list<EventID> eventIDs;
};


/**
 * \class EventCounter
 * \brief A class to fill a cut-flow table
 * 
 * An instance of the class allows to define several selection steps, each with an associated
 * counter and an optional list of event IDs. After the counters are filled, the cut-flow table is
 * written to the standard output and the lists of passed events are stored in a set of files (one
 * file for each selection step).
 * 
 * The class is not thread-safe. It is expected to be used for synchronisation only and not intended
 * to be exploited in a production version of a program. Given its transient usage, the foreseen way
 * to insert the class in the code is to make a global object in a dedicated branch. After the
 * synchronisation is over, the branch should be removed.
 */
class EventCounter
{
    public:
        /// Constructor from a title to indentify an instance of the class
        EventCounter(std::string const &title);
        
        /// Copy constructor is deleted
        EventCounter(EventCounter const &) = delete;
        
        /// Assignment operator is deleted
        EventCounter &operator=(EventCounter const &) = delete;
    
    public:
        /// Registers a new selection step
        void RegisterSelectionStep(std::string const &label, std::string const &description);
        
        /// Increases the counter for a given selection step
        void AddEvent(std::string const &label);
        
        /// Increases the counter for a given selection step and saves event ID
        void AddEvent(std::string const &label, EventID const &id);
        
        /**
         * \brief Writes down accumulated information about selection steps
         * 
         * Cut-flow table is written to the standard output, and lists of event IDs are saved to
         * ASCII files.
         */
        void WriteResults() const;
    
    private:
        /// A title to identify an instance of the class
        std::string const title;
        
        /// A map from labels of selection step (needed to access them) to SelectionStep objets
        std::map<std::string, SelectionStep> steps;
};