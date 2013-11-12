#include <EventCounter.hpp>

#include <stdexcept>
#include <iostream>
#include <fstream>


using namespace std;


SelectionStep::SelectionStep():
    description(""),
    nPassed(0)
{}


SelectionStep::SelectionStep(string const &description_):
    description(description_),
    nPassed(0)
{}


SelectionStep::SelectionStep(SelectionStep &&src):
    description(move(src.description)),
    eventIDs(move(src.eventIDs))
{}


EventCounter::EventCounter(string const &title_):
    title(title_)
{}


void EventCounter::RegisterSelectionStep(string const &label, string const &description)
{
    // Make sure the label is not assigned to some selection step already
    if (steps.find(label) not_eq steps.end())
        throw logic_error(string("EventCounter::RegisterSelectionStep: Cannot register a selection "
         "step with label \"" + label + "\" as it is already in use."));
    
    
    steps[label].description = description;
}


void EventCounter::AddEvent(string const &label)
{
    // Make sure a selection step with the given label has been defined
    if (steps.find(label) == steps.end())
        throw logic_error(string("EventCounter::AddEvent: Selection step with label \"" +
         label + "\" is not defined."));
    
    
    ++steps[label].nPassed;
}


void EventCounter::AddEvent(string const &label, EventID const &id)
{
    AddEvent(label);
    steps[label].eventIDs.emplace_back(id);
}


void EventCounter::WriteResults() const
{
    // Write down the cut-flow table
    cout << "Cut-flow table for counter \"" << title << "\"\n";
    
    for (auto const &s: steps)
        cout << s.first << " (" << s.second.description << ")\n " << s.second.nPassed << '\n';
    
    cout << '\n';
    
    
    // Create files with event IDs
    for (auto const &s: steps)
    {
        // Skip the step if no IDs are saved
        if (s.second.eventIDs.size() == 0)
            continue;
        
        // Make sure the user provided (or did not) event IDs for this step consistently
        if (s.second.nPassed not_eq s.second.eventIDs.size())
        {
            cout << "Warning: Number of events that passed selection step \"" << s.first <<
             "\" does not match the number of saved IDs. This step will be skipped.\n";
            continue;
        }
        
        
        // Open a file to keep the IDs and write a header
        ofstream file(title + '_' + s.first + ".txt");
        file << "Selection step \"" << s.first << "\": " << s.second.description << "\n\n";
        file << "Run   Luminosity block   Event number\n";
        
        for (auto const &id: s.second.eventIDs)
            file << id.Run() << "  " << id.LumiBlock() << "  " << id.Event() << '\n';
        
        file.close();
    }
}