#include <TriggerSelection.hpp>

#include <TTree.h>

#include <algorithm>


using namespace std;


TriggerSelectionData::TriggerSelectionData(vector<TriggerRange const *> const &ranges_):
    TriggerSelectionInterface(),
    ranges(ranges_),
    currentRange(nullptr)
{
    // A sanity check
    if (ranges_.size() == 0)
        throw std::logic_error("TriggerSelectionData::TriggerSelectionData: The provided "
         "collection of pointers to TriggerRange objects is empty.");
}


void TriggerSelectionData::UpdateTree(TTree *triggerTree_, bool)
{
    // Update the tree pointer and counters
    triggerTree = triggerTree_;
    nEntriesTree = triggerTree->GetEntries();
    nextEntryTree = 0;
    
    triggerTree->SetBranchStatus("*", false);
    
    
    // Invalidate the current range
    currentRange = nullptr;
}


bool TriggerSelectionData::ReadNextEvent(EventID const &eventID)
{
    // A sanity check
    if (not triggerTree)
        throw logic_error("TriggerSelectionData::ReadNextEvent: Attempting to read an unspecified "
         "trigger tree.");
    
    
    // Check if there are still events to read
    if (nextEntryTree == nEntriesTree)
        return false;
    
    
    // Check if the current trigger range accommodates the given event ID and update it if needed
    if (not currentRange or not currentRange->InRange(eventID))
    {
        // Find the range that contains the event with the given ID
        auto res = find_if(ranges.begin(), ranges.end(),
         [&eventID](TriggerRange const *r){return (r->InRange(eventID));});
        
        if (res == ranges.end())
        // No range contains the given event ID
        {
            // Set the state of the object such that the event will be rejected
            currentRange = nullptr;
            eventAccepted = false;
            
            // There is no need to actually read the event. Just update the counter
            ++nextEntryTree;
            return true;
        }
        
        
        // A valid trigger range has been found. Get the corresponding branch of the tree and assign
        //the buffer to it
        currentRange = *res;
        triggerTree->SetBranchStatus("*", false);
        
        TBranch *branch =
         triggerTree->GetBranch((currentRange->GetDataTriggerPattern() + "__accept").c_str());
        
        if (not branch)
            throw runtime_error(string("TriggerSelectionData::ReadNextEvent: State of the "
             "trigger \"HLT_") + currentRange->GetDataTriggerPattern() + "_v*\" is not stored in "
             "the source tree.");
        
        branch->SetStatus(true);
        branch->SetAddress(&eventAccepted);
    }
    
    
    // Finally, read the event
    triggerTree->GetEntry(nextEntryTree);
    ++nextEntryTree;
    
    return true;
}


bool TriggerSelectionData::PassTrigger() const
{
    // Everything has already been done in ReadNextEvent
    return eventAccepted;
}


double TriggerSelectionData::GetWeight(PECReader const &reader) const
{
    return (currentRange->PassEventSelection(reader)) ? 1. : 0.;
}


TriggerSelectionInterface *TriggerSelectionData::Clone() const
{
    return new TriggerSelectionData(ranges);
}



TriggerSelectionMC::TriggerSelectionMC(std::vector<TriggerRange const *> const &ranges_):
    TriggerSelectionInterface(),
    buffer(nullptr)
{
    // A sanity check
    if (ranges_.size() == 0)
        throw std::logic_error("TriggerSelectionMC::TriggerSelectionMC: The provided "
         "collection of pointers to TriggerRange objects is empty.");
    
    
    // Save pointers to the TriggerRange objects
    ranges.reserve(ranges_.size());
    
    for (auto const &r: ranges_)
        ranges.emplace_back(r, nullptr);

   TH1::AddDirectory(kFALSE);
    

    sf_file = new TFile("/afs/cern.ch/work/n/ntsirova/public/PECLocalFwk/data/Muon/muon_sf.root","read");
    
    //ROOTLock::Lock();
    
    // Delete the previous set of SFs
    //delete hist_id;
       
    // Read new SFs
    hist_id = dynamic_cast<TH1D *>(sf_file->Get("muon_id"));
    hist_iso = dynamic_cast<TH1D *>(sf_file->Get("muon_iso"));
    hist_trig = dynamic_cast<TH1D *>(sf_file->Get("muon_trig"));
}


TriggerSelectionMC::~TriggerSelectionMC()
{
    sf_file->Close();
    delete [] buffer;
}


void TriggerSelectionMC::UpdateTree(TTree *triggerTree_, bool)
{
    // Update the tree pointer and counters
    triggerTree = triggerTree_;
    nEntriesTree = triggerTree->GetEntries();
    nextEntryTree = 0;
    
    
    // Create buffers into which the trigger decision will be read. They are allocated in a single
    //block. The number of buffers cannot exceed ranges.size(), but it might be smaller if some
    //MC triggers are repeated in several TriggerRange objects. For the sake of simplicity and
    //performance, no attempt to count distinctive triggers is made, and an array of size
    //ranges.size() is allocated.
    buffer = new Bool_t[ranges.size()];
    
    
    // Set statuses and addresses of relevant branches of the trigger tree. A certain complication
    //is caused by the fact that the same trigger might be specified in several trigger ranges
    triggerTree->SetBranchStatus("*", false);
    unsigned curBufferIndex = 0;
    
    for (unsigned i = 0; i < ranges.size(); ++i)
    {
        string const &curTriggerName = ranges.at(i).first->GetMCTriggerPattern();
        
        
        // Try to find name of the current trigger in previous ranges
        unsigned iPrev = 0;
        
        for (; iPrev < i; ++iPrev)
            if (curTriggerName == ranges.at(iPrev).first->GetMCTriggerPattern())
                break;
        
        
        if (iPrev == i)
        //^ It means that the current trigger was not encountered before
        {
            // Get the corresponding branch of the tree
            TBranch *branch =
             triggerTree->GetBranch((curTriggerName + "__accept").c_str());
            
            if (not branch)
                throw runtime_error(string("TriggerSelectionMC::UpdateTree: State of the "
                 "trigger \"HLT_") + curTriggerName + "_v*\" is not stored in the source tree.");
            
            
            // Set the status and address of the branch
            branch->SetStatus(true);
            branch->SetAddress(buffer + curBufferIndex);
            
            // Store the address also in ranges
            ranges.at(i).second = buffer + curBufferIndex;
            
            
            // Increment the buffer counter
            ++curBufferIndex;
        }
        else
        //^ The trigger has already been encountered. The tree had been adjusted
        {
            // Only need to set the pointer to buffer in ranges
            ranges.at(i).second = ranges.at(iPrev).second;
        }
    }
}


bool TriggerSelectionMC::ReadNextEvent(EventID const &)
{
    // A sanity check
    if (not triggerTree)
        throw logic_error("TriggerSelectionMC::ReadNextEvent: Attempting to read an unspecified "
         "trigger tree.");
    
    
    // Check if there are still events to read
    if (nextEntryTree == nEntriesTree)
        return false;
    
    
    // Read the next entry
    triggerTree->GetEntry(nextEntryTree);
    ++nextEntryTree;
    
    return true;
}


bool TriggerSelectionMC::PassTrigger() const
{
    // Check all the requested triggers
    for (auto const &r: ranges)
        if (*r.second)  // the event is accepted by the MC trigger in this range
            return true;
    
    // The workflow reaches this point if only no trigger accepts the event
    return false;
}


double TriggerSelectionMC::GetWeight(PECReader const &reader) const
{
    double weight = 0.;
    
    for (auto const &r: ranges)
    {
        if (*r.second and r.first->PassEventSelection(reader))
            weight += r.first->GetLuminosity() * ScaleFactor(r.first, reader);
    }
    
    
    return weight;
}


double TriggerSelectionMC::ScaleFactor(TriggerRange const *, PECReader const &reader) const
{
    if (reader.GetLeptons().size() != 1)
        throw logic_error("There is not 1 lepton in the event.");
    auto const lepton = reader.GetLeptons().at(0);
    
    double sf_id, sf_iso, sf_trig, sf;
    sf_id = hist_id->GetBinContent(hist_id->FindFixBin(fabs(lepton.Eta())));
    sf_iso = hist_iso->GetBinContent(hist_iso->FindFixBin(fabs(lepton.Eta())));
    sf_trig = hist_trig->GetBinContent(hist_trig->FindFixBin(fabs(lepton.Eta())));
    sf = sf_id * sf_iso * sf_trig;
    
    //cout << lepton.Eta() << ' ' << sf_id << ' ' << sf_iso << ' ' << sf_trig << ' ' << sf << endl;
    return sf;
}


TriggerSelectionInterface *TriggerSelectionMC::Clone() const
{
    // Make a vector of pointers to the TriggerRange objects
    vector<TriggerRange const *> pureRanges;
    pureRanges.reserve(ranges.size());
    
    for (auto const &r: ranges)
        pureRanges.emplace_back(r.first);
    
    
    // Feed it to the appropriate constructor
    return new TriggerSelectionMC(pureRanges);
}



TriggerSelection::TriggerSelection(vector<TriggerRange const *> const &ranges_):
    TriggerSelectionInterface(),
    ranges(ranges_)
{}


void TriggerSelection::UpdateTree(TTree *triggerTree_, bool isData)
{
    // Create a new object to perform the trigger selection
    if (isData)
        selection.reset(new TriggerSelectionData(ranges));
    else
        selection.reset(new TriggerSelectionMC(ranges));
    
    
    // Set the trigger tree in the selection object
    selection->UpdateTree(triggerTree_, isData);
}


bool TriggerSelection::ReadNextEvent(EventID const &eventID)
{
    // A sanity check
    if (not selection)
        throw logic_error("TriggerSelection::ReadNextEvent: Attempting to read an unspecified "
         "trigger tree.");
    
    return selection->ReadNextEvent(eventID);
}


bool TriggerSelection::PassTrigger() const
{
    return selection->PassTrigger();
}


double TriggerSelection::GetWeight(PECReader const &reader) const
{
    return selection->GetWeight(reader);
}


TriggerSelectionInterface *TriggerSelection::Clone() const
{
    return new TriggerSelection(ranges);
}