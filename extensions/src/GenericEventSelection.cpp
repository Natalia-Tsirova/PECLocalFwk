#include <GenericEventSelection.hpp>

#include <algorithm>


using namespace std;


GenericEventSelection::JetTagBin::JetTagBin(unsigned nJets_, unsigned nTags_):
    nJets(nJets_), nTags(nTags_)
{}


GenericEventSelection::GenericEventSelection(double jetPtThreshold_, BTagger const &bTagger_):
    EventSelectionInterface(),
    jetPtThreshold(jetPtThreshold_), bTagger(bTagger_)
{
    // Fill the lepton-flavour map
    flavourMap[Lepton::Flavour::Electron] = 0;
    flavourMap[Lepton::Flavour::Muon] = 1;
    flavourMap[Lepton::Flavour::Tau] = 2;
}


bool GenericEventSelection::PassLeptonStep(vector<Lepton> const &tightLeptons,
 vector<Lepton> const &looseLeptons) const
{
    // Both the tight leptons collection and the thresholds of each flavour are sorted in the
    //decreasing order in pt. The algorithm checks that nth lepton of a given flavour has a greater
    //pt than the nth threshold for this flavour. If there are more leptons than thresholds, all the
    //additional leptons must fail the softest threshold. But the method also vetoes the additional
    //loose leptons, therefore one can just reject an event which contains more leptons than there
    //are thresholds specified. This is beacause tight leptons are a strict subset of loose ones. On
    //the other hand, if the number of leptons is smaller than the number of thresholds, the event
    //lacks tight leptons and must be rejected, too.
    
    // Initialize the iterators for all the lepton flavours
    for (unsigned i = 0; i < 3; ++i)
        leptonThresholdIts.at(i) = leptonPtThresholds.at(i).cbegin();
    
    
    // Loop over the tight leptons
    for (auto const &lep: tightLeptons)
    {
        // Choose the iterator corresponding to the flavour of the current lepton
        unsigned const flavourIndex = flavourMap.at(lep.GetFlavour());
        auto &thresholdIt = leptonThresholdIts.at(flavourIndex);
        
        // Make sure we have not yet exhausted the allowed number of tight leptons of such flavour.
        //However, remember that PECReader fills the collection of tight leptons with the same pt
        //treshold as for loose leptons. For this reason the condition below can be satisfied either
        //because there are, indeed, too many high-pt leptons or because there are additional
        //leptons that should be vetoed. In any case the event should be rejected
        if (thresholdIt == leptonPtThresholds.at(flavourIndex).cend())
            return false;
        
        // Compare the lepton's pt to the threshold. For the valid logic of this comparison both the
        //tight leptons and the thresholds must be sorted in the decreasing order in pt
        if (lep.Pt() < *thresholdIt)
            return false;
        
        // Move this iterator to the next (lower) threshold
        ++thresholdIt;
    }
    
    
    // Make sure all the iterators have reached the ends of the corresponding lists. If it is not
    //the case, the event contains fewer leptons than requested
    for (unsigned i = 0; i < 3; ++i)
        if (leptonThresholdIts.at(i) != leptonPtThresholds.at(i).cend())
            return false;
    
    
    // Veto the additional loose leptons. Since they are required to include the tight leptons, it
    //is sufficient to simply check their number
    if (tightLeptons.size() != looseLeptons.size())
        return false;
    
    
    // If the algorithm has reached this point, the event-selection requirements are satisfied
    return true;
}


bool GenericEventSelection::PassJetStep(vector<Jet> const &jets) const
{
    // Calculate the jet and the b-tagged jet multiplicities
    unsigned nJets = 0, nTags = 0;
    
    for (auto const &j: jets)
        //if (IsAnalysisJet(j))  // the input collection is required to be already filtered
        {
            ++nJets;
            
            if (bTagger(j))
                ++nTags;
        }
    
    
    // Check against the allowed jet-tag bins
    for (auto const &bin: jetBins)
        if (bin.nJets == nJets and bin.nTags == nTags)
            return true;
    //^ The brute-force looping over the allowed bins is employed as their number is small, and the
    //effect of an overhead from an optimization will probably be bigger than the possible gain in
    //performance
    return false;
}


bool GenericEventSelection::IsAnalysisJet(Jet const &jet) const
{
    return (jet.Pt() > jetPtThreshold);
}


void GenericEventSelection::AddLeptonThreshold(Lepton::Flavour flavour, double ptThreshold)
{
    unsigned const i = flavourMap.at(flavour);
    
    // Find the first element that is smaller than the new threshold. It might be end().
    auto const it = find_if(leptonPtThresholds.at(i).begin(), leptonPtThresholds.at(i).end(),
     [ptThreshold](double pt){return (pt < ptThreshold);});
    
    // Insert the new threshold just before the found element. Thanks to this, the list is always
    //sorted in the decreasing order
    leptonPtThresholds.at(i).insert(it, ptThreshold);
}


void GenericEventSelection::AddJetTagBin(unsigned nJets, unsigned nTags)
{
    jetBins.emplace_back(nJets, nTags);
}


EventSelectionInterface *GenericEventSelection::Clone() const
{
    return new GenericEventSelection(*this);
}
