#include <UntagEventSelection.hpp>

#include <algorithm>


using namespace std;


UntagEventSelection::JetBin::JetBin(unsigned nJets_):
    nJets(nJets_)
{}


UntagEventSelection::UntagEventSelection(double jetPtThreshold_):
    EventSelectionInterface(),
    jetPtThreshold(jetPtThreshold_)
{
    // Fill the lepton-flavour map
    flavourMap[Lepton::Flavour::Electron] = 0;
    flavourMap[Lepton::Flavour::Muon] = 1;
    flavourMap[Lepton::Flavour::Tau] = 2;
}


bool UntagEventSelection::PassLeptonStep(vector<Lepton> const &tightLeptons,
 vector<Lepton> const &looseLeptons) const
{
    // Both the tight leptons collection and the thresholds of each flavour are sorted in the
    //decreasing order in pt. The algorithm checks that nth lepton of a given flavour has a greater
    //pt than the nth threshold for this flavour. If there are more leptons than thresholds, all the
    //additional leptons must fail the softest threshold. But the method also vetoes the additional
    //loose leptons, therefore one can just reject an event which contains more leptons than there
    //are thresholds specified. On the other hand, if the number of leptons is smaller than the
    //number of thresholds, the event lacks tight leptons and must be rejected, too.
    
    // Initialize the iterators for all the lepton flavours
    for (unsigned i = 0; i < 3; ++i)
        leptonThresholdIts.at(i) = leptonPtThresholds.at(i).cbegin();
    
    
    // Loop over the tight leptons
    for (auto const &lep: tightLeptons)
    {
        // Choose the iterator corresponding to the flavour of the current lepton
        unsigned const flavourIndex = flavourMap.at(lep.GetFlavour());
        auto &thresholdIt = leptonThresholdIts.at(flavourIndex);
        
        // Make sure we have not yet exhausted the allowed number of tight leptons of such flavour
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


bool UntagEventSelection::PassJetStep(vector<Jet> const &jets) const
{
    // Calculate the jet and the b-tagged jet multiplicities
    unsigned nJets = 0;
    
    for (auto const &j: jets)
        //if (IsAnalysisJet(j))  // the input collection is required to be already filtered
        {
            ++nJets;
            
        }
    
    
    // Check against the allowed jet bins
    for (auto const &bin: jetBins)
        if (bin.nJets == nJets)
            return true;
    //^ The brute-force looping over the allowed bins is employed as their number is small, and the
    //effect of an overhead from an optimization will probably be bigger than the possible gain in
    //performance
    return false;
}


bool UntagEventSelection::IsAnalysisJet(Jet const &jet) const
{
    return (jet.Pt() > jetPtThreshold);
}


void UntagEventSelection::AddLeptonThreshold(Lepton::Flavour flavour, double ptThreshold)
{
    unsigned const i = flavourMap.at(flavour);
    
    // Find the first element that is smaller than the new threshold. It might be end().
    auto const it = find_if(leptonPtThresholds.at(i).begin(), leptonPtThresholds.at(i).end(),
     [ptThreshold](double pt){return (pt < ptThreshold);});
    
    // Insert the new threshold just before the found element. Thanks to this, the list is always
    //sorted in the decreasing order
    leptonPtThresholds.at(i).insert(it, ptThreshold);
}


void UntagEventSelection::AddJetBin(unsigned nJets)
{
    jetBins.emplace_back(nJets);
}


EventSelectionInterface *UntagEventSelection::Clone() const
{
    return new UntagEventSelection(*this);
}
