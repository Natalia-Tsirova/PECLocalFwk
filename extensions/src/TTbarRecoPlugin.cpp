#include <TTbarRecoPlugin.hpp>

#include <Processor.hpp>

#include <deque>


using namespace std;


TTbarRecoPlugin::TTbarRecoPlugin(BTagger const &bTagger_):
    Plugin("TTbarReco"),
    bTagger(bTagger_)
{}


Plugin *TTbarRecoPlugin::Clone() const
{
    return new TTbarRecoPlugin(bTagger);
}


void TTbarRecoPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


void TTbarRecoPlugin::EndRun()
{}


bool TTbarRecoPlugin::ProcessEvent()
{
    // Make sure the event contains reasonable physics objects
    if ((*reader)->GetLeptons().size() not_eq 1 or
     (*reader)->GetJets().size() + (*reader)->GetAdditionalJets().size() < 4)
        return false;
    
    
    // Put all the reconstructed jets into a single vector
    allJets.clear();
    allJets.insert(allJets.end(), (*reader)->GetJets().begin(), (*reader)->GetJets().end());
    allJets.insert(allJets.end(), (*reader)->GetAdditionalJets().begin(),
     (*reader)->GetAdditionalJets().end());
    
    
    // Define some short-cuts
    auto const &lepton = (*reader)->GetLeptons().front();
    auto const &met = (*reader)->GetMET();
    
    
    // Precalculate variables that do not depend on the choice of event interpretation
    double Ht = lepton.Pt() + met.Pt();
    
    for (auto const &j: allJets)
        Ht += j.Pt();
    
    TLorentzVector p4RecoWLep(lepton.P4() + (*reader)->GetNeutrino().P4());
    
    
    // Variables to keep trace of the best event interpretation
    Interpretation bestInterpretation;
    double bestBNNScore = -100.;
    
    
    // There are four jets in semileptonic ttbar. Consider all the way to choose them without
    //ordering. The jets are chosen with the help of a mask
    deque<bool> mask(allJets.size(), false);
    fill(mask.end() - 4, mask.end(), true);
    
    // Loop over the different permutations of bits in the mask
    do
    {
        // Store the indices of the jets selected in the current combination
        unmaskedJetIndices.clear();
        
        for (unsigned i = 0; i < mask.size(); ++i)
            if (mask.at(i))
                unmaskedJetIndices.push_back(i);
        
        
        // For a given combination there are 4*3 = 12 ways to choose the b-jets; we do not
        //distinguish the two light-flavour jets
        for (unsigned i = 0; i < 4; ++i)  // chooses b-jet from the semileptonic top-quark decay
            for (unsigned j = 1; j < 4; ++j)  // chooses b-jet from the hadronic top-quark decay
            {
                // Define interpretations for jets given by unmaskedJetIndices. The variables
                //below are indices in unmaskedJetIndices vector
                unsigned const bTopLepIndex = i;
                unsigned const bTopHadIndex = (i + j) % 4;
                //^ It's easy to see that the above two numbers are always different
                
                unsigned const q1TopHadIndex = ((bTopLepIndex - bTopHadIndex + 4) % 4 >= 2) ?
                 (bTopHadIndex + 1) % 4 : (bTopHadIndex + 3) % 4;
                //^ Numbers modulus n can be thought of as segments of a circumference. Then
                //one can define an oriented distance d(a,b) as a number of segments that
                //should be added to a on the conterclockwise side in order to arrive to b.
                //It can be expressed as d(a,b) = (n + b - a) % n. We have two segments
                //occupied by bTopIndex and qIndex, so there are two more free segments, and
                //we need need to find index of one of them. Let d(a,b) >= d(b,a). Then the
                //first segment to anticlockwise direction from a is free. Its index is
                //given by (a + 1) % n. If d(a,b) < d(b,a), the first segment to clockwise
                //direction from a is free, which is (a + n - 1) % n.
                //This observation motivates the expression for b1HiggsIndex. One can also
                //rewrite it without a conditional operator:
                // (n + a + (b - a + n) % n * 2 / n * 2 - 1) % n,
                //but it looks even more cryptic
                
                unsigned const q2TopHadIndex = 6 - (bTopLepIndex + bTopHadIndex + q1TopHadIndex);
                //^ Since sum of all the four indices must be 0 + ... + 3 = 6
                
                
                // Fill the current interpretation
                Interpretation interpr;
                interpr.bTopLep = unmaskedJetIndices.at(bTopLepIndex);
                interpr.bTopHad = unmaskedJetIndices.at(bTopHadIndex);
                interpr.q1TopHad = unmaskedJetIndices.at(q1TopHadIndex);
                interpr.q2TopHad = unmaskedJetIndices.at(q2TopHadIndex);
                
                // Save 4-momenta of the reconstructed objects
                interpr.p4RecoTopLep = p4RecoWLep + allJets.at(interpr.bTopLep).P4();
                interpr.p4RecoWHad = allJets.at(interpr.q1TopHad).P4() +
                 allJets.at(interpr.q2TopHad).P4();
                interpr.p4RecoTopHad = allJets.at(interpr.bTopHad).P4() + interpr.p4RecoWHad;
                
                
                // Calculate input variables for the current interpretation
                CalculateRecoVars(interpr, lepton, p4RecoWLep, Ht);
                
                // Apply the BNN and check if the current interpretation is better than the best one
                //found so far (the list of input variables is provided in the header file with
                //description of the BNN)
                double const bnnScore = bnnReco(DeltaRBJetTopLepWLep, DeltaRLightJets,
                 log(MassTopHad), log(MassTopLep), log(MassWHad), MaxEtaTop, log(MinPtBJet),
                 PassBTagTopHad, PassBTagTopLep, RelHt);
                
                if (bnnScore > bestBNNScore)
                {
                    bestInterpretation = interpr;
                    bestBNNScore = bnnScore;
                }
            }
    }
    while (next_permutation(mask.begin(), mask.end()));
    
    
    // All the possible interpretations of the event have been checked, and the best one has been
    //found. Set the reconstructed four-moment so the user can read them
    recoTopQuarkLep.SetP4(bestInterpretation.p4RecoTopLep);
    recoTopQuarkHad.SetP4(bestInterpretation.p4RecoTopHad);
    recoWBosonHad.SetP4(bestInterpretation.p4RecoWHad);
    
    
    return true;
}


Candidate const &TTbarRecoPlugin::GetRecoTopQuarkLep() const
{
    return recoTopQuarkLep;
}


Candidate const &TTbarRecoPlugin::GetRecoTopQuarkHad() const
{
    return recoTopQuarkHad;
}


Candidate const &TTbarRecoPlugin::GetRecoWBosonHad() const
{
    return recoWBosonHad;
}


void TTbarRecoPlugin::CalculateRecoVars(Interpretation const &interpr, Lepton const &lepton,
 TLorentzVector const &p4RecoWLep, double Ht)
{
    MassTopLep = interpr.p4RecoTopLep.M();
    PtTopLep = interpr.p4RecoTopLep.Pt();
    EtaTopLep = interpr.p4RecoTopLep.Eta();
    
    MassTopHad = interpr.p4RecoTopHad.M();
    PtTopHad = interpr.p4RecoTopHad.Pt();
    EtaTopHad = interpr.p4RecoTopHad.Eta();
    
    MassWHad = interpr.p4RecoWHad.M();
    PtWHad = interpr.p4RecoWHad.Pt();
    EtaWHad = interpr.p4RecoWHad.Eta();
    
    DeltaRTopTop = interpr.p4RecoTopLep.DeltaR(interpr.p4RecoTopHad);
    DeltaRTopLepWLep = interpr.p4RecoTopLep.DeltaR(p4RecoWLep);
    DeltaRTopHadWHad = interpr.p4RecoTopHad.DeltaR(interpr.p4RecoWHad);
    DeltaRBJetTopLepWLep = p4RecoWLep.DeltaR(allJets.at(interpr.bTopLep).P4());
    DeltaRBJetTopHadWHad = interpr.p4RecoWHad.DeltaR(allJets.at(interpr.bTopHad).P4());
    DeltaRLightJets =
     allJets.at(interpr.q1TopHad).P4().DeltaR(allJets.at(interpr.q2TopHad).P4());
    
    MinEtaTop = min(fabs(interpr.p4RecoTopLep.Eta()), fabs(interpr.p4RecoTopHad.Eta()));
    MaxEtaTop = max(fabs(interpr.p4RecoTopLep.Eta()), fabs(interpr.p4RecoTopHad.Eta()));
    DEtaTopTop = fabs(interpr.p4RecoTopLep.Eta() - interpr.p4RecoTopHad.Eta());
    
    RelHt = (interpr.p4RecoTopLep.Pt() + interpr.p4RecoTopHad.Pt()) / Ht;
    
    
    PassBTagTopLep = (bTagger(allJets.at(interpr.bTopLep))) ? 1. : 0.;
    PassBTagTopHad = (bTagger(allJets.at(interpr.bTopHad))) ? 1. : 0.;
    NLightPassBTagTopHad = 0. + (bTagger(allJets.at(interpr.q1TopHad))) +
     (bTagger(allJets.at(interpr.q2TopHad)));
    
    
    MinPtBJet = min(allJets.at(interpr.bTopLep).Pt(), allJets.at(interpr.bTopHad).Pt());
    MinPtLightJet = min(allJets.at(interpr.q1TopHad).Pt(), allJets.at(interpr.q2TopHad).Pt());
    
    
    TLorentzVector p4Lep(lepton.P4());
    TLorentzVector p4Top(interpr.p4RecoTopLep);
    
    TVector3 b = p4RecoWLep.BoostVector();
    p4Lep.Boost(-b);
    p4Top.Boost(-b);
    
    CosLepTopLepWLep = p4Lep.Vect().Dot(p4Top.Vect()) / p4Lep.Vect().Mag() /
     p4Top.Vect().Mag();
}