#include <THRecoPlugin.hpp>

#include <Processor.hpp>

#include <deque>

#include <sys/stat.h>


using namespace std;


THRecoPlugin::THRecoPlugin(BTagger const &bTagger_):
    Plugin("THReco"),
    bTagger(bTagger_)
{}


Plugin *THRecoPlugin::Clone() const
{
    return new THRecoPlugin(bTagger);
}


void THRecoPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


void THRecoPlugin::EndRun()
{}


bool THRecoPlugin::ProcessEvent()
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
    
    TLorentzVector p4RecoW(lepton.P4() + (*reader)->GetNeutrino().P4());
    
    
    // Variables to keep trace of the best event interpretation
    Interpretation bestInterpretation;
    double bestBNNScore = -100.;
    
    
    // There are four jets in thq. Consider all the way to choose them without ordering. The jets
    //are chosen with the help of a mask
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
        
        
        // For a given combination there are 4*3 = 12 ways to choose the b-jet from
        //top-quark decay and the recoil jet; we do not distinguish the two b-jets from
        //Higgs boson decay
        for (unsigned i = 0; i < 4; ++i)  // chooses b-jet from the top-quark
            for (unsigned j = 1; j < 4; ++j)  // chooses jet from the recoil quark
            {
                // Define interpretations for jets given by unmaskedJetIndices. The variables
                //below are indices in unmaskedJetIndices vector
                unsigned const bTopIndex = i;
                unsigned const qIndex = (i + j) % 4;
                //^ It's easy to see that the above two numbers are always different
                
                unsigned const b1HiggsIndex = ((bTopIndex - qIndex + 4) % 4 >= 2) ?
                 (qIndex + 1) % 4 : (qIndex + 3) % 4;
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
                
                unsigned const b2HiggsIndex = 6 - (bTopIndex + qIndex + b1HiggsIndex);
                //^ Since sum of all the four indices must be 0 + ... + 3 = 6
                
                
                // Fill the current interpretation
                Interpretation interpr;
                interpr.bTop = unmaskedJetIndices.at(bTopIndex);
                interpr.qRecoil = unmaskedJetIndices.at(qIndex);
                interpr.b1Higgs = unmaskedJetIndices.at(b1HiggsIndex);
                interpr.b2Higgs = unmaskedJetIndices.at(b2HiggsIndex);
                
                // Save 4-momenta of the reconstructed objects
                interpr.p4RecoTop = p4RecoW + allJets.at(interpr.bTop).P4();
                interpr.p4RecoHiggs = allJets.at(interpr.b1Higgs).P4() +
                 allJets.at(interpr.b2Higgs).P4();
                
                
                // Calculate input variables for the current interpretation
                CalculateRecoVars(interpr, lepton, p4RecoW, Ht);
                
                // Apply the BNN and check if the current interpretation is better than the best one
                //found so far (the list of input variables is provided in the header file with
                //description of the BNN)
                double const bnnScore = bnnReco(fabs(EtaHiggs), fabs(EtaLJet), DeltaEtaLepTop,
                 DeltaRBJetsHiggs, log(MassHiggs), log(MassTop), log(MinPtBJet), NPassBTagHiggs,
                 PassBTagTop, RelHt);
                
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
    recoTopQuark.SetP4(bestInterpretation.p4RecoTop);
    recoHiggsBoson.SetP4(bestInterpretation.p4RecoHiggs);
    
    
    return true;
}


Candidate const &THRecoPlugin::GetRecoTopQuark() const
{
    return recoTopQuark;
}


Candidate const &THRecoPlugin::GetRecoHiggsBoson() const
{
    return recoHiggsBoson;
}


void THRecoPlugin::CalculateRecoVars(Interpretation const &interpr, Lepton const &lepton,
 TLorentzVector const &p4RecoW, double Ht)
{
    MassTop = interpr.p4RecoTop.M();
    PtTop = interpr.p4RecoTop.Pt();
    EtaTop = interpr.p4RecoTop.Eta();
    
    MassHiggs = interpr.p4RecoHiggs.M();
    PtHiggs = interpr.p4RecoHiggs.Pt();
    EtaHiggs = interpr.p4RecoHiggs.Eta();
    
    PtLJet = allJets.at(interpr.qRecoil).Pt();
    EtaLJet = allJets.at(interpr.qRecoil).Eta();
    
    DeltaRTopHiggs = interpr.p4RecoTop.DeltaR(interpr.p4RecoHiggs);
    
    DeltaRTopW = interpr.p4RecoTop.DeltaR(p4RecoW);
    DeltaRBJetTopW = p4RecoW.DeltaR(allJets.at(interpr.bTop).P4());
    
    DeltaEtaLepTop = fabs(lepton.Eta() - EtaTop);
    
    DeltaRBJetsHiggs =
     allJets.at(interpr.b1Higgs).P4().DeltaR(allJets.at(interpr.b2Higgs).P4());
    
    RelHt = (interpr.p4RecoTop.Pt() + interpr.p4RecoHiggs.Pt()) / Ht;
    
    MinPtBJet = min({allJets.at(interpr.bTop).Pt(), allJets.at(interpr.b1Higgs).Pt(),
     allJets.at(interpr.b2Higgs).Pt()});
        
    PassBTagTop = (bTagger(allJets.at(interpr.bTop))) ? 1. : 0.;
    NPassBTagHiggs = 0. + (bTagger(allJets.at(interpr.b1Higgs))) +
     (bTagger(allJets.at(interpr.b2Higgs)));
}