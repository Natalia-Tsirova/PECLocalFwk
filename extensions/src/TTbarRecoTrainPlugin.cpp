#include <TTbarRecoTrainPlugin.hpp>

#include <Processor.hpp>
#include <ROOTLock.hpp>

#include <deque>

#include <sys/stat.h>


using namespace std;


TTbarRecoTrainPlugin::TTbarRecoTrainPlugin(string const &outDirectory_, BTagger const &bTagger_,
 bool pruned_ /*= true*/):
    Plugin("TTbarRecoTrain"),
    bTagger(bTagger_), outDirectory(outDirectory_), pruned(pruned_)
{
    // Make sure the directory path ends with a slash
    if (outDirectory.back() != '/')
        outDirectory += '/';
    
    // Create the output directory if it does not exist
    struct stat dirStat;
    
    if (stat(outDirectory.c_str(), &dirStat) != 0)  // the directory does not exist
        mkdir(outDirectory.c_str(), 0755);
    
    
    // Initialise random-number generator if needed
    if (pruned)
        rGen.reset(new TRandom3(0));
}


Plugin *TTbarRecoTrainPlugin::Clone() const
{
    return new TTbarRecoTrainPlugin(outDirectory, bTagger, pruned);
}


void TTbarRecoTrainPlugin::BeginRun(Dataset const &dataset)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    
    
    // Creation of ROOT objects is not thread-safe and must be protected
    ROOTLock::Lock();
    
    // Create the output file
    file = new TFile((outDirectory + dataset.GetFiles().front().GetBaseName() + ".root").c_str(),
     "recreate");
    
    // Create the tree
    tree = new TTree("Vars", "Observables for ttbar MVA reconstruction");
    
    // End of critical block
    ROOTLock::Unlock();
    
    
    // Assign branch addresses
    tree->Branch("run", &runNumber);
    tree->Branch("event", &eventNumber);
    tree->Branch("lumiSection", &lumiSection);
    
    tree->Branch("NJets30", &NJets30);
    tree->Branch("NTags30", &NTags30);
    
    tree->Branch("InterpretationRank", &InterpretationRank);
    tree->Branch("Distance", &Distance);
    
    tree->Branch("MassTopLep", &MassTopLep);
    tree->Branch("PtTopLep", &PtTopLep);
    tree->Branch("EtaTopLep", &EtaTopLep);
    
    tree->Branch("MassTopHad", &MassTopHad);
    tree->Branch("PtTopHad", &PtTopHad);
    tree->Branch("EtaTopHad", &EtaTopHad);
    
    tree->Branch("MassWHad", &MassWHad);
    tree->Branch("PtWHad", &PtWHad);
    tree->Branch("EtaWHad", &EtaWHad);
    
    tree->Branch("DeltaRTopTop", &DeltaRTopTop);
    
    tree->Branch("DeltaRTopLepWLep", &DeltaRTopLepWLep);
    tree->Branch("DeltaRTopHadWHad", &DeltaRTopHadWHad);
    tree->Branch("DeltaRBJetTopLepWLep", &DeltaRBJetTopLepWLep);
    tree->Branch("DeltaRBJetTopHadWHad", &DeltaRBJetTopHadWHad);
    
    tree->Branch("DeltaRLightJets", &DeltaRLightJets);
    
    tree->Branch("MinEtaTop", &MinEtaTop);
    tree->Branch("MaxEtaTop", &MaxEtaTop);
    tree->Branch("DEtaTopTop", &DEtaTopTop);
    
    tree->Branch("CosLepTopLepWLep", &CosLepTopLepWLep);
    
    tree->Branch("RelHt", &RelHt);
    
    tree->Branch("MinPtBJet", &MinPtBJet);
    tree->Branch("MinPtLightJet", &MinPtLightJet);
    
    tree->Branch("PassBTagTopLep", &PassBTagTopLep);
    tree->Branch("PassBTagTopHad", &PassBTagTopHad);
    
    tree->Branch("NLightPassBTagTopHad", &NLightPassBTagTopHad);
    
    tree->Branch("CSVBJetTopLep", &CSVBJetTopLep);
    tree->Branch("CSVBJetTopHad", &CSVBJetTopHad);
    tree->Branch("MaxCSVLightJetsTopHad", &MaxCSVLightJetsTopHad);
    
    if (dataset.IsMC())
        tree->Branch("weight", &weight);
}


void TTbarRecoTrainPlugin::EndRun()
{
    // Operations with ROOT objects performed here are not thread-safe and must be guarded
    ROOTLock::Lock();
    
    // Write the tree and close the file
    file->cd();
    tree->Write("", TObject::kOverwrite);
    
    // Delete the objects
    delete tree;
    delete file;
    
    ROOTLock::Unlock();
}


bool TTbarRecoTrainPlugin::ProcessEvent()
{
    // Make sure the event contains reasonable physics objects
    if ((*reader)->GetLeptons().size() not_eq 1 or
     (*reader)->GetJets().size() + (*reader)->GetAdditionalJets().size() < 4)
        return false;
    
    
    // Find four-momenta of relevant generator-level particles
    TLorentzVector p4GenTopLep, p4GenTopHad, p4GenWHad;
    
    // Loop over the generator particles
    for (auto const &p: (*reader)->GetHardGenParticles())
    {
        int const absPdgId = abs(p.GetPdgId());
        
        if (absPdgId == 6)  // a top-quark
        {
            if (p.FindFirstDaughterRecursive({11, 13, 15, -11, -13, -15}) not_eq nullptr)
                p4GenTopLep = p.P4();
            else
                p4GenTopHad = p.P4();
        }
        else if (absPdgId == 24)  // a W-boson
        {
            if (p.FindFirstDaughterRecursive({11, 13, 15, -11, -13, -15}) == nullptr)
                p4GenWHad = p.P4();
        }
    }
    
    // A sanity check
    if (p4GenTopLep.Pt() == 0. or p4GenTopHad.Pt() == 0. or p4GenWHad.Pt() == 0.)
        throw runtime_error("TTbarRecoTrainPlugin::ProcessEvent: One of the required "
         "generator-level particles has not been found.");
    
    
    // Put all the reconstructed jets into a vector
    allJets.clear();
    allJets.insert(allJets.end(), (*reader)->GetJets().begin(), (*reader)->GetJets().end());
    allJets.insert(allJets.end(), (*reader)->GetAdditionalJets().begin(),
     (*reader)->GetAdditionalJets().end());
    
    
    // Define some short-cuts
    auto const &lepton = (*reader)->GetLeptons().front();
    auto const &met = (*reader)->GetMET();
    
    
    // Precalculate variables that do not depend on the choice of event interpretation
    auto const &eventID = (*reader)->GetEventID();
    runNumber = eventID.Run();
    eventNumber = eventID.Event();
    lumiSection = eventID.LumiBlock();
    
    NJets30 = (*reader)->GetJets().size();
    NTags30 = 0;
    
    for (auto const &j: (*reader)->GetJets())
        if (bTagger(j))
            ++NTags30;
    
    weight = (*reader)->GetCentralWeight();
    
    double Ht = lepton.Pt() + met.Pt();
    
    for (auto const &j: allJets)
        Ht += j.Pt();
    
    TLorentzVector p4RecoWLep(lepton.P4() + (*reader)->GetNeutrino().P4());
    
    
    // Clear a vector that will store all interpretations of the current event
    interpretations.clear();
    
    
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
                
                // Calculate distance to generator-level objets ($\Delta R \oplus \Delta p_T^{rel}$
                //metrics is used)
                interpr.distance = interpr.p4RecoTopLep.DeltaR(p4GenTopLep) +
                 interpr.p4RecoTopHad.DeltaR(p4GenTopHad) +
                 interpr.p4RecoWHad.DeltaR(p4GenWHad) +
                 fabs(interpr.p4RecoTopLep.Pt() - p4GenTopLep.Pt()) / p4GenTopLep.Pt() +
                 fabs(interpr.p4RecoTopHad.Pt() - p4GenTopHad.Pt()) / p4GenTopHad.Pt() +
                 fabs(interpr.p4RecoWHad.Pt() - p4GenWHad.Pt()) / p4GenWHad.Pt();
                
                
                // Store the interpretation
                interpretations.push_back(interpr);
            }
    }
    while (next_permutation(mask.begin(), mask.end()));
    
    
    // All the interpretations of the current event have been found. Now sort them in the distance
    sort(interpretations.begin(), interpretations.end(),
     [](Interpretation const &lh, Interpretation const &rh){return (rh.distance > lh.distance);});
    
    
    // Calculate input variables
    if (pruned)  // only one interpretation will be saved
    {
        unsigned interprIndex;
        
        if (rGen->Integer(2) == 0)  // 50% chance to save the best interpretation has been taken
        {
            InterpretationRank = 0;
            interprIndex = 0;
        }
        else  // 50% chance to take a background interpretation has been taken
        {
            // Choose one non-best interpretation randomly
            interprIndex = 1 + rGen->Integer(interpretations.size() - 1);
            
            if (interprIndex == interpretations.size() - 1)  // the worst interpretation
                InterpretationRank = 2;
            else
                InterpretationRank = 1;
        }
        
        Distance = interpretations.at(interprIndex).distance;
        CalculateRecoVars(interpretations.at(interprIndex), lepton, p4RecoWLep, Ht);
        
        tree->Fill();
    }
    else  // save all the interpretations
    {
        for (unsigned interprIndex = 0; interprIndex < interpretations.size(); ++interprIndex)
        {
            if (interprIndex == 0)  // the best interpretation
                InterpretationRank = 0;
            else if (interprIndex == interpretations.size() - 1)  // the worst interpretation
                InterpretationRank = 2;
            else
                InterpretationRank = 1;
            
            Distance = interpretations.at(interprIndex).distance;
            CalculateRecoVars(interpretations.at(interprIndex), lepton, p4RecoWLep, Ht);
            
            tree->Fill();
        }
    }
    
    
    return true;
}


void TTbarRecoTrainPlugin::CalculateRecoVars(Interpretation const &interpr, Lepton const &lepton,
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
    
    CSVBJetTopLep = max(0., allJets.at(interpr.bTopLep).CSV());
    CSVBJetTopHad = max(0., allJets.at(interpr.bTopHad).CSV());
    MaxCSVLightJetsTopHad = max({0., allJets.at(interpr.q1TopHad).CSV(),
     allJets.at(interpr.q2TopHad).CSV()});
    
    
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