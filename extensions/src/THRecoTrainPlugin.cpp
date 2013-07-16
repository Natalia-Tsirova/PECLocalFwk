#include <THRecoTrainPlugin.hpp>

#include <Processor.hpp>
#include <ROOTLock.hpp>

#include <deque>

#include <sys/stat.h>


using namespace std;


THRecoTrainPlugin::THRecoTrainPlugin(string const &outDirectory_, BTagger const &bTagger_,
 bool pruned_ /*= true*/):
    Plugin("THRecoTrain"),
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


Plugin *THRecoTrainPlugin::Clone() const
{
    return new THRecoTrainPlugin(outDirectory, bTagger, pruned);
}


void THRecoTrainPlugin::BeginRun(Dataset const &dataset)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    
    
    // Creation of ROOT objects is not thread-safe and must be protected
    ROOTLock::Lock();
    
    // Create the output file
    file = new TFile((outDirectory + dataset.GetFiles().front().GetBaseName() + ".root").c_str(),
     "recreate");
    
    // Create the tree
    tree = new TTree("Vars", "Basic kinematical variables");
    
    // End of critical block
    ROOTLock::Unlock();
    
    
    // Assign branch addresses
    tree->Branch("run", &runNumber);
    tree->Branch("event", &eventNumber);
    tree->Branch("lumiSection", &lumiSection);
    
    tree->Branch("InterpretationRank", &InterpretationRank);
    tree->Branch("Distance", &Distance);
    
    tree->Branch("MassTop", &MassTop);
    tree->Branch("PtTop", &PtTop);
    tree->Branch("EtaTop", &EtaTop);
    
    tree->Branch("MassHiggs", &MassHiggs);
    tree->Branch("PtHiggs", &PtHiggs);
    tree->Branch("EtaHiggs", &EtaHiggs);
    
    tree->Branch("PtLJet", &PtLJet);
    tree->Branch("EtaLJet", &EtaLJet);
    
    tree->Branch("DeltaRTopHiggs", &DeltaRTopHiggs);
    tree->Branch("DeltaRTopW", &DeltaRTopW);
    tree->Branch("DeltaRBJetTopW", &DeltaRBJetTopW);
    tree->Branch("DeltaEtaLepTop", &DeltaEtaLepTop);
    tree->Branch("DeltaRBJetsHiggs", &DeltaRBJetsHiggs);
    
    tree->Branch("CosLepTopW", &CosLepTopW);
    
    tree->Branch("RelHt", &RelHt);
    
    tree->Branch("MinPtBJet", &MinPtBJet);
    
    tree->Branch("PassBTagTop", &PassBTagTop);
    tree->Branch("PassBTagLJet", &PassBTagLJet);
    tree->Branch("NPassBTagHiggs", &NPassBTagHiggs);
    
    tree->Branch("CSVBJetTop", &CSVBJetTop);
    tree->Branch("CSVLJet", &CSVLJet);
    tree->Branch("MinCSVBJetsHiggs", &MinCSVBJetsHiggs);
    
    if (dataset.IsMC())
        tree->Branch("weight", &weight);
}


void THRecoTrainPlugin::EndRun()
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


bool THRecoTrainPlugin::ProcessEvent()
{
    // Make sure the event contains reasonable physics objects
    if ((*reader)->GetLeptons().size() not_eq 1 or (*reader)->GetJets().size() < 2)
        return false;
    
    
    // Find four-momenta of relevant generator-level particles
    TLorentzVector p4GenTop, p4GenHiggs, p4GenRecoil;
    
    for (auto const &p: (*reader)->GetHardGenParticles())
    {
        int const absPdgId = abs(p.GetPdgId());
        
        if (absPdgId == 6)  // a top-quark
            p4GenTop = p.P4();
        else if (absPdgId == 25)  // a Higgs boson
            p4GenHiggs = p.P4();
        else if (absPdgId <= 4)  // a candidate for the recoil quark
        {
            // The recoil quark is the only light-flavour quark (though, charm is also considered as
            //light) in the final state. A particle in the final state can be identified by a
            //requirement of having at least one great-grandmother. This condition is checked
            //below
            
            list<GenParticle const *> const *mothers = &p.GetMothers();
            
            if (mothers->size() > 0)
            {
                mothers = &mothers->front()->GetMothers();
                
                if (mothers->size() > 0)
                    p4GenRecoil = p.P4();
            }
        }
    }
    
    // A sanity check
    if (p4GenTop.Pt() == 0. or p4GenHiggs.Pt() == 0. or p4GenRecoil.Pt() == 0.)
        throw runtime_error("THRecoTrainPlugin::ProcessEvent: One of the required generator-level "
         "particles has not been found.");
    
    
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
    
    weight = (*reader)->GetCentralWeight();
    
    double Ht = lepton.Pt() + met.Pt();
    
    for (auto const &j: allJets)
        Ht += j.Pt();
    
    TLorentzVector p4RecoW(lepton.P4() + (*reader)->GetNeutrino().P4());
    
    
    // Clear a vector that will store all interpretations of the current event
    interpretations.clear();
    
    
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
                
                // Calculate distance to generator-level objets ($\Delta R \oplus \Delta p_T^{rel}$
                //metrics is used)
                interpr.distance = interpr.p4RecoTop.DeltaR(p4GenTop) +
                 interpr.p4RecoHiggs.DeltaR(p4GenHiggs) +
                 allJets.at(interpr.qRecoil).P4().DeltaR(p4GenRecoil) +
                 fabs(interpr.p4RecoTop.Pt() - p4GenTop.Pt()) / p4GenTop.Pt() +
                 fabs(interpr.p4RecoHiggs.Pt() - p4GenHiggs.Pt()) / p4GenHiggs.Pt() +
                 fabs(allJets.at(interpr.qRecoil).Pt() - p4GenRecoil.Pt()) /
                  p4GenRecoil.Pt();
                
                
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
        CalculateRecoVars(interpretations.at(interprIndex), lepton, p4RecoW, Ht);
        
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
            CalculateRecoVars(interpretations.at(interprIndex), lepton, p4RecoW, Ht);
            
            tree->Fill();
        }
    }
    
    
    return true;
}


void THRecoTrainPlugin::CalculateRecoVars(Interpretation const &interpr, Lepton const &lepton,
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
    PassBTagLJet = (bTagger(allJets.at(interpr.qRecoil))) ? 1. : 0.;
    NPassBTagHiggs = 0. + (bTagger(allJets.at(interpr.b1Higgs))) +
     (bTagger(allJets.at(interpr.b2Higgs)));
    
    CSVBJetTop = max(0., allJets.at(interpr.bTop).CSV());
    CSVLJet = max(0., allJets.at(interpr.qRecoil).CSV());
    MinCSVBJetsHiggs = max(0.,
     min(allJets.at(interpr.b1Higgs).CSV(), allJets.at(interpr.b2Higgs).CSV()));
    
    
    TLorentzVector p4Lep(lepton.P4());
    TLorentzVector p4Top(interpr.p4RecoTop);
    
    TVector3 b = p4RecoW.BoostVector();
    p4Lep.Boost(-b);
    p4Top.Boost(-b);
    
    CosLepTopW = -p4Lep.Vect().Dot(p4Top.Vect()) / p4Lep.Vect().Mag() /
     p4Top.Vect().Mag();
}