#include <THInputVarsPlugin.hpp>

#include <Processor.hpp>
#include <ROOTLock.hpp>

#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <deque>

#include <sys/stat.h>


using namespace std;


THInputVarsPlugin::THInputVarsPlugin(string const &outDirectory_, BTagger const &bTagger_):
    Plugin("THInputVars"),
    bTagger(bTagger_), outDirectory(outDirectory_)
{
    // Make sure the directory path ends with a slash
    if (outDirectory.back() != '/')
        outDirectory += '/';
    
    // Create the output directory if it does not exist
    struct stat dirStat;
    
    if (stat(outDirectory.c_str(), &dirStat) != 0)  // the directory does not exist
        mkdir(outDirectory.c_str(), 0755);
}


Plugin *THInputVarsPlugin::Clone() const
{
    return new THInputVarsPlugin(outDirectory, bTagger);
}


void THInputVarsPlugin::BeginRun(Dataset const &dataset)
{
    // Save pointers to other required plugins
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    thqReconstructor =
     dynamic_cast<THRecoPlugin const *>(processor->GetPluginBefore("THReco", name));
    ttbarReconstructor =
     dynamic_cast<TTbarRecoPlugin const *>(processor->GetPluginBefore("TTbarReco", name));
    
    
    // Creation of ROOT objects is not thread-safe and must be protected
    ROOTLock::Lock();
    
    // Create the output file
    file = new TFile((outDirectory + dataset.GetFiles().front().GetBaseName() + ".root").c_str(),
     "recreate");
    
    // Create the tree
    tree = new TTree("Vars", "Observables for thq extraction");
    
    // End of critical block
    ROOTLock::Unlock();
    
    
    // Assign branch addresses
    tree->Branch("run", &runNumber);
    tree->Branch("event", &eventNumber);
    tree->Branch("lumiSection", &lumiSection);
    
    tree->Branch("NJets30", &NJets30);
    tree->Branch("NTags30", &NTags30);
    
    tree->Branch("thq_MassHiggs", &thq_MassHiggs);
    tree->Branch("thq_PtHiggs", &thq_PtHiggs);
    tree->Branch("thq_EtaHiggs", &thq_EtaHiggs);
    
    tree->Branch("thq_PtLJet", &thq_PtLJet);
    tree->Branch("thq_EtaLJet", &thq_EtaLJet);
    
    tree->Branch("thq_DeltaRTopHiggs", &thq_DeltaRTopHiggs);
    tree->Branch("thq_DeltaRBJetsHiggs", &thq_DeltaRBJetsHiggs);
    
    tree->Branch("thq_CosLepLJetTH", &thq_CosLepLJetTH);
    
    tree->Branch("thq_MassTopHiggs", &thq_MassTopHiggs);
    
    tree->Branch("tt_MassTopHad", &tt_MassTopHad);
    tree->Branch("tt_PtTopHad", &tt_PtTopHad);
    tree->Branch("tt_EtaTopHad", &tt_EtaTopHad);
    
    tree->Branch("tt_MassWHad", &tt_MassWHad);
    tree->Branch("tt_PtWHad", &tt_PtWHad);
    tree->Branch("tt_EtaWHad", &tt_EtaWHad);
    
    tree->Branch("tt_RelHt", &tt_RelHt);
    tree->Branch("tt_DeltaRLightJets", &tt_DeltaRLightJets);
    tree->Branch("tt_MaxMassBHadQ", &tt_MaxMassBHadQ);
    
    tree->Branch("glb_PtJ1", &glb_PtJ1);
    tree->Branch("glb_PtJ2", &glb_PtJ2);
    
    tree->Branch("glb_SqrtSHat", &glb_SqrtSHat);
    tree->Branch("glb_Sphericity", &glb_Sphericity);
    
    if (dataset.IsMC())
        tree->Branch("weight", &weight);
}


void THInputVarsPlugin::EndRun()
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


bool THInputVarsPlugin::ProcessEvent()
{
    // Put all the reconstructed jets into a single vector
    allJets.clear();
    allJets.insert(allJets.end(), (*reader)->GetJets().begin(), (*reader)->GetJets().end());
    allJets.insert(allJets.end(), (*reader)->GetAdditionalJets().begin(),
     (*reader)->GetAdditionalJets().end());
    
    
    // Define some short-cuts
    auto const &lepton = (*reader)->GetLeptons().front();
    auto const &met = (*reader)->GetMET();
    
    
    // Calculate service and global variables
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
    
    TLorentzVector p4RecoW(lepton.P4() + (*reader)->GetNeutrino().P4());
    
    
    glb_PtJ1 = allJets.at(0).Pt();
    glb_PtJ2 = allJets.at(1).Pt();
    
    TLorentzVector p4AllJets;
    
    for (auto const &j: allJets)
        p4AllJets += j.P4();
    
    glb_SqrtSHat = (p4AllJets + p4RecoW).M();
    
    
    // Calculate sphericity
    TMatrixDSym sphericityTensor(3);
    double norm = 0.;
    
    for (auto const &p: {lepton.P4(), (*reader)->GetNeutrino().P4()})
    {
        TVector3 p3(p.Vect());
        norm += p3.Mag2();
        
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                sphericityTensor(i, j) += p3[i] * p3[j];
    }
    
    for (auto const &jet: allJets)
    {
        TVector3 p3(jet.P4().Vect());
        norm += p3.Mag2();
        
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                sphericityTensor(i, j) += p3[i] * p3[j];
    }
    
    sphericityTensor *= 1. / norm;
    
    TMatrixDSymEigen eigenValCalc(sphericityTensor);
    TVectorD eigenVals(eigenValCalc.GetEigenValues());
    
    glb_Sphericity = 1.5 * (eigenVals[1] + eigenVals[2]);
    
    
    // Calculate variables reconstructed under thq hypothesis
    auto const &higgs = thqReconstructor->GetRecoHiggsBoson();
    auto const &top = thqReconstructor->GetRecoTopQuark();
    auto const &recoilQuark = thqReconstructor->GetRecoRecoilQuark();
    
    thq_MassHiggs = higgs.M();
    thq_PtHiggs = higgs.Pt();
    thq_EtaHiggs = higgs.Eta();
    
    thq_PtLJet = recoilQuark.Pt();
    thq_EtaLJet = recoilQuark.Eta();
    
    thq_DeltaRTopHiggs = higgs.P4().DeltaR(top.P4());
    thq_DeltaRBJetsHiggs = allJets.at(thqReconstructor->GetInterpretation().b1Higgs).P4().DeltaR(
     allJets.at(thqReconstructor->GetInterpretation().b2Higgs).P4());
    
    thq_MassTopHiggs = (higgs.P4() + top.P4()).M();
    
    
    // Calculate thq_CosLepLJetTH
    TVector3 b((higgs.P4() + top.P4()).BoostVector());
    
    TLorentzVector boostedLepton(lepton.P4());
    boostedLepton.Boost(-b);
    TVector3 const p3Lepton(boostedLepton.Vect());
    
    TLorentzVector boostedLJet(recoilQuark.P4());
    boostedLJet.Boost(-b);
    TVector3 const p3LJet(boostedLJet.Vect());
    
    thq_CosLepLJetTH = p3Lepton.Dot(p3LJet) / (p3Lepton.Mag() * p3LJet.Mag());
    
    
    // Finally calculate observables constructed under ttbar hypothesis
    auto const &topHad = ttbarReconstructor->GetRecoTopQuarkHad();
    auto const &wHad = ttbarReconstructor->GetRecoWBosonHad();
    
    tt_MassTopHad = topHad.M();
    tt_PtTopHad = topHad.Pt();
    tt_EtaTopHad = topHad.Eta();
    
    tt_MassWHad = wHad.M();
    tt_PtWHad = wHad.Pt();
    tt_EtaWHad = wHad.Eta();
    
    tt_RelHt = (topHad.Pt() + ttbarReconstructor->GetRecoTopQuarkLep().Pt()) / Ht;
    
    
    auto const &q1TopHad = allJets.at(ttbarReconstructor->GetInterpretation().q1TopHad);
    auto const &q2TopHad = allJets.at(ttbarReconstructor->GetInterpretation().q2TopHad);
    auto const &bTopHad = allJets.at(ttbarReconstructor->GetInterpretation().bTopHad);
    
    tt_DeltaRLightJets = q1TopHad.P4().DeltaR(q2TopHad.P4());
    tt_MaxMassBHadQ = max((bTopHad.P4() + q1TopHad.P4()).M(), (bTopHad.P4() + q2TopHad.P4()).M());
    
    
    return true;
}
