#include <SingleTopTChanPlugin_zerotag.hpp>

#include <Processor.hpp>
#include <ROOTLock.hpp>

#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <sys/stat.h>
#include <iostream>
#include <boost/filesystem.hpp>

using namespace std;


SingleTopTChanPlugin_zerotag::SingleTopTChanPlugin_zerotag(string const &outDirectory_, std::shared_ptr<BTagger const> &bTagger_, bool const isWeightSyst_):
    Plugin("SingleTop"),
    bTagger(bTagger_), outDirectory(outDirectory_), isWeightSyst(isWeightSyst_)
{
    // Make sure the directory path ends with a slash
    if (outDirectory.back() != '/')
        outDirectory += '/';
    
    // Create the output directory if it does not exist
    boost::filesystem::create_directories(outDirectory);
}


Plugin *SingleTopTChanPlugin_zerotag::Clone() const
{
    //return new SingleTopTChanPlugin(outDirectory, bTagger, syst);
//return new SingleTopTChanPlugin(outDirectory, bTagger);
//return new SingleTopTChanPlugin(outDirectory, bTagger, isWeightSyst);
    return new SingleTopTChanPlugin_zerotag(*this);
}


void SingleTopTChanPlugin_zerotag::BeginRun(Dataset const &dataset)
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
    
    tree->Branch("Pt_Lep", &Pt_Lep);
    tree->Branch("Eta_Lep", &Eta_Lep);
    tree->Branch("RelIso_Lep", &RelIso_Lep);
    tree->Branch("Charge_Lep", &Charge_Lep);
    tree->Branch("MET", &MET);
    tree->Branch("MtW", &MtW);
    tree->Branch("Phi_MET", &Phi_MET);
    tree->Branch("DPhi_LepNu", &DPhi_LepNu);
    
    tree->Branch("Pt_J1", &Pt_J1);
    tree->Branch("Eta_J1", &Eta_J1);
    tree->Branch("Pt_J2", &Pt_J2);
    tree->Branch("Eta_J2", &Eta_J2);
    tree->Branch("Pt_LJ", &Pt_LJ);
    tree->Branch("Eta_LJ", &Eta_LJ);
    tree->Branch("Pt_BJ1", &Pt_BJ1);
    tree->Branch("Pt_BJ2", &Pt_BJ2);
    
    tree->Branch("M_J1J2", &M_J1J2);
    tree->Branch("DR_J1J2", &DR_J1J2);
    tree->Branch("Pt_J1J2", &Pt_J1J2);
    tree->Branch("Ht_J1J2", &Ht_J1J2);
    
    tree->Branch("DR_LepJ1", &DR_LepJ1);
    tree->Branch("DR_LepJ2", &DR_LepJ2);
    tree->Branch("DPhi_LepJ1", &DPhi_LepJ1);
    
    tree->Branch("N_J", &N_J);
    tree->Branch("N_BJ", &N_BJ);
    tree->Branch("N_LJ", &N_LJ);
    tree->Branch("Ht", &Ht);
    tree->Branch("Ht_J", &Ht_J);
    tree->Branch("Ht_JNotBest", &Ht_JNotBest);
    tree->Branch("M_J", &M_J);
    tree->Branch("M_JNotBest", &M_JNotBest);
    tree->Branch("Pt_JNotBest", &Pt_JNotBest);
    tree->Branch("M_JW", &M_JW);
    tree->Branch("Pt_W", &Pt_W);
    
    tree->Branch("DPhi_LepW", &DPhi_LepW);
    tree->Branch("DPhi_LepBJ1", &DPhi_LepBJ1);
    tree->Branch("DPhi_WNu", &DPhi_WNu);
    tree->Branch("DPhi_WBJ1", &DPhi_WBJ1);
    tree->Branch("DR_LepBJ1", &DR_LepBJ1);
    tree->Branch("DR_WBJ1", &DR_WBJ1);
    
    tree->Branch("Mtop_BJ1", &Mtop_BJ1);
    tree->Branch("Mtop_BestJ", &Mtop_BestJ);
    tree->Branch("Pttop_BJ1", &Pttop_BJ1);
    tree->Branch("Cos_LepLJ_BJ1", &Cos_LepLJ_BJ1);
    tree->Branch("Cos_WLJ_BJ1", &Cos_WLJ_BJ1);
    tree->Branch("Cos_LepJ1", &Cos_LepJ1);
    tree->Branch("Cos_LepW_W", &Cos_LepW_W);
    
    tree->Branch("Sphericity", &Sphericity);
    tree->Branch("Aplanarity", &Aplanarity);
    tree->Branch("Planarity", &Planarity);
    
    tree->Branch("nPV", &nPV);
    
    if (dataset.IsMC())
    {
        tree->Branch("weight", &weight);
	if (isWeightSyst)
	{
	    tree->Branch("weight_PileUpUp", &weight_PileUpUp);
	    tree->Branch("weight_PileUpDown", &weight_PileUpDown);
	    tree->Branch("weight_TagRateUp", &weight_TagRateUp);
	    tree->Branch("weight_TagRateDown", &weight_TagRateDown);
	    tree->Branch("weight_MistagRateUp", &weight_MistagRateUp);
	    tree->Branch("weight_MistagRateDown", &weight_MistagRateDown);
	}
    }
}


void SingleTopTChanPlugin_zerotag::EndRun()
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


bool SingleTopTChanPlugin_zerotag::ProcessEvent()
{

    // Make sure the event contains reasonable physics objects
    if ((*reader)->GetLeptons().size() not_eq 1 or (*reader)->GetJets().size() < 2)
        return false;
 
    // Save event ID
    auto const &eventID = (*reader)->GetEventID();
    runNumber = eventID.Run();
    eventNumber = eventID.Event();
    lumiSection = eventID.LumiBlock();
    
    
    // Define some short-cuts
    auto const &lepton = (*reader)->GetLeptons().front();
    auto const &jets = (*reader)->GetJets();
    auto const &met = (*reader)->GetMET();
    
    // Calculate lepton-only variables
    Pt_Lep = lepton.Pt();
    Eta_Lep = lepton.Eta();
    RelIso_Lep = lepton.RelIso();
    Charge_Lep = lepton.Charge();
    MET = met.Pt();
    Phi_MET = met.Phi();
    DPhi_LepNu = fabs(lepton.Phi() - met.Phi());
    if (DPhi_LepNu > M_PI)
        DPhi_LepNu = 2 * M_PI - DPhi_LepNu;
    
    MtW = sqrt(pow(lepton.Pt() + met.Pt(), 2) - pow(lepton.P4().Px() + met.P4().Px(), 2) -
     pow(lepton.P4().Py() + met.P4().Py(), 2));
    
    // Find the light-flavour jet, the hardest b-jet and second b-jet
    unsigned index = 0;
    int b2jetIndex = -1;
    N_BJ = 0; //the first assumption
    Eta_LJ = 0.;
    
    for (unsigned i = 0; i < jets.size(); ++i)
        if (not (*bTagger)(jets.at(i)) and fabs(jets.at(i).Eta()) > fabs(Eta_LJ))
        {
            index = i;
            Eta_LJ = jets.at(i).Eta();
        }
    auto const &lJet = jets.at(index);
    
    //Find the first b-jet
    for (index = 0; index < jets.size(); ++index)
        if ((*bTagger)(jets.at(index)))
        {
            ++N_BJ;
            break;
        }
    
    if (index == jets.size())  // there are no tagged jets
    {
        index = 0;
        double maxCSV = -100.;
        
        for (unsigned i = 0; i < jets.size(); ++i)
            if (jets.at(i).CSV() > maxCSV)
            {
                index = i;
                maxCSV = jets.at(i).CSV();
            }
    }
    auto const &bJet = jets.at(index);
    
    // Find the best jet (the best for the top reconstruction)
    
    double massDelta = fabs((lepton.P4() + met.P4() + jets.at(0).P4()).M() - 172.5);
    unsigned bestJetIndex = 0;
    
    for (unsigned i = 1; i < jets.size(); ++i)
    {
         double massDeltaCand = fabs((lepton.P4() + met.P4() + jets.at(0).P4()).M() - 172.5);
         if (massDeltaCand < massDelta)
         {
             bestJetIndex = i;
             massDelta = massDeltaCand;
         }
    }
    
    auto const &bestJet = jets.at(bestJetIndex);
    
    // Calculate single-jet variables
    Pt_J1 = jets.at(0).Pt();
    Eta_J1 = jets.at(0).Eta();
    Pt_J2 = jets.at(1).Pt();
    Eta_J2 = jets.at(1).Eta();
    Pt_BJ1 = bJet.Pt();
    Pt_BJ2 =1.;
    Pt_LJ = lJet.Pt();
    
    // Calculate dijet variables
    M_J1J2 = (jets.at(0).P4() + jets.at(1).P4()).M();
    DR_J1J2 = jets.at(0).P4().DeltaR(jets.at(1).P4());
    Pt_J1J2 = (jets.at(0).P4() + jets.at(1).P4()).Pt();
    Ht_J1J2 = jets.at(0).P4().Pt() + jets.at(1).P4().Pt();
    
    
    // Calculate multi-jet variables
    N_J = jets.size();
    N_LJ = N_J - N_BJ;
    
    TLorentzVector p4Jets;
    Ht_J = 0.;
    Ht = 0.;
    
    for (auto const &j: jets)
    {
        p4Jets += j.P4();
        Ht_J += j.Pt();
        Ht += j.Pt();
    }
    
    for (auto const &j: (*reader)->GetAdditionalJets())
    {
        p4Jets += j.P4();
        Ht_J += j.Pt();
        Ht += j.Pt();
    }
    
    Ht_JNotBest = Ht_J - bestJet.Pt();
    M_J = p4Jets.M();
    M_JNotBest = (N_J > 2) ? (p4Jets - bestJet.P4()).M() : 1.;
    Pt_JNotBest = (p4Jets - bestJet.P4()).Pt();
    
    // Calculate lepton-jet variables
    Ht += lepton.Pt();
    Ht += met.Pt();
    DR_LepJ1 = lepton.P4().DeltaR(jets.at(0).P4());
    DR_LepJ2 = lepton.P4().DeltaR(jets.at(1).P4());
    DR_LepBJ1 = lepton.P4().DeltaR(bJet.P4());
    DPhi_LepJ1 = fabs(lepton.Phi() - jets.at(0).Phi());
    if (DPhi_LepJ1 > M_PI)
        DPhi_LepJ1 = 2 * M_PI - DPhi_LepJ1;
    
    // Reconstruct W-boson
    TLorentzVector const p4W((*reader)->GetNeutrino().P4() + lepton.P4());
    
    M_JW = (p4W + p4Jets).M();
    Pt_W = p4W.Pt();
    
    DPhi_LepW = fabs(lepton.Phi() - p4W.Phi());
    if (DPhi_LepW > M_PI)
        DPhi_LepW = 2 * M_PI - DPhi_LepW;
    
    DPhi_WNu = fabs(p4W.Phi() - met.Phi());
    if (DPhi_WNu > M_PI)
        DPhi_WNu = 2 * M_PI - DPhi_WNu;
    
    DPhi_WBJ1 = fabs(p4W.Phi() - bJet.Phi());
    if (DPhi_WBJ1 > M_PI)
        DPhi_WBJ1 = 2 * M_PI - DPhi_WBJ1;
    
    DPhi_LepBJ1 = fabs(lepton.Phi() - bJet.Phi());
    if (DPhi_LepBJ1 > M_PI)
        DPhi_LepBJ1 = 2 * M_PI - DPhi_LepBJ1;
    
    DR_WBJ1 = p4W.DeltaR(bJet.P4());
    
    // Reconstruct the top-quark
    TLorentzVector const p4Top(p4W + bJet.P4()); //with BJ1
    TLorentzVector const p4Top_Best(p4W + bestJet.P4()); //with BestJet
    
    Mtop_BJ1 = p4Top.M();
    Pttop_BJ1 = p4Top.Pt();
    Mtop_BestJ = p4Top_Best.M();
    
    
    // Calculate cos(theta)
    TVector3 b(p4Top.BoostVector());
    
    TLorentzVector boostedLepton(lepton.P4());
    boostedLepton.Boost(-b);
    TVector3 p3Lepton(boostedLepton.Vect());
    
    TLorentzVector boostedLJet(lJet.P4());
    boostedLJet.Boost(-b);
    TVector3 const p3LJet(boostedLJet.Vect());
    
    Cos_LepLJ_BJ1 = p3Lepton.Dot(p3LJet) / (p3Lepton.Mag() * p3LJet.Mag());
    
    // cos(theta) a la arXiv:1208.6006
    TLorentzVector boostedW = p4W;
    boostedW.Boost(-b);
    TVector3 p3W(boostedW.Vect());
    Cos_WLJ_BJ1 = p3W.Dot(p3LJet) / (p3W.Mag() * p3LJet.Mag());
    
    // cos(theta*) a la arXiv:1208.6006
    b = p4W.BoostVector();
    boostedLepton = lepton.P4();
    boostedLepton.Boost(-b);
    p3Lepton = boostedLepton.Vect();
    Cos_LepW_W = p3Lepton.Dot(-b) / (p3Lepton.Mag() * b.Mag());
    
    // cos(lepton, J1)_Lab
    p3Lepton = lepton.P4().Vect();
    TVector3 p3Jet = jets.at(0).P4().Vect();
    Cos_LepJ1 = p3Lepton.Dot(p3Jet) / (p3Lepton.Mag() * p3Jet.Mag());    
    
    // Calculate sphericity
    TMatrixDSym sphericityTensor(3);
    double norm = 0.;
    
    TVector3 p3(p4W.Vect());
    norm += p3.Mag2();
    
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            sphericityTensor(i, j) = p3[i] * p3[j];
    
    for (auto const &j: jets)
    {
        p3 = j.P4().Vect();
        norm += p3.Mag2();
        
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                sphericityTensor(i, j) += p3[i] * p3[j];
    }
    
    sphericityTensor *= 1. / norm;
    
    TMatrixDSymEigen eigenValCalc(sphericityTensor);
    TVectorD eigenVals(eigenValCalc.GetEigenValues());
    
    Sphericity = 1.5 * (eigenVals[1] + eigenVals[2]);
    Aplanarity = 1.5 * eigenVals[2];
    Planarity = eigenVals[1] - eigenVals[2];
    
    // Number of reconstructed primary vertices
    nPV = (*reader)->GetNPrimaryVertices();
    
    
    // Event weight
    weight = (*reader)->GetCentralWeight();

    if (isWeightSyst)
    {
	vector<WeightPair> const weights_PileUp = (*reader)->GetSystWeight(SystTypeWeight::PileUp);
	vector<WeightPair> const weights_TagRate = (*reader)->GetSystWeight(SystTypeWeight::TagRate);
	vector<WeightPair> const weights_MistagRate = (*reader)->GetSystWeight(SystTypeWeight::MistagRate);
	
	weight_PileUpUp = weights_PileUp.at(0).up;
	weight_PileUpDown = weights_PileUp.at(0).down;
	weight_TagRateUp = weights_TagRate.at(0).up;
	weight_TagRateDown = weights_TagRate.at(0).down;
	weight_MistagRateUp = weights_MistagRate.at(0).up;
	weight_MistagRateDown = weights_MistagRate.at(0).down;
    }
    
    tree->Fill();
    return true;
}
