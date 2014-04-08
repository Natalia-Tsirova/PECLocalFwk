#include <PECReader.hpp>

#include <PECReaderConfig.hpp>

#include <CalculatePzNu.hpp>
#include <ROOTLock.hpp>
#include <Logger.hpp>

#include <TVector3.h>
#include <TObjString.h>
#include <TVectorD.h>

#include <iostream>
#include <algorithm>
#include <stdexcept>


using namespace std;
using namespace logging;


// Definition of a static data member
unsigned const PECReader::maxSize;


PECReader::PECReader(Dataset const &dataset_):
    dataset(dataset_),
    isInitialized(false),
    triggerSelection(nullptr), eventSelection(nullptr),
    bTagReweighter(nullptr), puReweighter(nullptr),
    readHardParticles(false), readGenJets(false), readPartonShower(false),
    sourceFile(nullptr),
    eventIDTree(nullptr), triggerTree(nullptr), generalTree(nullptr)
{}


PECReader::PECReader(Dataset const &dataset, PECReaderConfig const &config):
    PECReader(dataset)
{
    Configure(config);
}


void PECReader::Configure(PECReaderConfig const &config)
{
    if (config.IsSetTriggerSelection())
        SetTriggerSelection(config.GetTriggerSelection());
    
    if (config.IsSetEventSelection())
        SetEventSelection(config.GetEventSelection());
    
    if (config.IsSetBTagReweighter())
        SetBTagReweighter(config.GetBTagReweighter());
    
    if (config.IsSetPileUpReweighter())
        SetPileUpReweighter(config.GetPileUpReweighter());
    
    SetReadHardInteraction(config.GetReadHardInteraction());
    SetReadGenJets(config.GetReadGenJets());
    SetReadPartonShower(config.GetReadPartonShower());
    SetSystematics(config.GetSystematics());
}


void PECReader::SetTriggerSelection(TriggerSelectionInterface *triggerSelection_)
{
    triggerSelection = triggerSelection_;
}


void PECReader::SetEventSelection(EventSelectionInterface const *eventSelection_)
{
    eventSelection = eventSelection_;
}


void PECReader::SetBTagReweighter(WeightBTagInterface *bTagReweighter_)
{
    bTagReweighter = bTagReweighter_;
}


void PECReader::SetPileUpReweighter(WeightPileUpInterface const *puReweighter_)
{
    puReweighter = puReweighter_;
}


void PECReader::SetReadHardInteraction(bool flag /*= true*/)
{
    readHardParticles = flag;
}


void PECReader::SetReadGenJets(bool flag /*= true*/)
{
    readGenJets = flag;
}


void PECReader::SetReadPartonShower(bool flag /*= true*/) noexcept
{
    readPartonShower = flag;
}


void PECReader::SetSystematics(SystTypeAlgo type, int direction /*= 0*/)
{
    syst.Set(type, direction);
}


void PECReader::SetSystematics(SystVariation const &syst_)
{
    syst = syst_;
}


bool PECReader::NextSourceFile()
{
    // Perform initialization
    if (not isInitialized)
        Initialize();
    
    
    // Close the currently opened file
    CloseSourceFile();
    
    
    // Check if there are files available and open the next one
    if (sourceFileIt == dataset.GetFiles().end())
        return false;
    
    OpenSourceFile();
    
    
    // Move to the next file
    ++sourceFileIt;
    
    return true;
}


bool PECReader::NextEvent()
{
    // Make sure there is a valid source file opened
    if (sourceFile == nullptr)
        throw logic_error("RECReader::NextEvent: No valid source file has been opened. Probably, "
         " PECReader::NextSourceFile has never been called.");
    
    
    while (true)
    {
        if (curEventTree == nEventsTree)  // there are no more events in the source files
            return false;
        
        
        // Read the event ID
        eventIDTree->GetEntry(curEventTree);
        eventID.Set(runNumber, lumiSection, eventNumber);
        
        
        // Update the event in the trigger-selection object and check if it passes the selection
        if (triggerSelection)
        {
            triggerSelection->ReadNextEvent(eventID);
            
            if (not triggerSelection->PassTrigger())
            {
                ++curEventTree;
                continue;
            }
        }
        
        
        // Read the rest of event
        generalTree->GetEntry(curEventTree);
        
        
        ++curEventTree;
        
        if (BuildAndSelectEvent())  // an appropriate event has been read
        {
            CalculateEventWeights();
            
            if (weightCentral not_eq 0.)
            {
                if (readHardParticles)
                    ParseHardInteraction();
                
                if (readGenJets and dataset.IsMC())
                    BuildGenJets();
                
                if (readPartonShower and dataset.IsMC())
                    ReadPartonShower();
                
                break;
            }
        }
    }
    
    return true;
}


EventID const &PECReader::GetEventID() const
{
    return eventID;
}


vector<Lepton> const &PECReader::GetLeptons() const
{
    return tightLeptons;
}


vector<Jet> const &PECReader::GetJets() const
{
    return goodJets;
}


vector<Jet> const &PECReader::GetAdditionalJets() const
{
    return additionalJets;
}


Candidate const &PECReader::GetMET() const
{
    return correctedMET;
}


Candidate const &PECReader::GetNeutrino() const
{
    return neutrino;
}


unsigned PECReader::GetNPrimaryVertices() const
{
    return pvSize;
}


double PECReader::GetRho() const
{
    return puRho;
}


double PECReader::GetCentralWeight() const
{
    return weightCentral;
}

double PECReader::GetXSecWeight() const
{
    return weightCrossSection;
}

double PECReader::GetBTagWeight() const
{
    return weightBTag;
}

double PECReader::GetTriggerWeight() const
{
    return weightTrigger;
}

double PECReader::GetPUWeight() const
{
    return weightPU;
}


vector<WeightPair> const &PECReader::GetSystWeight(SystTypeWeight type) const
{
    if (syst.type != SystTypeAlgo::WeightOnly)
        throw logic_error("PECReader::GetSystWeight: Trying to access systematical shifts in event "
         "weight without requesting them.");
    
    switch (type)
    {
        case SystTypeWeight::PileUp:
            return systWeightPileUp;
        
        case SystTypeWeight::TagRate:
            return systWeightTagRate;
        
        case SystTypeWeight::MistagRate:
            return systWeightMistagRate;
        
        default:
            throw logic_error("PECReader::GetSystWeight: Requested variation is not supported.");
    }
}


vector<GenParticle> const &PECReader::GetHardGenParticles() const
{
    if (not readHardParticles)
        throw runtime_error("PECReader::GetHardGenParticles: In order to access the list of "
         "generator particles associated to the hard interaction, this functionality must first "
         "be requested via PECReader::SetReadHardInteraction.");
    
    return hardParticles;
}


vector<GenJet> const &PECReader::GetGenJets() const
{
    if (not readGenJets)
        throw runtime_error("PECReader::GetGenJets: In order to access generator-level jets, this "
         "functionality must first be requested via PECReader::SetReadGenJets.");
    
    if (not dataset.IsMC())
        throw runtime_error("PECReader::GetGenJets: Trying to get generatol-level jets in a "
         "real collision event.");
    
    return genJets;
}


vector<ShowerParton> const &PECReader::GetShowerPartons() const
{
    if (not readPartonShower)
        throw runtime_error("PECReader::GetShowerPartons: In order to access partons from parton "
         "shower, this functionality must first be requested via PECReader::SetReadPartonShower.");
    
    if (not dataset.IsMC())
        throw runtime_error("PECReader::GetShowerPartons: Trying to read partons from parton "
         "shower  for a real collision event.");
    
    return psPartons;
}


void PECReader::Initialize()
{
    // Verify that all the needed configuration modules have been specified
    if (not triggerSelection)
        logger << "Warning in PECReader::Initialize: No trigger selection has been specified." <<
         eom;
    
    if (not eventSelection)
        logger << "Warning in PECReader::Initialize: No event selection has been specified." << eom;
    
    if (dataset.IsMC())
    {
        if (not bTagReweighter)
            logger << "Warning in PECReader::Initialize: No object to propagate b-tagging scale " <<
             "factors has been specified. Simulation will not be reweighted for this effect." <<
             eom;
        
        if (not puReweighter)
            logger << "Warning in PECReader::Initialize: No object to reweight simulation for " <<
             "pile-up has been specified. Simulation will not be reweighted for this effect." <<
             eom;
    }
    
    
    // Perform remaining initialization
    sourceFileIt = dataset.GetFiles().begin();
    
    
    // Indicate that initialization has been performed
    isInitialized = true;
}


void PECReader::OpenSourceFile()
{
    // Set the event weight due to the cross-section
    if (dataset.IsMC())
        weightCrossSection = sourceFileIt->xSec / sourceFileIt->nEvents;
    else
        weightCrossSection = 1.;
    
    
    // Start of a critical ROOT block
    ROOTLock::Lock();
    
    // Open the source file
    sourceFile = TFile::Open(sourceFileIt->name.c_str());
    
    if (not sourceFile)
        throw runtime_error(string("PECReader::OpenSourceFile: File \"") + sourceFileIt->name +
         "\" does not exist or is not a valid ROOT file.");
    
    
    // Get the trees
    eventIDTree = dynamic_cast<TTree *>(sourceFile->Get("eventContent/EventID"));
    
    if (triggerSelection)
    {
        triggerTree = dynamic_cast<TTree *>(sourceFile->Get("trigger/TriggerInfo"));
        triggerSelection->UpdateTree(triggerTree, not dataset.IsMC());
    }
    
    generalTree = dynamic_cast<TTree *>(sourceFile->Get("eventContent/BasicInfo"));
    //generalTree->AddFriend("eventContent/IntegralProperties");
    generalTree->AddFriend("eventContent/BasicInfo");
    generalTree->AddFriend("eventContent/PUInfo");
    
    
    // Add the MC-truth information and MC weights
    if (dataset.IsMC())
    {
        generalTree->AddFriend("eventContent/GeneratorInfo");
        
        if (readGenJets)
            generalTree->AddFriend("genJets/GenJets");
        
        if (readPartonShower)
            generalTree->AddFriend("heavyFlavours/PartonShowerInfo");
        
        /*
        // Determine name of the corresponding file with weights
        string weightsFileName;
        
        if (config.IsSetWeightFilesLocation())
            weightsFileName = config.GetWeightFilesLocation() + sourceFileIt->GetBaseName() +
             "_weights.root";
        else
            weightsFileName = sourceFileIt->name.substr(0,
             sourceFileIt->name.find_last_of('.')) + "_weights.root";
        
        // Add trees from the file with weights
        generalTree->AddFriend(config.GetPileUpTreeName().c_str(), weightsFileName.c_str());
        */
    }
    
    // The file is opened, all the trees are got. Can end the ROOT critical block
    ROOTLock::Unlock();
    
    
    // Initialize the counters
    nEventsTree = generalTree->GetEntries();  // all the trees have the same number of events
    curEventTree = 0;
    
    
    // Assign the branches to read
    eventIDTree->SetBranchAddress("run", &runNumber);
    eventIDTree->SetBranchAddress("lumi", &lumiSection);
    eventIDTree->SetBranchAddress("event", &eventNumber);
    
    generalTree->SetBranchAddress("eleSize", &eleSize);
    generalTree->SetBranchAddress("elePt", elePt);
    generalTree->SetBranchAddress("eleEta", eleEta);
    generalTree->SetBranchAddress("elePhi", elePhi);
    generalTree->SetBranchAddress("eleRelIso", eleRelIso);
    generalTree->SetBranchAddress("eleDB", eleDB);
    generalTree->SetBranchAddress("eleTriggerPreselection", eleTriggerPreselection);
    generalTree->SetBranchAddress("eleMVAID", eleMVAID);
    generalTree->SetBranchAddress("elePassConversion", elePassConversion);
    generalTree->SetBranchAddress("eleSelectionA", eleQuality);
    generalTree->SetBranchAddress("eleCharge", eleCharge);
    
    generalTree->SetBranchAddress("muSize", &muSize);
    generalTree->SetBranchAddress("muPt", muPt);
    generalTree->SetBranchAddress("muEta", muEta);
    generalTree->SetBranchAddress("muPhi", muPhi);
    generalTree->SetBranchAddress("muRelIso", muRelIso);
    generalTree->SetBranchAddress("muDB", muDB);
    generalTree->SetBranchAddress("muQualityTight", muQualityTight);
    generalTree->SetBranchAddress("muCharge", muCharge);
    
    generalTree->SetBranchAddress("jetSize", &jetSize);
    generalTree->SetBranchAddress("jetPt", jetPt);
    generalTree->SetBranchAddress("jetEta", jetEta);
    generalTree->SetBranchAddress("jetPhi", jetPhi);
    generalTree->SetBranchAddress("jetMass", jetMass);
    
    if (dataset.IsMC() and syst.type == SystTypeAlgo::JER)
    {
        if (syst.direction > 0)
            generalTree->SetBranchAddress("jerFactorUp", jerFactor);
        else
            generalTree->SetBranchAddress("jerFactorDown", jerFactor);
    }
    
    /*
    generalTree->SetBranchAddress("softJetPt", &softJetPt);
    generalTree->SetBranchAddress("softJetEta", &softJetEta);
    generalTree->SetBranchAddress("softJetPhi", &softJetPhi);
    generalTree->SetBranchAddress("softJetMass", &softJetMass);
    generalTree->SetBranchAddress("softJetHt", &softJetHt);
    
    if (dataset.IsMC() and syst.type == SystTypeAlgo::JER)
    {
        if (systDir > 0)
        {
            generalTree->SetBranchAddress("softJetPtJERUp", &softJetPt);
            generalTree->SetBranchAddress("softJetEtaJERUp", &softJetEta);
            generalTree->SetBranchAddress("softJetPhiJERUp", &softJetPhi);
            generalTree->SetBranchAddress("softJetMassJERUp", &softJetMass);
            generalTree->SetBranchAddress("softJetHtJERUp", &softJetHt);
        }
        else
        {
            generalTree->SetBranchAddress("softJetPtJERDown", &softJetPt);
            generalTree->SetBranchAddress("softJetEtaJERDown", &softJetEta);
            generalTree->SetBranchAddress("softJetPhiJERDown", &softJetPhi);
            generalTree->SetBranchAddress("softJetMassJERDown", &softJetMass);
            generalTree->SetBranchAddress("softJetHtJERDown", &softJetHt);
        }
    */
    
    generalTree->SetBranchAddress("jetCSV", jetCSV);
    generalTree->SetBranchAddress("jetTCHP", jetTCHP);
    
    generalTree->SetBranchAddress("jetCharge", jetCharge);
    generalTree->SetBranchAddress("jetPullAngle", jetPullAngle);
    
    generalTree->SetBranchAddress("metSize", &metSize);
    generalTree->SetBranchAddress("metPt", metPt);
    generalTree->SetBranchAddress("metPhi", metPhi);
    
    generalTree->SetBranchAddress("pvSize", &pvSize);
    generalTree->SetBranchAddress("rho", &puRho);
    
    
    if (dataset.IsMC())
    {
        generalTree->SetBranchAddress("jetFlavour", jetFlavour);
        generalTree->SetBranchAddress("processID", &processID);
        
        // Some systematics is encoded in weights only. These are added to the central samples only
        /*
        if (syst == Systematics::None)
        {
            generalTree->SetBranchAddress("PDF.nVars", &nWeight_PDF);
            generalTree->SetBranchAddress("PDF.up", weight_PDFUp);
        }
        */
        
        if (syst.type == SystTypeAlgo::JEC)
        {
            generalTree->SetBranchAddress("jecUncertainty", jecUncertainty);
            
            /*
            generalTree->SetBranchAddress("softJetPtJECUnc", &softJetPtJECUnc);
            generalTree->SetBranchAddress("softJetEtaJECUnc", &softJetEtaJECUnc);
            generalTree->SetBranchAddress("softJetPhiJECUnc", &softJetPhiJECUnc);
            generalTree->SetBranchAddress("softJetMassJECUnc", &softJetMassJECUnc);
            generalTree->SetBranchAddress("softJetHtJECUnc", &softJetHtJECUnc);
            */
        }
        
        
        // Generator jets
        if (readGenJets)
        {
            generalTree->SetBranchAddress("genJets/GenJets.jetSize", &genJetSize);
            generalTree->SetBranchAddress("genJets/GenJets.jetPt", &genJetPt);
            generalTree->SetBranchAddress("genJets/GenJets.jetEta", &genJetEta);
            generalTree->SetBranchAddress("genJets/GenJets.jetPhi", &genJetPhi);
            generalTree->SetBranchAddress("genJets/GenJets.jetMass", &genJetMass);
            //generalTree->SetBranchAddress("genJets/GenJets.bMult", &genJetBMult);
            //generalTree->SetBranchAddress("genJets/GenJets.cMult", &genJetCMult);
        }
        
        
        // Partons from parton shower
        if (readPartonShower)
        {
            generalTree->SetBranchAddress("psSize", &psSize);
            generalTree->SetBranchAddress("psPdgId", psPdgId);
            generalTree->SetBranchAddress("psOrigin", psOrigin);
            generalTree->SetBranchAddress("psPt", psPt);
            generalTree->SetBranchAddress("psEta", psEta);
            generalTree->SetBranchAddress("psPhi", psPhi);
        }
        
        
        // Pile-up information
        generalTree->SetBranchAddress("puTrueNumInteractions", &puTrueNumInteractions);
    }
    
    if (dataset.IsMC() and readHardParticles)
    {
        generalTree->SetBranchAddress("hardPartSize", &hardPartSize);
        generalTree->SetBranchAddress("hardPartPdgId", hardPartPdgId);
        generalTree->SetBranchAddress("hardPartFirstMother", hardPartFirstMother);
        generalTree->SetBranchAddress("hardPartLastMother", hardPartLastMother);
        generalTree->SetBranchAddress("hardPartPt", hardPartPt);
        generalTree->SetBranchAddress("hardPartEta", hardPartEta);
        generalTree->SetBranchAddress("hardPartPhi", hardPartPhi);
        generalTree->SetBranchAddress("hardPartMass", hardPartMass);
    }
}


void PECReader::CloseSourceFile()
{
    // Delete the source file and trees (it is a critical section)
    ROOTLock::Lock();
    
    delete eventIDTree;
    delete triggerTree;
    delete generalTree;
    delete sourceFile;
    
    ROOTLock::Unlock();
    
    
    // Set the above pointers to nulls to indicate that the file has been closed
    sourceFile = nullptr;
    eventIDTree = triggerTree = generalTree = nullptr;
}


bool PECReader::BuildAndSelectEvent()
{
    // Filter inclusive Wjets dataset if needed
    if (dataset.GetProcess() == Dataset::Process::Wjets and dataset.TestFlag("WjetsKeep0p1p") and
     processID % 5 > 1)
        return false;
    
    
    // Reset the containers used in the compact event description
    tightLeptons.clear();
    looseLeptons.clear();
    goodJets.clear();
    additionalJets.clear();
    
    
    // Loop over the electrons
    for (int i = 0; i < eleSize; ++i)
    {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(elePt[i], eleEta[i], elePhi[i], 0.511e-3);
        
        
        if (p4.Pt() < 20. or fabs(p4.Eta()) > 2.5 or eleRelIso[i] > 0.15)
            continue;
        
        // A loose electron is found
        Lepton lepton(Lepton::Flavour::Electron, p4);
        lepton.SetRelIso(eleRelIso[i]);
        lepton.SetDB(eleDB[i]);
        lepton.SetCharge((eleCharge[i]) ? -1 : 1);
        
        looseLeptons.push_back(lepton);
        
        
        if (p4.Pt() < 20. or not eleQuality[i] or eleRelIso[i] > 0.1 or not elePassConversion[i]
         or not eleTriggerPreselection[i] or eleMVAID[i] < 0.9)
        //^ The pt cut is the same as for the loose electrons
            continue;
        
        // A tight electron is found
        tightLeptons.push_back(lepton);
    }
    
    
    // Loop over the muons
    for (int i = 0; i < muSize; ++i)
    {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(muPt[i], muEta[i], muPhi[i], 0.105);
        
        
        if (p4.Pt() < 10. or fabs(p4.Eta()) > 2.5 or muRelIso[i] > 0.2)
                continue;
        
        // A loose muon is found
        Lepton lepton(Lepton::Flavour::Muon, p4);
        lepton.SetRelIso(muRelIso[i]);
        lepton.SetDB(muDB[i]);
        lepton.SetCharge((muCharge[i]) ? -1 : 1);
        
        looseLeptons.push_back(lepton);
        
        
        if (p4.Pt() < 10. or fabs(p4.Eta()) > 2.1 or not muQualityTight[i] or fabs(muDB[i]) > 0.2 or
         muRelIso[i] > 0.12)
        //^ The pt cut is the same as for loose muons
            continue;
        
        // A tight muon is found
        tightLeptons.push_back(lepton);
    }
    
    
    // Leptonic step of the event selection
    if (eventSelection and not eventSelection->PassLeptonStep(tightLeptons, looseLeptons))
        return false;
    
    
    // Loop over the jets
    for (int i = 0; i < jetSize; ++i)
    {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(jetPt[i], jetEta[i], jetPhi[i], jetMass[i]);
        
        
        // Vary jet four-momentum within JEC uncertainty
        if (syst.type == SystTypeAlgo::JEC)
            p4 *= 1. + syst.direction * jecUncertainty[i];
        
        // Rescale jet four-momentum to account for JER systematical variation
        if (syst.type == SystTypeAlgo::JER)
            p4 *= jerFactor[i];
        
        
        // Reject too soft or too forward jets
        if (p4.Pt() < 20. or fabs(p4.Eta()) > 4.7)
            continue;
        
        
        Jet jet(p4);
        
        jet.SetCSV(jetCSV[i]);
        jet.SetTCHP(jetTCHP[i]);
        
        jet.SetCharge(jetCharge[i]);
        jet.SetPullAngle(jetPullAngle[i]);
        
        if (dataset.IsMC())
            jet.SetParentID(jetFlavour[i]);
        
        if (not eventSelection or eventSelection->IsAnalysisJet(jet))
            goodJets.push_back(jet);
        else
            additionalJets.push_back(jet);
    }
    
        
    // Make sure the jets are ordered in pt (in decreasing order). The ordering might have been
    //broken after the JER smearing was performed
    sort(goodJets.rbegin(), goodJets.rend());
    sort(additionalJets.rbegin(), additionalJets.rend());
    
    // Event selection on the number of jets and tags
    if (eventSelection and not eventSelection->PassJetStep(goodJets))
        return false;
    
    
    // Several versions of MET are stored in a PEC file
    unsigned metIndex = 1;  // index of a corrected MET that is not varied for some systematics
    
    switch (syst.type)
    {
        case SystTypeAlgo::JEC:
            // Check [1] and around that line
            //[1] https://svnweb.cern.ch/trac/singletop/browser/tags/2012Alpha_v2/CMSSW/SingleTop/python/ObjectsDefinitions_cff.py#L342
            metIndex = (syst.direction > 0) ? 2 : 3;
            break;
        
        case SystTypeAlgo::JER:
            metIndex = (syst.direction > 0) ? 4 : 5;
            break;
        
        case SystTypeAlgo::METUnclustered:
            metIndex = (syst.direction > 0) ? 6 : 7;
            break;
        
        default:
            break;
    }
    
    
    // For unknown reason extremely rarely MET can be NaN. Check for it
    if (isnan(metPt[metIndex]) or isnan(metPhi[metIndex]))
    {
        logger << "Warning: MET is NaN in event #" << curEventTree << " in file \"" <<
         sourceFile->GetName() << "\" (ID " << runNumber << ":" << lumiSection << ":" <<
         eventNumber << "). The event is skipped." << eom;
        return false;
    }
    
    
    // Save MET to the dedicated variable
    correctedMET.SetPtEtaPhiM(metPt[metIndex], 0., metPhi[metIndex], 0.);
    
    
    // Reconstruct the neutrino with the leading tight lepton
    double const nuPz =
     Nu4Momentum(tightLeptons.front().P4(), metPt[metIndex], metPhi[metIndex]).Pz();
    double const nuEnergy = sqrt(metPt[metIndex] * metPt[metIndex] + nuPz * nuPz);
    neutrino.SetPtEtaPhiM(metPt[metIndex], 0.5 * log((nuEnergy + nuPz) / (nuEnergy - nuPz)),
     metPhi[metIndex], 0.);
    
        
    return true;
}


void PECReader::CalculateEventWeights()
{
    // Calculate weight due to trigger selection. This is the only event weight that can make
    //sense for real data (e.g. if there is an additional selection specified in a TriggerRange
    //object). However, with real data this weight is either 0. or 1.
    weightTrigger = (triggerSelection) ? triggerSelection->GetWeight(*this) : 1.;
    
    if (not dataset.IsMC())
    {
        // Don't forget to update the central weight
        weightCentral = weightTrigger;
        
        return;
    }
    
    
    // Reweighting for the pile-up
    WeightPileUpInterface::Weights weightPileUp;
    
    if (puReweighter)
        weightPileUp = puReweighter->GetWeights(puTrueNumInteractions);
    else
        weightPileUp.Set(1., 1., 1.);
    
    
    // Reweighting for b-tagging
    weightBTag =
     (bTagReweighter) ?
     bTagReweighter->CalcWeight(goodJets, WeightBTagInterface::Variation::Nominal) : 1.;
    
    weightPU = weightPileUp.central;
    
    // Calculate the central weight
    weightCentral = weightCrossSection * weightTrigger * weightPileUp.central * weightBTag;
    
    
    // Clear vectors with varied weights
    systWeightPileUp.clear();
    systWeightTagRate.clear();
    systWeightMistagRate.clear();
    
    
    // Fill the vectors with varied weights
    if (syst.type == SystTypeAlgo::WeightOnly and puReweighter)
    {
        systWeightPileUp.emplace_back();
        systWeightPileUp.back().up = weightCentral / weightPileUp.central * weightPileUp.up;
        systWeightPileUp.back().down = weightCentral / weightPileUp.central * weightPileUp.down;
    }
    
    if (syst.type == SystTypeAlgo::WeightOnly and bTagReweighter)
    {
        double const weightButBTagging = weightCentral / weightBTag;
        
        systWeightTagRate.emplace_back();
        systWeightTagRate.back().up = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTagInterface::Variation::TagRateUp);
        systWeightTagRate.back().down = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTagInterface::Variation::TagRateDown);
        
        systWeightMistagRate.emplace_back();
        systWeightMistagRate.back().up = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTagInterface::Variation::MistagRateUp);
        systWeightMistagRate.back().down = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTagInterface::Variation::MistagRateDown);
    }
}


void PECReader::ParseHardInteraction()
{
    // Reset the vector
    hardParticles.clear();
    
    // Make sure the vector will not be reallocated (otherwise the pointers will be broken)
    hardParticles.reserve(hardPartSize);
    
    
    // Loop over the arrays with properties of the particles from the hard interaction
    for (unsigned i = 0; i < unsigned(hardPartSize); ++i)
    {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(hardPartPt[i], hardPartEta[i], hardPartPhi[i], hardPartMass[i]);
        
        hardParticles.emplace_back(p4, hardPartPdgId[i]);
        
        
        // Set pointers to mothers and daughters
        int iMother = hardPartFirstMother[i];
        
        if (iMother >= 0 and iMother < hardPartSize)
        {
            hardParticles.at(i).AddMother(&hardParticles.at(iMother));
            hardParticles.at(iMother).AddDaughter(&hardParticles.at(i));
        }
        
        if (hardPartFirstMother[i] != hardPartLastMother[i])
        //^ Second mother is added if only it is different from the first one
        {
            int iMother = hardPartLastMother[i];
            
            if (iMother >= 0 and iMother < hardPartSize)
            {
                hardParticles.at(i).AddMother(&hardParticles.at(iMother));
                hardParticles.at(iMother).AddDaughter(&hardParticles.at(i));
            }
        }
    }
}


void PECReader::BuildGenJets()
{
    // Reset the vector
    genJets.clear();
    
    
    // Loop over the generator jets and fill the vector
    for (unsigned i = 0; i < unsigned(genJetSize); ++i)
    {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(genJetPt[i], genJetEta[i], genJetPhi[i], genJetMass[i]);
        
        genJets.emplace_back(p4);
        //genJets.back().SetMultiplicities(genJetBMult[i], genJetCMult[i]);
    }
}


void PECReader::ReadPartonShower()
{
    // Reset the vector
    psPartons.clear();
    
    
    // Loop over the buffers and fill the vector
    for (unsigned i = 0; i < psSize; ++i)
    {
        ShowerParton::Origin origin = ShowerParton::Origin::Unknown;
        
        if (psOrigin[i] == 1)
            origin = ShowerParton::Origin::ISR;
        else if (psOrigin[i] == 2)
            origin = ShowerParton::Origin::FSR;
        else if (psOrigin[i] == 3)
            origin = ShowerParton::Origin::Proton;
        
        
        psPartons.emplace_back(psPt[i], psEta[i], psPhi[i], psPdgId[i], origin);
    }
}
