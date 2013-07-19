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


PECReader::PECReader(Dataset const &dataset_):
    dataset(dataset_),
    isInitialized(false),
    triggerSelection(nullptr), eventSelection(nullptr), puReweighter(nullptr),
    readHardParticles(false),
    bTagReweighter(nullptr),
    sourceFile(nullptr),
    eventIDTree(nullptr), triggerTree(nullptr), generalTree(nullptr),
    triggerNames(nullptr)
{}


PECReader::PECReader(Dataset const &dataset, PECReaderConfig const &config):
    PECReader(dataset)
{
    Configure(config);
}


PECReader::~PECReader()
{
    delete triggerNames;
}


void PECReader::Configure(PECReaderConfig const &config)
{
    if (config.IsSetTriggerSelection())
        SetTriggerSelection(config.GetTriggerSelection());
    
    if (config.IsSetEventSelection())
        SetEventSelection(config.GetEventSelection());
    
    if (config.IsSetBTagger() and config.IsSetBTagDatabase())
        SetBTaggingConfig(config.GetBTagger(), config.GetBTagDatabase());
    
    if (config.IsSetPileUpReweighter())
        SetPileUpReweighter(config.GetPileUpReweighter());
    
    SetReadHardInteraction(config.GetReadHardInteraction());
    SetSystematics(config.GetSystematics());
}


void PECReader::SetTriggerSelection(TriggerSelectionInterface const *triggerSelection_)
{
    triggerSelection = triggerSelection_;
}


void PECReader::SetEventSelection(EventSelectionInterface const *eventSelection_)
{
    eventSelection = eventSelection_;
}


void PECReader::SetBTaggingConfig(BTagger const *bTagger, BTagDatabase const *bTagDatabase)
{
    bTagReweighter.reset(new WeightBTag(*bTagger, *bTagDatabase));
}


void PECReader::SetPileUpReweighter(WeightPileUpInterface const *puReweighter_)
{
    puReweighter = puReweighter_;
}


void PECReader::SetReadHardInteraction(bool flag /*= true*/)
{
    readHardParticles = flag;
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
        
        
        // Read the event ID and trigger information
        eventIDTree->GetEntry(curEventTree);
        triggerTree->GetEntry(curEventTree);
        
        eventID.Set(runNumber, lumiSection, eventNumber);
        
        // In case of MC disable the branch with the trigger names: it is needed only once
        if (curEventTree == 0 and dataset.IsMC())
            triggerTree->SetBranchStatus("names", false);

        
        
        // Check if the event meets the trigger selection
        if (not triggerSelection or
         not triggerSelection->PassTrigger(eventID, triggerNames, hasFired))
        {
            ++curEventTree;
            continue;
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
    return PVSize;
}


double PECReader::GetCentralWeight() const
{
    return weightCentral;
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
             "factots has been specified. Simulation will not be reweighted for this effect." <<
             eom;
        
        if (not puReweighter)
            logger << "Warning in PECReader::Initialize: No object to reweight simulation for " <<
             "pile-up has been specified. Simulation will not be reweighted for this effect." <<
             eom;
    }
    
    
    // Create dynamically-allocated ROOT objects. Because of design of ROOT memory system, this
    //block of code must be atomic
    ROOTLock::Lock();
    
    triggerNames = new TClonesArray("TObjString");
    
    ROOTLock::Unlock();
    
    
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
    
    // Get the trees
    eventIDTree = dynamic_cast<TTree *>(sourceFile->Get("eventContent/EventID"));
    triggerTree = dynamic_cast<TTree *>(sourceFile->Get("trigger/TriggerInfo"));
    generalTree = dynamic_cast<TTree *>(sourceFile->Get("eventContent/BasicInfo"));
    //generalTree->AddFriend("eventContent/IntegralProperties");
    generalTree->AddFriend("eventContent/BasicInfo");
    generalTree->AddFriend("eventContent/PUInfo");
    
    // Add the MC-truth information and MC weights
    if (dataset.IsMC())
    {
        generalTree->AddFriend("eventContent/GeneratorInfo");
        
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
    
    triggerTree->SetBranchAddress("size", &triggerSize);
    triggerTree->SetBranchAddress("names", &triggerNames);
    triggerTree->SetBranchAddress("hasFired", hasFired);
    
    generalTree->SetBranchAddress("eleSize", &eleSize);
    generalTree->SetBranchAddress("elePt", elePt);
    generalTree->SetBranchAddress("eleEta", eleEta);
    generalTree->SetBranchAddress("elePhi", elePhi);
    generalTree->SetBranchAddress("eleRelIso", eleRelIso);
    generalTree->SetBranchAddress("eleDB", eleDB);
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
    generalTree->SetBranchAddress("muSelectionA", muQualityLoose);
    generalTree->SetBranchAddress("muSelectionB", muQualityTight);
    generalTree->SetBranchAddress("muCharge", muCharge);
    
    generalTree->SetBranchAddress("jetSize", &jetSize);
    generalTree->SetBranchAddress("jetEta", jetEta);
    generalTree->SetBranchAddress("jetPhi", jetPhi);
    
    if (dataset.IsMC() and syst.type == SystTypeAlgo::JER)
    {
        if (syst.direction > 0)
        {
            generalTree->SetBranchAddress("jetPtJERUp", jetPt);
            generalTree->SetBranchAddress("jetMassJERUp", jetMass);
        }
        else
        {
            generalTree->SetBranchAddress("jetPtJERDown", jetPt);
            generalTree->SetBranchAddress("jetMassJERDown", jetMass);
        }
    }
    else
    {
        generalTree->SetBranchAddress("jetPt", jetPt);
        generalTree->SetBranchAddress("jetMass", jetMass);
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
    //generalTree->SetBranchAddress("jetTCHE", jetTCHE);
    //generalTree->SetBranchAddress("jetJP", jetJP);
    
    generalTree->SetBranchAddress("metSize", &metSize);
    generalTree->SetBranchAddress("metPt", metPt);
    generalTree->SetBranchAddress("metPhi", metPhi);
    
    generalTree->SetBranchAddress("PVSize", &PVSize);
    
    
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
            generalTree->SetBranchAddress("PDF.down", weight_PDFDown);
            
            if (WjetsMG)
            {
                generalTree->SetBranchAddress("WjetQ2Weights/Q2Weights.up", &weight_WjetsQ2Up);
                generalTree->SetBranchAddress("WjetQ2Weights/Q2Weights.down", &weight_WjetsQ2Down);
            }
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
        
        /*
        if (WjetsMG)
            generalTree->SetBranchAddress("simpleClass", &WHFClass);
        */
        
        
        // Pile-up information
        generalTree->SetBranchAddress("PUTrueNumInteractions", &puTrueNumInteractions);
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
    
    
    // Notify the trigger object that a new file has been opened
    if (triggerSelection)
        triggerSelection->NewFile(not dataset.IsMC());
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
        
        
        if (p4.Pt() < 20. or fabs(p4.Eta()) > 2.5 or eleRelIso[i] > 0.15 or eleMVAID[i] < 0.)
            continue;
        
        // A loose electron is found
        Lepton lepton(Lepton::Flavour::Electron, p4);
        lepton.SetRelIso(eleRelIso[i]);
        lepton.SetDB(eleDB[i]);
        lepton.SetCharge((eleCharge[i]) ? -1 : 1);
        
        looseLeptons.push_back(lepton);
        
        
        if (p4.Pt() < 20. or !eleQuality[i] or eleRelIso[i] > 0.1 or not elePassConversion[i]
         or eleMVAID[i] < 0.5)
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
        
        
        if (p4.Pt() < 10. or fabs(p4.Eta()) > 2.5 or !muQualityLoose[i] or muRelIso[i] > 0.2)
                continue;
        
        // A loose muon is found
        Lepton lepton(Lepton::Flavour::Muon, p4);
        lepton.SetRelIso(muRelIso[i]);
        lepton.SetDB(muDB[i]);
        lepton.SetCharge((muCharge[i]) ? -1 : 1);
        
        looseLeptons.push_back(lepton);
        
        
        if (p4.Pt() < 10. or fabs(p4.Eta()) > 2.1 or !muQualityTight[i] or fabs(muDB[i]) > 0.2 or
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
        
        if (syst.type == SystTypeAlgo::JEC)
            p4 *= 1. + syst.direction * jecUncertainty[i];
        
        
        if (p4.Pt() < 20. or fabs(p4.Eta()) > 4.7)
            continue;
        
        
        Jet jet(p4);
        jet.SetCSV(jetCSV[i]);
        jet.SetTCHP(jetTCHP[i]);
        
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
    
    // In 2012Alpha_v2 campaign in case of MC wrong MET is stored under index 1. Below, a hack is
    //inserted to correct for it. In muon channel an event contains to loose electrons. Therefore,
    //variations of MET due to electron energy scale reproduce the correct central value for MET.
    //Similarly for electron channel
    if (dataset.IsMC())
    {
        if (tightLeptons.front().GetFlavour() == Lepton::Flavour::Muon)
        //^ I assume there is exactly one tight lepton
            metIndex = 8;  // electron energy up
        else
            metIndex = 10;  // muon energy up
    }
    
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
        logger << "WARNING: MET is NaN in event #" << curEventTree << " in file \"" <<
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
    // Calculate different components of an event weight
    double const weightTrigger = (triggerSelection) ? triggerSelection->GetWeight(*this) : 1.;
    
    
    WeightPileUpInterface::Weights weightPileUp;
    
    if (puReweighter)
        weightPileUp = puReweighter->GetWeights(puTrueNumInteractions);
    else
        weightPileUp.Set(1., 1., 1.);
    
    
    double const weightBTagging =
     (bTagReweighter) ? bTagReweighter->CalcWeight(goodJets, WeightBTag::Variation::Central) : 1.;
    
    
    // Calculate the central weight
    weightCentral = weightCrossSection * weightTrigger * weightPileUp.central * weightBTagging;
    
    
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
        double const weightButBTagging = weightCentral / weightBTagging;
        
        systWeightTagRate.emplace_back();
        systWeightTagRate.back().up = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTag::Variation::TagRateUp);
        systWeightTagRate.back().down = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTag::Variation::TagRateDown);
        
        systWeightMistagRate.emplace_back();
        systWeightMistagRate.back().up = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTag::Variation::MistagRateUp);
        systWeightMistagRate.back().down = weightButBTagging *
         bTagReweighter->CalcWeight(goodJets, WeightBTag::Variation::MistagRateDown);
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
        
        
        for (int iMother: {hardPartFirstMother[i], hardPartLastMother[i]})
        //^ Isn't C++11 just awesome? ;)
            if (iMother >= 0 and iMother < hardPartSize)
            {
                hardParticles.at(i).AddMother(&hardParticles.at(iMother));
                hardParticles.at(iMother).AddDaughter(&hardParticles.at(i));
            }
    }
}
