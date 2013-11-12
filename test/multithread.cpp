#include <GenericEventSelection.hpp>
#include <Dataset.hpp>
#include <PECReader.hpp>
#include <BTagger.hpp>
#include <BTagDatabase.hpp>
#include <TriggerSelection.hpp>
#include <WeightPileUp.hpp>
#include <RunManager.hpp>
#include <BasicKinematicsPlugin.hpp>

#include <iostream>


using namespace std;


int main()
{
    // Define the b-tagging objects
    BTagger bTagger(BTagger::Algorithm::CSV, BTagger::WorkingPoint::Medium);
    BTagDatabase bTagDatabase(bTagger);
    
    
    // Define the event selection
    GenericEventSelection sel(30., bTagger);
    sel.AddLeptonThreshold(Lepton::Flavour::Muon, 26.);
    sel.AddJetTagBin(2, 1);
    sel.AddJetTagBin(3, 1);
    sel.AddJetTagBin(3, 2);
    
    
    // Define datasets
    string const filePrefix("/afs/cern.ch/user/a/aapopov/workspace/data/2012Alpha/");
    double const brWlnu = 3 * 0.1080;
    list<Dataset> datasets;
    
    // t-channel single top
    datasets.emplace_back(Dataset::Process::ttchan, Dataset::Generator::POWHEG,
     Dataset::ShowerGenerator::Undefined);
    datasets.back().AddFile(filePrefix + "t-tchan-pw_rev468_QJd.root", 56.4 * brWlnu, 3915598);
    datasets.back().AddFile(filePrefix + "tbar-tchan-pw_rev468_koy.root", 30.7 * brWlnu, 1711403);
    //^ The SM x-sections are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
    
    // th with exotic coupling
    datasets.emplace_back(Dataset::Process::thqExotic, Dataset::Generator::MadGraph,
     Dataset::ShowerGenerator::Pythia);
    datasets.back().AddFile(filePrefix + "tqh-nc-mg-v2_rev468_Goi.root", 34.1e-3, 499958);
    
     
    // Define the triggers
    list<TriggerRange> triggerRanges;
    triggerRanges.emplace_back(190000, 210000, "HLT_IsoMu24_v", 20e3, "HLT_IsoMu24_v");
    
    TriggerSelection triggerSel(triggerRanges);
    
    
    // An object to reweight for pile-up
    WeightPileUp weigtPileUp("SingleMu2012ABCD_Alpha-v2_pixelLumi.pileupTruth_finebin.root", 0.06);
    
    
    // Construct the run manager
    RunManager manager(datasets.begin(), datasets.end());
    
    // Set the configuration for PECReader
    auto &config = manager.GetPECReaderConfig();
    config.SetModule(&triggerSel);
    config.SetModule(&sel);
    config.SetModule(&bTagger);
    config.SetModule(&bTagDatabase);
    config.SetModule(&weigtPileUp);
    
    // Register a plugin
    manager.RegisterPlugin(new BasicKinematicsPlugin("basicTuples"));
    
    // Process the datasets
    manager.Process(3);
    
    
    return 0;
}