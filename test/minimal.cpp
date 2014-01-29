#include <GenericEventSelection.hpp>
#include <Dataset.hpp>
#include <PECReader.hpp>
#include <BTagger.hpp>
#include <BTagDatabase.hpp>
#include <TriggerSelection.hpp>

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
    string const filePrefix("/afs/cern.ch/user/a/aapopov/workspace/data/2012Bravo/");
    double const brWlnu = 3 * 0.1080;
    list<Dataset> datasets;
    
    // ttbar
    datasets.emplace_back(Dataset::Process::ttbar, Dataset::Generator::MadGraph,
     Dataset::ShowerGenerator::Pythia);
    datasets.back().AddFile(filePrefix + "ttbar-semilep-mg-p1_53X.02.01_VNz_p1.root",
     234. * brWlnu * (1. - brWlnu) * 2, 24953451);
    datasets.back().AddFile(filePrefix + "ttbar-semilep-mg-p1_53X.02.01_VNz_p2.root",
     234. * brWlnu * (1. - brWlnu) * 2, 24953451);
    datasets.back().AddFile(filePrefix + "ttbar-semilep-mg-p1_53X.02.01_VNz_p3.root",
     234. * brWlnu * (1. - brWlnu) * 2, 24953451);
    datasets.back().AddFile(filePrefix + "ttbar-semilep-mg-p1_53X.02.01_VNz_p4.root",
     234. * brWlnu * (1. - brWlnu) * 2, 24953451);
    datasets.back().AddFile(filePrefix + "ttbar-dilep-mg_53X.02.01_FFe_p1.root",
     234. * brWlnu * brWlnu, 11991428);
    datasets.back().AddFile(filePrefix + "ttbar-dilep-mg_53X.02.01_FFe_p2.root",
     234. * brWlnu * brWlnu, 11991428);
    
    // t-channel single top
    datasets.emplace_back(Dataset::Process::ttchan, Dataset::Generator::POWHEG,
     Dataset::ShowerGenerator::Undefined);
    datasets.back().AddFile(filePrefix + "t-tchan-pw_53X.02.01_PIN.root ", 56.4 * brWlnu, 3915598);
    datasets.back().AddFile(filePrefix + "tbar-tchan-pw_53X.02.01_VcT.root ", 30.7 * brWlnu, 1711403);
    //^ The SM x-sections are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
    
    // tth
    datasets.emplace_back(Dataset::Process::tth, Dataset::Generator::Pythia,
     Dataset::ShowerGenerator::Undefined);
    datasets.back().AddFile(filePrefix + "tth_53X.02.01_bVJ.root", 0.1302, 995697);
    //^ The cross-section for tth is taken form https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV#ttH_Process
    
    // th with exotic coupling
    datasets.emplace_back(Dataset::Process::thqExotic, Dataset::Generator::MadGraph,
     Dataset::ShowerGenerator::Pythia);
    datasets.back().AddFile(filePrefix + "thq-nc-mg_53X.02.01_Wyg.root", 36.4e-3, 4847334);
    
     
    // Define the triggers
    list<TriggerRange> triggerRanges;
    triggerRanges.emplace_back(0, -1, "IsoMu24_eta2p1", 19.7e3, "IsoMu24_eta2p1");
    
    TriggerSelection triggerSel(triggerRanges);
    
    
    // Perform dataset-specific initialization
    bTagDatabase.SetDataset(datasets.front());
    
    
    // Build an instance of PECReader
    PECReader reader(datasets.front());
    reader.SetTriggerSelection(&triggerSel);
    reader.SetEventSelection(&sel);
    reader.SetBTaggingConfig(&bTagger, &bTagDatabase);
    //reader.SetReadHardInteraction();
    reader.NextSourceFile();
    
    
    // Loop over few events
    for (unsigned i = 0; i < 10; ++i)
    {
        cout << "Event " << i << '\n';
        reader.NextEvent();
        
        cout << "Tight leptons' pts:";
        
        for (auto const &l: reader.GetLeptons())
            cout << " pt: " << l.Pt() << ", iso: " << l.RelIso() << ", dB: " << l.DB();
        
        cout << "\nAnalysis jets' pts:";
        
        for (auto const &j: reader.GetJets())
            cout << " " << j.Pt();
        
        cout << "\nAdditional jets' pts:";
        
        for (auto const &j: reader.GetAdditionalJets())
            cout << " " << j.Pt();
        
        cout << "\nNeutrino's pt: " << reader.GetNeutrino().Pt() << ", pz: " <<
         reader.GetNeutrino().P4().Pz();
        
        /*
        cout << "\nHard generator-level particles:";
        
        for (auto const &p: reader.GetHardGenParticles())
        {
            cout << "\n PDG ID: " << p.GetPdgId() <<
             "\n # mothers: " << p.GetMothers().size() <<
             "\n # daughters: " << p.GetDaughters().size();
        }
        */
        
        cout << "\nEvent weight: " << reader.GetCentralWeight();
        
        cout << "\n\n";
    }
    
    
    return 0;
}