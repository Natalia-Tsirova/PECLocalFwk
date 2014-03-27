#include <SingleTopTChanBTagEffPlugin_wjets.hpp>

#include <Processor.hpp>
#include <ROOTLock.hpp>

#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <sys/stat.h>


using namespace std;


SingleTopTChanBTagEffPlugin_wjets::SingleTopTChanBTagEffPlugin_wjets(string const &outDirectory_, string const hftype_str_):
    Plugin("SingleTopBTagEff"), outDirectory(outDirectory_), hftype_str(hftype_str_)
{
    // Make sure the directory path ends with a slash
    if (outDirectory.back() != '/')
        outDirectory += '/';
    
    // Create the output directory if it does not exist
    struct stat dirStat;
    
    if (stat(outDirectory.c_str(), &dirStat) != 0)  // the directory does not exist
        mkdir(outDirectory.c_str(), 0755);
    
    // Convert hftype_str to WjetsHF Type
    if (hftype_str == "W_qq") hftype = WjetsHFPlugin::Type::W_qq;
    else if (hftype_str == "W_c") hftype = WjetsHFPlugin::Type::W_c;
    else if (hftype_str == "W_other") hftype = WjetsHFPlugin::Type::W_other;
    else if (hftype_str == "W_light") hftype = WjetsHFPlugin::Type::W_light;
    else throw runtime_error("SingleTopTChanPlugin_wjets: Undefined heavy flavour type = " + hftype_str);
}


Plugin *SingleTopTChanBTagEffPlugin_wjets::Clone() const
{
    return new SingleTopTChanBTagEffPlugin_wjets(*this);
}


void SingleTopTChanBTagEffPlugin_wjets::BeginRun(Dataset const &dataset)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    WjetsHFClassifier = dynamic_cast<WjetsHFPlugin const *>(processor->GetPluginBefore("WjetsHF", name));
    
   // Creation of ROOT objects is not thread-safe and must be protected
    ROOTLock::Lock();
    
    // Create the output file
    file = new TFile((outDirectory + dataset.GetFiles().front().GetBaseName() + "_" + hftype_str + ".root").c_str(),
     "recreate");
    
    // Create histograms
    histTagB = new TH2D("histTagB", "", 10000, 30., 200., 10, -2.4, 2.4);
    histB = new TH2D("histB", "", 10000, 30., 200., 10, -2.4, 2.4);
    histTagC = new TH2D("histTagC", "", 10000, 30., 200., 10, -2.4, 2.4);
    histC = new TH2D("histC", "", 10000, 30., 200., 10, -2.4, 2.4);
    histTagUDS = new TH2D("histTagUDS", "", 10000, 30., 200., 10, -2.4, 2.4);
    histUDS = new TH2D("histUDS", "", 10000, 30., 200., 10, -2.4, 2.4);
    histTagG = new TH2D("histTagG", "", 10000, 30., 200., 10, -2.4, 2.4);
    histG = new TH2D("histG", "", 10000, 30., 200., 10, -2.4, 2.4);
    
    histB->Sumw2();
    histTagB->Sumw2();
    histC->Sumw2();
    histTagC->Sumw2();
    histUDS->Sumw2();
    histTagUDS->Sumw2();
    histG->Sumw2();
    histTagG->Sumw2();
    
    // End of critical block
    ROOTLock::Unlock();
}


void SingleTopTChanBTagEffPlugin_wjets::EndRun()
{
    // Operations with ROOT objects performed here are not thread-safe and must be guarded
    ROOTLock::Lock();
    
    // Write the tree and close the file
    file->cd();
    TH2D EffB(*histTagB);
    EffB.SetName("EffB");
    EffB.Divide(histB);
    TH2D EffC(*histTagC);
    EffC.SetName("EffC");
    EffC.Divide(histC);
    TH2D EffUDS(*histTagUDS);
    EffUDS.SetName("EffUDS");
    EffUDS.Divide(histUDS);
    TH2D EffG(*histTagG);
    EffG.SetName("EffG");
    EffG.Divide(histG);
    EffB.Write("", TObject::kOverwrite);
    EffC.Write("", TObject::kOverwrite);
    EffUDS.Write("", TObject::kOverwrite);
    EffG.Write("", TObject::kOverwrite);
    histTagB->Write("", TObject::kOverwrite);
    histB->Write("", TObject::kOverwrite);
    histTagC->Write("", TObject::kOverwrite);
    histC->Write("", TObject::kOverwrite);
    histTagUDS->Write("", TObject::kOverwrite);
    histUDS->Write("", TObject::kOverwrite);
    histTagG->Write("", TObject::kOverwrite);
    histG->Write("", TObject::kOverwrite);
    
    // Delete the objects
    /*delete histL;
    delete histTagL;
    delete histC;
    delete histTagC;
    delete histB;
    delete histTagB;
    
    delete file;*/
    
    ROOTLock::Unlock();
}


bool SingleTopTChanBTagEffPlugin_wjets::ProcessEvent()
{
    // Classify the Wjets event
    if (WjetsHFClassifier->GetDecision() != hftype) return false;
    
    // Make sure the event contains reasonable physics objects
    if ((*reader)->GetLeptons().size() not_eq 1 or (*reader)->GetJets().size() < 2)
        return false;
    
    // Save event ID
    auto const &eventID = (*reader)->GetEventID();
    runNumber = eventID.Run();
    eventNumber = eventID.Event();
    lumiSection = eventID.LumiBlock();
    
    weight = (*reader)->GetCentralWeight();
    
    // Define some short-cuts
    auto const &jets = (*reader)->GetJets();
    
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        if (abs(jets.at(i).GetParentID()) == 5)
        {
            histB->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
            
            if (jets.at(i).CSV() > 0.898)
                histTagB->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
        }
        
        else if (abs(jets.at(i).GetParentID()) == 4)
        {
            histC->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
            
            if (jets.at(i).CSV() > 0.898)
                histTagC->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
        }
        
        else if (abs(jets.at(i).GetParentID()) < 4)
        {
            histUDS->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
            
            if (jets.at(i).CSV() > 0.898)
                histTagUDS->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
        }
        else if (jets.at(i).GetParentID() == 21)
        {
            histG->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
            
            if (jets.at(i).CSV() > 0.898)
                histTagG->Fill(jets.at(i).Pt(),jets.at(i).Eta(), weight);
        }
        
    }

    return true;
}
