#include <BTagEfficiencies.hpp>

#include <FileInPath.hpp>
#include <ROOTLock.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept>


using namespace std;


BTagEfficiencies::BTagEfficiencies(string const &fileName, string const &directory /*= ""*/):
    BTagEffInterface(),
    inFileDirectory(directory)
{
    // Guard creation of a ROOT file
    ROOTLock::Lock();
    
    // Resolve path to the source file and open it. If the file is missing, pathBuilder will throw
    //an exception
    FileInPath pathBuilder;
    TFile *fp = TFile::Open(pathBuilder.Resolve("BTag", fileName).c_str());
    srcFile.reset(fp);
    
    ROOTLock::Unlock();
    
    
    // Make sure the in-file directory path is either empty or terminates with a slash
    if (inFileDirectory.length() > 0 and inFileDirectory[inFileDirectory.length() - 1] != '/')
        inFileDirectory += '/';
}


BTagEfficiencies::BTagEfficiencies(BTagEfficiencies const &src):
    BTagEffInterface(src),
    srcFile(src.srcFile),
    inFileDirectory(src.inFileDirectory),
    processLabels(src.processLabels),
    processMap(src.processMap),
    defaultProcessLabel(src.defaultProcessLabel),
    effHists(src.effHists)
{}


BTagEfficiencies::BTagEfficiencies(BTagEfficiencies &&src) noexcept:
    BTagEffInterface(move(src)),
    srcFile(move(src.srcFile)),
    inFileDirectory(move(src.inFileDirectory)),
    processLabels(move(src.processLabels)),
    processMap(move(src.processMap)),
    defaultProcessLabel(move(src.defaultProcessLabel)),
    effHists(move(src.effHists))
{}


BTagEfficiencies &BTagEfficiencies::operator=(BTagEfficiencies const &rhs)
{
    BTagEffInterface::operator=(rhs);
    
    srcFile = rhs.srcFile;
    inFileDirectory = rhs.inFileDirectory;
    processLabels = rhs.processLabels;
    processMap = rhs.processMap;
    defaultProcessLabel = rhs.defaultProcessLabel;
    effHists = rhs.effHists;
    
    return *this;
}


BTagEfficiencies::~BTagEfficiencies() noexcept
{}


void BTagEfficiencies::SetProcessLabel(Dataset::Process code, string const &label)
{
    // Try to find the given process code among already registered ones
    auto mapRuleIt = find_if(processMap.begin(), processMap.end(),
     [code](decltype(*processMap.cbegin()) &rule){return (rule.first == code);});
    
    
    // Try to find the given label among already registered ones
    auto labelIt = find(processLabels.begin(), processLabels.end(), label);
    int labelIndex = distance(processLabels.begin(), labelIt);
    
    // Add label to the container if it was not found there. After this operation labelIndex refers
    //to the correct position
    if (labelIt == processLabels.end())
        processLabels.push_back(label);
    
    
    // Check if the given process code is a new one
    if (mapRuleIt == processMap.end())
    {
        // Add a new map rule
        processMap.emplace_back(code, labelIndex);
    }
    else
    {
        // Reset label registered for the process code
        mapRuleIt->second = labelIndex;
    }
}


void BTagEfficiencies::SetProcessLabel(list<Dataset::Process> const &codes, string const &label)
{
    for (auto const &code: codes)
        SetProcessLabel(code, label);
}


void BTagEfficiencies::SetDefaultProcessLabel(string const &label)
{
    defaultProcessLabel = label;
}


BTagEffInterface *BTagEfficiencies::Clone() const
{
    return new BTagEfficiencies(*this);
}


void BTagEfficiencies::LoadPayload(Dataset const &dataset)
{
    // Clear the map with efficiency histograms
    effHists.clear();
    
    
    // Find the text label corresponding to the new process
    string curProcessLabel;
    auto mapRuleIt = find_if(processMap.begin(), processMap.end(),
     [&dataset](decltype(*processMap.cbegin()) &rule){return dataset.TestProcess(rule.first);});
    
    if (mapRuleIt != processMap.end())
        curProcessLabel = processLabels.at(mapRuleIt->second);
    else
    {
        // No label is defined for this dataset. Check if the default one is available
        if (defaultProcessLabel.length() == 0)
        {
            ostringstream ost;
            ost << "BTagEfficiencies::LoadPayload: Cannot find which histogram with b-tagging " <<
             "efficiencies should be used with the dataset containing file \"" <<
             dataset.GetFiles().front().GetBaseName() << ".root\". " <<
             "No rule to define the histogram name is available for the dataset, and no default " <<
             "name is specified.";
            
            throw runtime_error(ost.str());
        }
        
        // Apparently, there is a default label. Use it
        curProcessLabel = defaultProcessLabel;
    }
    
    
    // Loop over possible working points
    for (auto const &wp: {BTagger::WorkingPoint::Tight, BTagger::WorkingPoint::Medium,
     BTagger::WorkingPoint::Loose})
    {
        string const wpCode(WorkingPointToText(wp));
        
        
        // Read histograms for all jet flavours. This is not a thread-safe operation
        ROOTLock::Lock();
        
        shared_ptr<TH2> bHist(dynamic_cast<TH2 *>(srcFile->Get(
         (inFileDirectory + curProcessLabel + "_b_" + wpCode).c_str())));
        shared_ptr<TH2> cHist(dynamic_cast<TH2 *>(srcFile->Get(
         (inFileDirectory + curProcessLabel + "_c_" + wpCode).c_str())));
        shared_ptr<TH2> udsHist(dynamic_cast<TH2 *>(srcFile->Get(
         (inFileDirectory + curProcessLabel + "_uds_" + wpCode).c_str())));
        shared_ptr<TH2> gHist(dynamic_cast<TH2 *>(srcFile->Get(
         (inFileDirectory + curProcessLabel + "_g_" + wpCode).c_str())));
        
        // Make sure the histograms are not associated with a file
        for (auto const &p: {bHist, cHist, udsHist, gHist})
        {
            if (p)
                p->SetDirectory(nullptr);
        }
        
        ROOTLock::Unlock();
        
        
        // Add the histograms to the histogram map if the pointers are not null
        if (bHist)
            effHists[make_pair(wp, 5)] = bHist;
        
        if (cHist)
            effHists[make_pair(wp, 4)] = cHist;
        
        if (udsHist)
        {
            effHists[make_pair(wp, 1)] = udsHist;
            effHists[make_pair(wp, 2)] = udsHist;
            effHists[make_pair(wp, 3)] = udsHist;
        }
        
        if (gHist)
        {
            effHists[make_pair(wp, 21)] = gHist;
            effHists[make_pair(wp, 0)] = gHist;
            //^ Jets with unidentified flavour are considered along with gluons. The motivation is
            //that pile-up jets usually obtain a flavour of either 0 or 21
        }
    }
}


double BTagEfficiencies::GetEfficiency(BTagger::WorkingPoint wp, Candidate const &jet, int flavour)
 const
{
    // Find the appropriate efficiency histogram
    auto histIt = effHists.find(make_pair(wp, abs(flavour)));
    
    if (histIt == effHists.end())
    {
        ostringstream ost;
        ost << "BTagEfficiencies::GetEfficiency: Failed to find an efficiency histogram for " <<
         "working point " << int(wp) << " and jet flavour " << flavour << ".";
        
        throw runtime_error(ost.str());
    }
    
    
    // Find bin that contains the jet
    int const bin = histIt->second->FindFixBin(jet.Pt(), fabs(jet.Eta()));
    
        
    return histIt->second->GetBinContent(bin);
}


string BTagEfficiencies::WorkingPointToText(BTagger::WorkingPoint wp)
{
    switch (wp)
    {
        case BTagger::WorkingPoint::Tight:
            return "T";
        
        case BTagger::WorkingPoint::Medium:
            return "M";
        
        case BTagger::WorkingPoint::Loose:
            return "L";
        
        default:
        {
            ostringstream ost;
            ost << "BTagEfficiencies::WorkingPointToText: Unknown working point " << int(wp) <<
             " is specified.";
            
            throw runtime_error(ost.str());
        }
    }
}
