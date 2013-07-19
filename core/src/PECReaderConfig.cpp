#include <PECReaderConfig.hpp>

#include <stdexcept>


using namespace std;


PECReaderConfig::PECReaderConfig():
    readHardInteraction(false)
{}


PECReaderConfig::PECReaderConfig(PECReaderConfig const &src):
    triggerSelection(src.triggerSelection->Clone()),
    eventSelection(src.eventSelection->Clone()),
    bTagger(src.bTagger),
    bTagDatabase(new BTagDatabase(*src.bTagDatabase)),
    puReweighter(src.puReweighter->Clone()),
    weightFilesLocation(src.weightFilesLocation),
    readHardInteraction(src.readHardInteraction),
    syst(src.syst)
{}


void PECReaderConfig::SetTriggerSelection(TriggerSelectionInterface const *triggerSelection_)
{
    triggerSelection.reset(triggerSelection_->Clone());
}


void PECReaderConfig::SetModule(TriggerSelectionInterface const *triggerSelection)
{
    SetTriggerSelection(triggerSelection);
}


void PECReaderConfig::SetEventSelection(EventSelectionInterface const *eventSelection_)
{
    eventSelection.reset(eventSelection_->Clone());
}


void PECReaderConfig::SetModule(EventSelectionInterface const *eventSelection)
{
    SetEventSelection(eventSelection);
}


void PECReaderConfig::SetBTagger(BTagger const *bTagger_)
{
    bTagger.reset(new BTagger(*bTagger_));
}


void PECReaderConfig::SetModule(BTagger const *bTagger)
{
    SetBTagger(bTagger);
}


void PECReaderConfig::SetBTagDatabase(BTagDatabase const *bTagDatabase_)
{
    bTagDatabase.reset(new BTagDatabase(*bTagDatabase_));
}


void PECReaderConfig::SetModule(BTagDatabase const *bTagDatabase)
{
    SetBTagDatabase(bTagDatabase);
}


void PECReaderConfig::SetPileUpReweighter(WeightPileUpInterface const *puReweighter_)
{
    puReweighter.reset(puReweighter_->Clone());
}


void PECReaderConfig::SetModule(WeightPileUpInterface const *puReweighter)
{
    SetPileUpReweighter(puReweighter);
}


void PECReaderConfig::SetWeightFilesLocation(string const &weightFilesLocation_)
{
    weightFilesLocation = weightFilesLocation_;
    
    // Make sure the last symbol is a slash
    if (weightFilesLocation.back() not_eq '/')
        weightFilesLocation += '/';
}


void PECReaderConfig::SetReadHardInteraction(bool readHardInteraction_)
{
    readHardInteraction = readHardInteraction_;
}


void PECReaderConfig::SetSystematics(SystVariation const &syst_)
{
    syst = syst_;
}


bool PECReaderConfig::IsSetTriggerSelection() const
{
    return bool(triggerSelection);
}


TriggerSelectionInterface const *PECReaderConfig::GetTriggerSelection() const
{
    if (not triggerSelection)
        throw logic_error("PECReaderConfig::GetTriggerSelection: Accessing an undefined "
         "configuration parameter.");
    
    return triggerSelection.get();
}


bool PECReaderConfig::IsSetEventSelection() const
{
    return bool(eventSelection);
}


EventSelectionInterface const *PECReaderConfig::GetEventSelection() const
{
    if (not eventSelection)
        throw logic_error("PECReaderConfig::GetEventSelection: Accessing an undefined "
         "configuration parameter.");
    
    return eventSelection.get();
}


bool PECReaderConfig::IsSetBTagger() const
{
    return bool(bTagger);
}


BTagger const *PECReaderConfig::GetBTagger() const
{
    if (not bTagger)
        throw logic_error("PECReaderConfig::GetBTagger: Accessing an undefined "
         "configuration parameter.");
    
    return bTagger.get();
}


bool PECReaderConfig::IsSetBTagDatabase() const
{
    return bool(bTagDatabase);
}


BTagDatabase const *PECReaderConfig::GetBTagDatabase() const
{
    if (not bTagDatabase)
        throw logic_error("PECReaderConfig::GetBTagDatabase: Accessing an undefined "
         "configuration parameter.");
    
    return bTagDatabase.get();
}


BTagDatabase *PECReaderConfig::GetBTagDatabase()
{
    if (not bTagDatabase)
        throw logic_error("PECReaderConfig::GetBTagDatabase: Accessing an undefined "
         "configuration parameter.");
    
    return bTagDatabase.get();
}


bool PECReaderConfig::IsSetPileUpReweighter() const
{
    return bool(puReweighter);
}


WeightPileUpInterface const *PECReaderConfig::GetPileUpReweighter() const
{
    if (not puReweighter)
        throw logic_error("PECReaderConfig::GetPileUpReweigher: Accessing an undefined "
         "configuration parameter.");
    
    return puReweighter.get();
}


WeightPileUpInterface *PECReaderConfig::GetPileUpReweighter()
{
    if (not puReweighter)
        throw logic_error("PECReaderConfig::GetPileUpReweigher: Accessing an undefined "
         "configuration parameter.");
    
    return puReweighter.get();
}


bool PECReaderConfig::IsSetWeightFilesLocation() const
{
    return (weightFilesLocation.length() not_eq 0);
}


string const &PECReaderConfig::GetWeightFilesLocation() const
{
    if (not IsSetWeightFilesLocation())
        throw logic_error("PECReaderConfig::GetWeightFilesLocation: Accessing an undefined "
         "configuration parameter.");
    
    return weightFilesLocation;
}


bool PECReaderConfig::GetReadHardInteraction() const
{
    return readHardInteraction;
}


SystVariation const &PECReaderConfig::GetSystematics() const
{
    return syst;
}
