#include <FilterEventIDPlugin.hpp>

#include <Processor.hpp>

#include <stdexcept>
#include <fstream>
#include <sstream>


using namespace std;


FilterEventIDPlugin::FilterEventIDPlugin(string const &name_, string const &eventIDsFileName,
 bool rejectKnownEvent_ /*= true*/):
    Plugin(name_),
    rejectKnownEvent(rejectKnownEvent_),
    eventIDsCurFile(nullptr)
{
    // Open the text file with event IDs
    ifstream eventIDsFile(eventIDsFileName);
    
    if (not eventIDsFile.good())
        throw runtime_error(string("FilterEventIDPlugin::FilterEventIDPlugin: Cannot open file \"" +
         eventIDsFileName + "\"."));
    
    
    // Read IDs from the file
    string line;
    bool stop = false;
    
    while (not stop)
    {
        // Read until a file name is found
        while (true)
        {
            getline(eventIDsFile, line);
            
            if (line.find("# Name of the file") not_eq string::npos)
                break;
            
            if (eventIDsFile.eof())  // no more files mentioned in eventIDsFile
            {
                stop = true;
                break;
            }
        }
        
        if (stop)
            // Skip the rest of the loop
            continue;
        
        // The next line is the file name
        getline(eventIDsFile, line);
        auto &eventIDsCurFile = eventIDsAllFiles[line];
        
        
        // Skip two lines and read the number of events
        for (unsigned i = 0; i < 3; ++i)
            getline(eventIDsFile, line);
        
        istringstream ist(line);
        unsigned long nEntries;
        ist >> nEntries;
        
        eventIDsCurFile.reserve(nEntries);
        
        
        // Skip two lines
        getline(eventIDsFile, line);
        getline(eventIDsFile, line);
        
        
        // Read IDs
        while (true)
        {
            getline(eventIDsFile, line);
            
            if (eventIDsFile.eof() or line.length() == 0)
                break;
            
            ist.clear();
            ist.str(line);  // should contain three numbers separated by semicolon
            
            unsigned long run, lumiSection, event;
            ist >> run;
            ist.ignore();
            ist >> lumiSection;
            ist.ignore();
            ist >> event;
            
            eventIDsCurFile.emplace_back(run, lumiSection, event);
        }
    }
    
    eventIDsFile.close();
}


FilterEventIDPlugin::FilterEventIDPlugin(string const &name_,
 map<string, vector<EventID>> const &eventIDsAllFiles_, bool rejectKnownEvent_):
    Plugin(name_),
    rejectKnownEvent(rejectKnownEvent_),
    eventIDsAllFiles(eventIDsAllFiles_),
    eventIDsCurFile(nullptr)
{}


Plugin *FilterEventIDPlugin::Clone() const
{
    return new FilterEventIDPlugin(name, eventIDsAllFiles, rejectKnownEvent);
}


void FilterEventIDPlugin::BeginRun(Dataset const &dataset)
{
    // Lists of event IDs are provided on per-file basis. Since the plugin is not notified when
    //a new file in the current dataset is started, it makes sense only to filer atomic datasets
    if (dataset.GetFiles().size() not_eq 1)
        throw logic_error("TTbarRecoPlugin::BeginRun: The plugin can filter atomic datasets only.");
    
    
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    
    
    // Make a short-cut for list of event IDs for the new atomic dataset
    string const fileName(dataset.GetFiles().front().name);
    string const shortFileName = fileName.substr(fileName.find_last_of('/') + 1);
    
    auto const resIt = eventIDsAllFiles.find(shortFileName);
    eventIDsCurFile = (resIt == eventIDsAllFiles.end()) ? nullptr : &resIt->second;
}


void FilterEventIDPlugin::EndRun()
{}


bool FilterEventIDPlugin::ProcessEvent()
{
    // Check if there is a list of event IDs for the current ROOT file
    if (eventIDsCurFile == nullptr)
        return not rejectKnownEvent;
    
    // Search for ID of the current event in the vector
    auto const &id = (*reader)->GetEventID();
    bool const eventFound = (find(eventIDsCurFile->begin(), eventIDsCurFile->end(),
     EventID(id.Run(), id.LumiBlock(), id.Event())) not_eq eventIDsCurFile->end());
    
    
    return (rejectKnownEvent) ? not eventFound : eventFound;
}