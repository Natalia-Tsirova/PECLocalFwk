/**
 * \file PECReader.hpp
 * \author Andrey Popov, Natalia Tsirova
 * 
 * The module defines a class to read files in PlainEventContent (PEC) format.
 */

#pragma once

#include <PECReaderForward.hpp>
 
#include <PECReaderConfigForward.hpp>
#include <PhysicsObjects.hpp>
#include <Dataset.hpp>
#include <EventID.hpp>
#include <GenParticle.hpp>
#include <TriggerSelectionInterface.hpp>
#include <EventSelectionInterface.hpp>
#include <WeightBTagInterface.hpp>
#include <WeightPileUpInterface.hpp>
#include <SystDefinition.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector2.h>

#include <vector>
#include <list>
#include <string>
#include <memory>


/**
 * \class PECReader
 * \brief Class to read files in PlainEventContent (PEC) format
 * 
 * This is the core class of the package. It reads a set of files in PEC format [1] as specified by
 * an instance of class Dataset. It performs an event selection requested by the user with the help
 * of instances of classes TriggerSelectionInterface and EventSelectionInterface. For each event the
 * user is provided a set of collections of different physics objects form the event; the objects
 * are represented by dedicated classes described in module PhysicsObjects.
 * [1] https://twiki.cern.ch/twiki/bin/view/CMS/PlainEventContentTuples
 * 
 * In case of simulated events, reweighting for b-tagging scale factors and pile-up is also
 * performed with the help of dedicated classes.
 * 
 * Both trigger-bit and event selection as well as reweighting for b-tagging and pile-up are not
 * mandatory. If one of these modules is not provided, the class prints a warning and falls back to
 * a reasonable default behaviour.
 * 
 * The user can configure an instance of class PECReader using methods that modify one parameter at
 * a time or provide a complete (or partial) configuration described in an instance of class
 * PECReaderConfig. An instance of PECReader does not own the configuration parameters.
 * 
 * Quality criteria to identify physics objects are hard-coded in the class and are not expected to
 * be accessed by the user; instead, they are fixed to CMS-wide recommendations.
 * 
 * The class is non-copyable. No move constructor is implemented.
 */
class PECReader
{
public:
    /**
     * \brief Constructor from a dataset
     * 
     * The user must set configuration parameters afterwards.
     */
    PECReader(Dataset const &dataset);
    
    /**
     * \brief Constructor from a dataset and a configuration
     * 
     * A local copy of provided configuration is kept in the object. The user can modify the
     * copy later.
     */
    PECReader(Dataset const &dataset, PECReaderConfig const &config);
    
    /// Copy constructor is deleted
    PECReader(PECReader const &) = delete;
    
    /// Move constructor is deleted
    PECReader(PECReader &&) = delete;
    
    /// Assignment operator is deleted
    PECReader &operator=(PECReader const &) = delete;
    
    /// Destructor
    ~PECReader() = default;

public:
    /// Configures this from a configuration object
    void Configure(PECReaderConfig const &config);
    
    /// Sets the trigger selection
    void SetTriggerSelection(TriggerSelectionInterface *triggerSelection);
    
    /// Sets event selection
    void SetEventSelection(EventSelectionInterface const *eventSelection);
    
    /**
     * \brief Sets b-tagging configuration
     * 
     * This action has no effect on real data.
     */
    void SetBTagReweighter(WeightBTagInterface *bTagReweighter);
    
    /**
     * \brief Sets object to reweight simulation for pile-up
     * 
     * This action has no effect on real data.
     */
    void SetPileUpReweighter(WeightPileUpInterface const *puReweigher);
    
    /**
     * \brief Specifies whether information about the hard interaction is to be read
     * 
     * This action has no effect on real data.
     */
    void SetReadHardInteraction(bool flag = true);
    
    /**
     * \brief Specifies whether information on generator-level jets is to be read
     * 
     * This action has no effect on real data.
     */
    void SetReadGenJets(bool flag = true);
    
    /**
     * \brief Specifies whether information on parton shower should be read
     * 
     * This action has no effect on real data.
     */
    void SetReadPartonShower(bool flag = true) noexcept;
    
    /**
     * \brief Sets desired systematical variation
     * 
     * From the point of view of processing, there are two groups of sources of systematics.
     * Sources from the first group affect event weights only. They are calculated altogether
     * when the user requests systematical variation of type SystTypeAlgo::WeightOnly. Argument
     * direction is meaningless in this case and must be set to 0.
     * 
     * Sources from the second group change shapes of unweighted distributions (JEC uncertainty
     * is an example). Only one variation of such type can be evaluated at a time, and class
     * PECReader should be run several times. Parameter direction must equal +1 or (-1) to
     * choose "up" or "down" variation.
     * 
     * The user can instruct the class not to calculate any systematical variation by providing
     * type SystTypeAlgo::None.
     * 
     * Enumeration SystTypeAlgo is defined in file SystDefinition.hpp.
     */
    void SetSystematics(SystTypeAlgo type, int direction = 0);
    
    /// See documentation for SetSystematics(SystType, int)
    void SetSystematics(SystVariation const &syst);
    
    /**
     * \brief Opens a new file in the dataset
     * 
     * Opens the next file in the dataset for reading. Returns true in case of success and false
     * otherwise (when there are no more files in the dataset).
     */
    bool NextSourceFile();
    
    /**
     * \brief Reads the next event
     * 
     * Reads the next event that pass the event selection from the source files. If no such
     * event is found in current file, returns false, true otherwise.
     */
    bool NextEvent();
    
    /// Returns ID of the current event
    EventID const &GetEventID() const;
    
    /**
     * \brief Returns a list of tight leptons in the current event
     * 
     * The thresholds on transverse momenta is set the same as for loose leptons.
     */
    std::vector<Lepton> const &GetLeptons() const;
    
    /**
     * \brief Returns analysis-level jets in the current event
     * 
     * These juts meet requirements in eventSelection::IsAnalysisJet method.
     */
    std::vector<Jet> const &GetJets() const;
    
    /**
     * \brief Returns additional jets in the current events
     * 
     * These jets fail requirements of eventSelection::IsAnalysisJet method. Normally, they
     * are moderately soft jets needed for some observables.
     */
    std::vector<Jet> const &GetAdditionalJets() const;
    
    /// Returns MET
    Candidate const &GetMET() const;
    
    /**
     * \brief Returns reconstructed neutrino
     * 
     * The neutrino is reconstructed under the hypothesis that it originates from W-boson decay.
     * The accompanying charged lepton is identified with the leading tight lepton. The code
     * will crash if there are no tight leptons. Transverse component of neutrino momentum is
     * not affected by the reconstruction and is exactly the same as returned by GetMET method.
     * 
     * \note The method will turn obsolete in near future.
     */
    Candidate const &GetNeutrino() const;
    
    /// Returns number of reconstructed primary vertices (size of "offlinePrimaryVertices")
    unsigned GetNPrimaryVertices() const;
    
    /// Returns average angular energy density (rho)
    double GetRho() const;
        
    /**
     * \brief Returns central weight for the current event
     * 
     * In case of real data it is always 1.
     */
    double GetCentralWeight() const;
    
    /**
     * \brief Returns the XSec weight for the current event
     * 
     */
    double GetXSecWeight() const;
    
    /**
     * \brief Returns the btag weight for the current event
     * 
     */
    double GetBTagWeight() const;
    
    /**
     * \brief Returns the trigger weight for the current event
     * 
     */
    double GetTriggerWeight() const;
    
    /**
     * \brief Returns the PU weight for the current event
     * 
     */
    double GetPUWeight() const;
    
    /**
     * \brief Returns systematical variations of event weight for a specified source
     * 
     * For most of uncertainty sources the vector contains a single pair, but for some types
     * it aggregates several statistically independent variations (e.g. for PDF). The weights
     * should be used as is (no need to multiply them by the central weight). If variations for
     * a specified source cannot be calculated (for example, if the user has not provided the
     * corresponding module), the returned vector is empty.
     * 
     * If the user has not instructed this to calculate the varied weights by calling of
     * SetSystematics, the method throws an exception.
     */
    std::vector<WeightPair> const &GetSystWeight(SystTypeWeight type) const;
    
    /// Returns generator-level particles involved in the hard interaction
    std::vector<GenParticle> const &GetHardGenParticles() const;
    
    /// Returns generator-level jets
    std::vector<GenJet> const &GetGenJets() const;
    
    /**
     * \brief Returns partons from parton shower
     * 
     * Note that usually PEC samples keep only heavy-flavour partons if any. The method throws an
     * exception if called for real data.
     */
    std::vector<ShowerParton> const &GetShowerPartons() const;

private:
    /**
     * \brief Verifies that this is properly configured and performs final initializations
     * 
     * Motivation for this method is the fact that an instance of PECReader might not be fully
     * configured at construction time as the user might change some of parameters afterwards.
     * This method is run before the first file is opened.
     */
    void Initialize();
    
    /**
     * \brief Prepares to read the current ROOT file
     * 
     * Opens the current ROOT file, assigns the buffers to read the trees' branches, initializes
     * the event counters.
     */
    void OpenSourceFile();
    
    /**
     * \brief Closes the current ROOT file
     */
    void CloseSourceFile();
    
    /**
     * \brief Performs the event selection
     * 
     * The method performs the event selection and builds physical objects to be used by
     * plugins (jets, leptons, neutrino).
     */
    bool BuildAndSelectEvent();
    
    /**
     * \brief Calculate event weights (including systematics)
     * 
     * The method is called for both simulation and real data. However, only trigger weight is
     * evaluated in the latter case; it is needed to allow an additional event selection via
     * TriggerRange::PassEventSelection.
     */
    void CalculateEventWeights();
    
    /// Stores particles from the hard interaction in hardParticles collection
    void ParseHardInteraction();
    
    /// Stores generator-level jets in the dedicated vector
    void BuildGenJets();
    
    /// Reads information on parton shower and stores the partons in the dedicated vector
    void ReadPartonShower();

private:
    /// A copy of dataset to be processed
    Dataset const dataset;
    
    /// Specifies whether the object is fully configured
    bool isInitialized;
    
    /// Pointer to an object to perform trigger selection
    TriggerSelectionInterface *triggerSelection;
    
    /// Pointer to an object to perform event selection
    EventSelectionInterface const *eventSelection;
    
    /// An object to reweight event for b-tagging scale factors
    WeightBTagInterface *bTagReweighter;
    
    /// Pointer to an object to perform pile-up reweighting (set for a simulation dataset only)
    WeightPileUpInterface const *puReweighter;
    
    /// A short-cut for PECReaderConfig::GetReadHardInteraction
    bool readHardParticles;
    
    /// Indicates whether generator-level jets should be read
    bool readGenJets;
    
    /// Indicates whether information on parton shower should be read
    bool readPartonShower;
    
    
    /// Systematical variation
    SystVariation syst;
    
    
    /// Central event weight (as opposed to systematical variations)
    double weightCentral;
    
    /// Weight due to cross-section
    double weightCrossSection;
    
    /// Weight due to btag
    double weightBTag;
    
    /// Weight due to trigger
    double weightTrigger;
    
    /// Weight due to PU
    double weightPU;
    
    /// Systematical variations in event weight due to uncertainty in pile-up
    std::vector<WeightPair> systWeightPileUp;
    
    /// Systematical variations in event weight due to uncertainty in b-tagging tag rate
    std::vector<WeightPair> systWeightTagRate;
    
    /// Systematical variations in event weight due to uncertainty in b-tagging mistag rate
    std::vector<WeightPair> systWeightMistagRate;
    
    
    /// Iterator to the current Dataset::File object
    std::list<Dataset::File>::const_iterator sourceFileIt;
    
    
    TFile *sourceFile;   ///< The current source file
    TTree *eventIDTree;  ///< The tree with the event ID information
    TTree *triggerTree;  ///< The tree with the trigger information
    TTree *generalTree;  ///< The tree with all the information but triggers and event ID
    unsigned long nEventsTree;  ///< The total number of events in the trees
    unsigned long curEventTree;  ///< The index of the current event in the trees
    EventID eventID;  ///< An aggregate to store the event ID
    
    /// Maximal length to allocate buffers to read trees
    static unsigned const maxSize = 64;
    
    // Input buffers
    ULong64_t runNumber, lumiSection, eventNumber;
    
    UChar_t eleSize;
    Float_t elePt[maxSize];
    Float_t eleEta[maxSize];
    Float_t elePhi[maxSize];
    Float_t eleRelIso[maxSize];
    Float_t eleDB[maxSize];
    Bool_t eleTriggerPreselection[maxSize];
    Float_t eleMVAID[maxSize];
    Bool_t elePassConversion[maxSize];
    Bool_t eleQuality[maxSize];
    Bool_t eleCharge[maxSize];
    
    UChar_t muSize;
    Float_t muPt[maxSize];
    Float_t muEta[maxSize];
    Float_t muPhi[maxSize];
    Float_t muRelIso[maxSize];
    Float_t muDB[maxSize];
    Bool_t muQualityTight[maxSize];
    Bool_t muCharge[maxSize];
    
    UChar_t jetSize;
    Float_t jetPt[maxSize];
    Float_t jetEta[maxSize];
    Float_t jetPhi[maxSize];
    Float_t jetMass[maxSize];
    Float_t jetCSV[maxSize];
    Float_t jetTCHP[maxSize];
    Char_t jetFlavour[maxSize];
    Float_t jecUncertainty[maxSize];
    Float_t jerFactor[maxSize];
    Float_t jetCharge[maxSize];
    Float_t jetPullAngle[maxSize];
    
    /*
    Float_t softJetPt;
    Float_t softJetEta;
    Float_t softJetPhi;
    Float_t softJetMass;
    Float_t softJetHt;
    Float_t softJetPtJECUnc;
    Float_t softJetEtaJECUnc;
    Float_t softJetPhiJECUnc;
    Float_t softJetMassJECUnc;
    Float_t softJetHtJECUnc;
    */
    
    UChar_t metSize;
    Float_t metPt[maxSize];
    Float_t metPhi[maxSize];
    
    Short_t processID;  // needed to split the inclusive W+jets
    
    // Buffers to read the hard interaction
    UChar_t hardPartSize;
    Char_t hardPartPdgId[maxSize];
    Char_t hardPartFirstMother[maxSize], hardPartLastMother[maxSize];
    Float_t hardPartPt[maxSize];
    Float_t hardPartEta[maxSize];
    Float_t hardPartPhi[maxSize];
    Float_t hardPartMass[maxSize];
    
    
    // Buffers to read generator jets
    UChar_t genJetSize;
    Float_t genJetPt[maxSize], genJetEta[maxSize], genJetPhi[maxSize], genJetMass[maxSize];
    //UChar_t genJetBMult[maxSize], genJetCMult[maxSize];
    
    
    // Buffers for information on parton shower. Consult documentation for [1] for details
    //[1] https://github.com/andrey-popov/single-top/blob/master/plugins/PartonShowerOutcome.h
    // Number of partons stored
    UChar_t psSize;
    
    // PDG ID of partons
    Short_t psPdgId[maxSize];
    
    // Origin of partons
    UChar_t psOrigin[maxSize];
    
    /// Three-momentum of partons
    Float_t psPt[maxSize], psEta[maxSize], psPhi[maxSize];
    
    
    // Pile-up truth information
    
    // Number of recontructed primary vertices
    UChar_t pvSize;
    
    // "True" number of pile-up interactions (available in simulation only)
    Float_t puTrueNumInteractions;
    
    // Mean angular energy density
    Float_t puRho;
    
    
    
    Int_t nWeight_PDF;
    Float_t weight_PDFUp[maxSize], weight_PDFDown[maxSize];
    
    
    // The compact event description
    
    /// The tight leptons. Normally there is only one
    std::vector<Lepton> tightLeptons;
    
    /// The loose leptons
    std::vector<Lepton> looseLeptons;
    
    /// The selected jets to be used in the analysis. Normally, they have pt > 30 GeV/c
    std::vector<Jet> goodJets;
    
    /// The selected soft jets. Normally, they have 20 < pt < 30 GeV/c
    std::vector<Jet> additionalJets;
    
    /// MET of the current event
    Candidate correctedMET;
    
    /// The reconstructed neutrino
    Candidate neutrino;
    
    /// The generator particles from the hard interaction
    std::vector<GenParticle> hardParticles;
    
    /// The generator-level jets
    std::vector<GenJet> genJets;
    
    /**
     * \brief Partons from parton shower
     * 
     * They are read if only readPartonShower flag is set to true.
     */
    std::vector<ShowerParton> psPartons;
};
