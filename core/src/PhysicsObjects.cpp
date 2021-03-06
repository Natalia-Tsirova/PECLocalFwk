#include <PhysicsObjects.hpp>

#include <limits>


using namespace std;


// Methods of class Candidate
Candidate::Candidate() noexcept
{}


Candidate::Candidate(TLorentzVector const &p4_) noexcept:
    p4(p4_)
{}


void Candidate::SetP4(TLorentzVector const &p4_) noexcept
{
    p4 = p4_;
}


void Candidate::SetPtEtaPhiM(double pt, double eta, double phi, double mass) noexcept
{
    p4.SetPtEtaPhiM(pt, eta, phi, mass);
}


void Candidate::SetPxPyPzE(double px, double py, double pz, double E) noexcept
{
    p4.SetPxPyPzE(px, py, pz, E);
}


TLorentzVector const &Candidate::P4() const noexcept
{
    return p4;
}


double Candidate::Pt() const noexcept
{
    return p4.Pt();
}


double Candidate::Eta() const noexcept
{
    return p4.Eta();
}


double Candidate::Phi() const noexcept
{
    return p4.Phi();
}


double Candidate::M() const noexcept
{
    return p4.M();
}


bool Candidate::operator<(Candidate const &rhs) const noexcept
{
    return (p4.Pt() < rhs.p4.Pt());
}


// Methods of class Lepton
Lepton::Lepton() noexcept:
    Candidate(),
    flavour(Flavour::Unknown),
    relIso(-1.), dB(0.), charge(0)
{}


Lepton::Lepton(Lepton::Flavour flavour_, TLorentzVector const &p4) noexcept:
    Candidate(p4),
    flavour(flavour_),
    relIso(-1.), dB(0.), charge(0)
{}


void Lepton::SetRelIso(double relIso_) noexcept
{
    relIso = relIso_;
}


void Lepton::SetDB(double dB_) noexcept
{
    dB = dB_;
}


void Lepton::SetCharge(int charge_) noexcept
{
    charge = charge_;
}


Lepton::Flavour Lepton::GetFlavour() const noexcept
{
    return flavour;
}


double Lepton::RelIso() const noexcept
{
    return relIso;
}


double Lepton::DB() const noexcept
{
    return dB;
}


int Lepton::Charge() const noexcept
{
    return charge;
}


// Methods of class Jet
Jet::Jet() noexcept:
    Candidate(),
    CSVValue(-numeric_limits<double>::infinity()),
    JPValue(-numeric_limits<double>::infinity()),
    TCHPValue(-numeric_limits<double>::infinity()),
    parentPDGID(0),
    charge(-10.), pullAngle(-10.)
{}


Jet::Jet(TLorentzVector const &p4) noexcept:
    Candidate(p4),
    CSVValue(-numeric_limits<double>::infinity()),
    JPValue(-numeric_limits<double>::infinity()),
    TCHPValue(-numeric_limits<double>::infinity()),
    parentPDGID(0),
    charge(-10.), pullAngle(-10.)
{}


void Jet::SetBTags(double CSV, double JP, double TCHP) noexcept
{
    CSVValue = CSV;
    JPValue = JP;
    TCHPValue = TCHP;
}


void Jet::SetCSV(double CSV) noexcept
{
    CSVValue = CSV;
}


void Jet::SetJP(double JP) noexcept
{
    JPValue = JP;
}


void Jet::SetTCHP(double TCHP) noexcept
{
    TCHPValue = TCHP;
}


void Jet::SetParentID(int pdgID) noexcept
{
    parentPDGID = pdgID;
}


void Jet::SetCharge(double charge_) noexcept
{
    charge = charge_;
}


void Jet::SetPullAngle(double pullAngle_) noexcept
{
    pullAngle = pullAngle_;
}


double Jet::CSV() const noexcept
{
    return CSVValue;
}


double Jet::JP() const noexcept
{
    return JPValue;
}


double Jet::TCHP() const noexcept
{
    return TCHPValue;
}


int Jet::GetParentID() const noexcept
{
    return parentPDGID;
}


double Jet::Charge() const noexcept
{
    return charge;
}


double Jet::GetPullAngle() const noexcept
{
    return pullAngle;
}


// Methods of class GenJet
GenJet::GenJet() noexcept:
    Candidate(),
    bMult(0), cMult(0)
{}


GenJet::GenJet(TLorentzVector const &p4) noexcept:
    Candidate(p4),
    bMult(0), cMult(0)
{}


void GenJet::SetMultiplicities(unsigned bMult_, unsigned cMult_) noexcept
{
    bMult = bMult_;
    cMult = cMult_;
}


unsigned GenJet::GetBMultiplicity() const noexcept
{
    return bMult;
}


unsigned GenJet::GetCMultiplicity() const noexcept
{
    return cMult;
}


// Methods of the ShowerParton class
ShowerParton::ShowerParton() noexcept:
    Candidate(),
    pdgId(0), origin(Origin::Unknown)
{}


ShowerParton::ShowerParton(TLorentzVector const &p4_, int pdgId_,
 Origin origin_ /*= Origin::Unknown*/) noexcept:
    Candidate(p4_),
    pdgId(pdgId_), origin(origin_)
{}


ShowerParton::ShowerParton(double pt, double eta, double phi, int pdgId_,
 Origin origin_ /*= Origin::Unknown*/) noexcept:
    Candidate(),
    pdgId(pdgId_), origin(origin_)
{
    SetPtEtaPhiM(pt, eta, phi, GuessMass(pdgId));
}


void ShowerParton::SetOrigin(Origin origin_) noexcept
{
    origin = origin_;
}


ShowerParton::Origin ShowerParton::GetOrigin() const noexcept
{
    return origin;
}


void ShowerParton::SetPdgId(int pdgId_) noexcept
{
    pdgId = pdgId_;
}


int ShowerParton::GetPdgId() const noexcept
{
    return pdgId;
}


double ShowerParton::GuessMass(int pdgId) noexcept
{
    // Masses for s, c, b are set to values used in Pythia in Summer12 datasets
    switch (abs(pdgId))
    {
        case 6:
            return 172.5;
        
        case 5:
            return 4.8;
        
        case 4:
            return 1.5;
        
        case 3:
            return 0.5;
        
        case 2:
            return 0.;
        
        case 1:
            return 0.;
        
        default:
            return 0.;
    }
}
