#include <PhysicsObjects.hpp>


// Methods of class Candidate
Candidate::Candidate()
{}


Candidate::Candidate(TLorentzVector const &p4_):
    p4(p4_)
{}


void Candidate::SetP4(TLorentzVector const &p4_)
{
    p4 = p4_;
}


void Candidate::SetPtEtaPhiM(double pt, double eta, double phi, double mass)
{
    p4.SetPtEtaPhiM(pt, eta, phi, mass);
}


void Candidate::SetPxPyPzE(double px, double py, double pz, double E)
{
    p4.SetPxPyPzE(px, py, pz, E);
}


TLorentzVector const &Candidate::P4() const
{
    return p4;
}


double Candidate::Pt() const
{
    return p4.Pt();
}


double Candidate::Eta() const
{
    return p4.Eta();
}


double Candidate::Phi() const
{
    return p4.Phi();
}


double Candidate::M() const
{
    return p4.M();
}


bool Candidate::operator<(Candidate const &rhs) const
{
    return (p4.Pt() < rhs.p4.Pt());
}


// Methods of class Lepton
Lepton::Lepton():
    Candidate(),
    flavour(Flavour::Unknown),
    relIso(-1.), dB(0.), charge(0)
{}


Lepton::Lepton(Lepton::Flavour flavour_, TLorentzVector const &p4):
    Candidate(p4),
    flavour(flavour_),
    relIso(-1.), dB(0.), charge(0)
{}


void Lepton::SetRelIso(double relIso_)
{
    relIso = relIso_;
}


void Lepton::SetDB(double dB_)
{
    dB = dB_;
}


void Lepton::SetCharge(int charge_)
{
    charge = charge_;
}


Lepton::Flavour Lepton::GetFlavour() const
{
    return flavour;
}


double Lepton::RelIso() const
{
    return relIso;
}


double Lepton::DB() const
{
    return dB;
}


int Lepton::Charge() const
{
    return charge;
}


// Methods of class Jet
Jet::Jet():
    Candidate(),
    CSVValue(-1e6), TCHPValue(-1e6),
    parentPDGID(0)
{}


Jet::Jet(TLorentzVector const &p4):
    Candidate(p4),
    CSVValue(-1e6), TCHPValue(-1e6),
    parentPDGID(0)
{}


void Jet::SetBTags(double CSV, double JP, double TCHP)
{
    CSVValue = CSV;
    JPValue = JP;
    TCHPValue = TCHP;
}


void Jet::SetCSV(double CSV)
{
    CSVValue = CSV;
}


void Jet::SetJP(double JP)
{
    JPValue = JP;
}


void Jet::SetTCHP(double TCHP)
{
    TCHPValue = TCHP;
}


void Jet::SetParentID(int pdgID)
{
    parentPDGID = pdgID;
}


void Jet::SetPassPileUpID(bool passesPileUpID_)
{
    passesPileUpID = passesPileUpID_;
}


double Jet::CSV() const
{
    return CSVValue;
}


double Jet::JP() const
{
    return JPValue;
}


double Jet::TCHP() const
{
    return TCHPValue;
}


int Jet::GetParentID() const
{
    return parentPDGID;
}


bool Jet::PassPileUpID() const
{
    return passesPileUpID;
}