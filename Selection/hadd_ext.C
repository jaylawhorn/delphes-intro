#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TLorentzVector.h"

using namespace std;

class TLepton : public TObject {
public:
  TLepton():
    PT(0), Eta(0), Phi(0), Charge(0), PID(0), IsolationVar(0)
  {}
  TLepton(float pt, float eta, float phi, int charge, int pid, float iso):
    PT(pt), Eta(eta), Phi(phi), Charge(charge), PID(pid), IsolationVar(iso)
  {}
  ~TLepton(){}

  float PT; float Eta; float Phi; int Charge; int PID; float IsolationVar;
};

class THadJet : public TObject {
public:
  THadJet():
    PT(0), Eta(0), Phi(0), Mass(0), BTag(0)
  {}
  THadJet(float pt, float eta, float phi, float mass, int btag):
    PT(pt), Eta(eta), Phi(phi), Mass(mass), BTag(btag)
  {}
  ~THadJet(){}

  float PT; float Eta; float Phi; float Mass; int BTag;
};

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 );


void hadd_ext(const TString inputfile="/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-2l-4j/TTTT.root") {

  Float_t eventWeight;
  Float_t eventWeightN;
  TFile *inFile = new TFile(inputfile, "update");

  TH1D *hTot = (TH1D*) inFile->Get("hTot");
  float nTot = hTot->GetBinContent(1);

  TTree *inTree = (TTree*) inFile->Get("Events");
  inTree->SetBranchAddress("eventWeight", &eventWeight);
  TBranch *bnew = inTree->Branch("eventWeightN", &eventWeightN, "eventWeightN/F");
  //inTree->Branch("eventWeightN", &eventWeightN, "eventWeightN/F");

  for (int i=0; i<inTree->GetEntries(); i++) {
    inTree->GetEntry(i);
    eventWeightN = eventWeight / nTot;
    bnew->Fill();
  }  

  inTree->Write();

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  return phiDiff;

}
