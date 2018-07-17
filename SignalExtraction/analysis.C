#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
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


void analysis(const TString inputfile="/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/ZZZ_hadd.root") {

  Float_t eventWeight;
  uint nEl, nMu, nJet, nBJet;
  Float_t metPt;
  Float_t metPhi;

  TClonesArray *vEl  = new TClonesArray("TLepton");
  TClonesArray *vMu  = new TClonesArray("TLepton");
  TClonesArray *vJet = new TClonesArray("THadJet");

  TFile *inFile = new TFile(inputfile, "READ");
  TTree *inTree = (TTree*) inFile->Get("Events");
  inTree->SetBranchAddress("eventWeight", &eventWeight);
  inTree->SetBranchAddress("nEl", &nEl);
  inTree->SetBranchAddress("nMu", &nMu);
  inTree->SetBranchAddress("nJet", &nJet);
  inTree->SetBranchAddress("nBJet", &nBJet);
  inTree->SetBranchAddress("metPt", &metPt);
  inTree->SetBranchAddress("metPhi", &metPhi);

  inTree->SetBranchAddress("vEl", &vEl);
  inTree->SetBranchAddress("vMu", &vMu);
  inTree->SetBranchAddress("vJet", &vJet);

  inTree->GetEntry(0);

  //      no spaces  here V or here V
  TH1D *hNLep = new TH1D("hNLep", "hNLep", 15, 0, 15);
  TH1D *hNJet = new TH1D("hNJet", "hNJet", 15, 0, 15);
  TH2D *hNvN  = new TH2D("hNvN",  "hNvN", 15, 0, 15, 15, 0, 15);

  TH1D *hPtEle1 = new TH1D("hPtEle1", "hPtEle1", 25, 0, 250);

  for (int i=0; i<inTree->GetEntries(); i++) {
    inTree->GetEntry(i);

    hNLep->Fill(vEl->GetEntries()+vMu->GetEntries());
    hNJet->Fill(vJet->GetEntries());
    hNvN->Fill(vEl->GetEntries()+vMu->GetEntries(), vJet->GetEntries());

    TLepton *ele1 = (TLepton*) (*vEl)[0];
    hPtEle1->Fill(ele1->PT);
  }  

  TCanvas *c = new TCanvas("c","c",800, 600);

  hNLep->Draw("");
  c->SaveAs("nLep.png");

  hNJet->Draw("");
  c->SaveAs("nJet.png");

  hNvN->Draw("colz");
  c->SaveAs("nvn.png");

  hPtEle1->Draw("");
  c->SaveAs("ptEle1.png");
  
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
