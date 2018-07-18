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
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TLorentzVector.h"
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

using namespace std;

class TLepton : public TObject {
public:
  TLepton(float pt, float eta, float phi, int charge, int pid, float iso): 
    PT(pt), Eta(eta), Phi(phi), 
    Charge(charge), PID(pid), IsolationVar(iso)
  {}
  ~TLepton(){}

  float PT;
  float Eta;
  float Phi;
  int Charge;
  int PID;
  float IsolationVar;

};

class THadJet : public TObject {
public:
  THadJet(float pt, float eta, float phi, float mass, int btag):
    PT(pt), Eta(eta), Phi(phi), 
    Mass(mass), BTag(btag)
  {}
  ~THadJet(){}
  
  float PT;
  float Eta;
  float Phi;
  float Mass;
  int BTag;
  
};


Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 );

void selectionWWZ(const TString inputfile="/eos/cms/store/group/upgrade/delphes_output/YR_Delphes/Delphes342pre15/WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_200PU/WWZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_10_0.root",
		  const Float_t xsec=40,
		  Int_t sampleNo=100,
		  const TString outputfile="test.root") {
  
  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;
  
  const Float_t MAX_MATCH_DIST = 0.4;
  
  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  if (!(branchJet)) {
    cout << "File broken" << endl;
    return;
  }
  
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("MuonLoose");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  
  //set up loop variables
  Muon *mu=0;
  Electron *ele=0;
  Jet *jet=0;
  MissingET *met=0;
  LHEFEvent *event=0;
  
  Float_t eventWeight;
  int nEl, nMu, nJet, nBJet;
  Float_t metPt;
  Float_t metPhi;
  
  TClonesArray *vEl  = new TClonesArray("TLepton"); TClonesArray &pEl = *vEl;
  TClonesArray *vMu  = new TClonesArray("TLepton"); TClonesArray &pMu = *vMu;
  TClonesArray *vJet = new TClonesArray("THadJet"); TClonesArray &pJet = *vJet;
  
  TFile *outFile = new TFile(outputfile, "RECREATE");
  TH1D *h = new TH1D("hTot","hTot",1,0,1);
  h->SetBinContent(1,numberOfEntries);
  
  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  
  outTree->Branch("eventWeight", &eventWeight, "eventWeight/f"); // event weight from cross-section and Event->Weight
  outTree->Branch("nEl",   &nEl,   "nEl/i");
  outTree->Branch("nMu",   &nMu,   "nMu/i");
  outTree->Branch("nJet",  &nJet,  "nJet/i");
  outTree->Branch("nBJet", &nBJet, "nBJet/i");
  outTree->Branch("metPt", &metPt, "metPt/f");
  outTree->Branch("metPhi", &metPhi, "metPhi/f");
  
  outTree->Branch("vEl", &vEl);
  outTree->Branch("vMu", &vMu);
  outTree->Branch("vJet", &vJet);
  
  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);
    pEl.Clear(); pMu.Clear(); 
    pJet.Clear();

    eventWeight = 1;
    //if (branchEvent && sampleNo!=100) {
    //event = (LHEFEvent*) branchEvent->At(0);
    //eventWeight*=event->Weight;
    //}
    eventWeight *= xsec;

    nMu=0;

    for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
      mu = (Muon*) branchMuon->At(iMuon);
      
      if (fabs(mu->Eta)>3.0) continue;
      if (mu->PT<10) continue;
      //if (mu->IsolationVar>0.2) continue;

      new(pMu[nMu]) TLepton(mu->PT, mu->Eta, mu->Phi, mu->Charge, 13, mu->IsolationVar);

      nMu++;
    }
    

    nEl=0;

    for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco electron loop
      ele = (Electron*) branchElectron->At(iEle);
      
      if (fabs(ele->Eta)>3.0) continue;
      if (ele->PT<10) continue;
      //if (ele->IsolationVar>0.2) continue;

      new(pEl[nEl]) TLepton(ele->PT, ele->Eta, ele->Phi, ele->Charge, 11, ele->IsolationVar);

      nEl++;
    }

    if (nEl+nMu<2) continue;

    nJet=0; nBJet=0;

    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reco jet loop
      jet = (Jet*) branchJet->At(iJet);
      
      if (fabs(jet->Eta)>4.0) continue;
      if (jet->PT<15) continue;

      int isOverlap=0;
      for (int iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) {
	mu = (Muon*) branchMuon->At(iMuon);

	if (fabs(mu->Eta)>3.0) continue;
	if (mu->PT<10) continue;
	if (mu->IsolationVar>0.2) continue;
	if (deltaR(mu->Eta, jet->Eta, mu->Phi, jet->Phi)<MAX_MATCH_DIST) {
	  isOverlap=1;
	  break;
	}
      }

      if (isOverlap==1) continue;

      for (int iEle=0; iEle<branchElectron->GetEntries(); iEle++) {
	ele = (Electron*) branchElectron->At(iEle);

	if (fabs(ele->Eta)>3.0) continue;
	if (ele->PT<10) continue;
	if (ele->IsolationVar>0.2) continue;
	if (deltaR(ele->Eta, jet->Eta, ele->Phi, jet->Phi)<MAX_MATCH_DIST) {
	  isOverlap=1;
	  break;
	}
      }

      if (isOverlap==1) continue;

      new(pJet[nJet]) THadJet(jet->PT, jet->Eta, jet->Phi, jet->Mass, jet->BTag);

      nJet++;

      if (jet->BTag==1) nBJet++;

    }

    met = (MissingET*) branchMET->At(0);

    metPt = met->MET;
    metPhi = met->Phi;
   
    outTree->Fill();
    
  }  

  outFile->Write();
  outFile->Save();
  
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
