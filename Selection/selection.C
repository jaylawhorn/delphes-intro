#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
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

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 );

void selection(const TString inputfile="/eos/cms/store/group/upgrade/delphes_output/YR_Delphes/Delphes342pre15/ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_200PU/ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_9_0.root",
	       const Float_t xsec=40,
	       Int_t sampleNo=100,
	       const TString outputfile="test.root") {

  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  //const Double_t ELE_MASS  = 0.000511;
  //const Double_t TAU_MASS  = 1.77682;

  //const Int_t ELE_ID_CODE = 11;
  //const Int_t MUON_ID_CODE = 13;
  //const Int_t TAU_ID_CODE = 15;
  //
  //const Int_t B_ID_CODE = 5;
  //const Int_t G_ID_CODE = 22;

  //const Float_t MAX_MATCH_DIST = 0.4;

  // event categories
  //enum { ZEE=0, ZMM, ZJJ };

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

  //TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("MuonLoose");
  //TClonesArray *branchMET =treeReader->UseBranch("MissingET");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  //TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  //TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  //set up loop variables
  Muon *mu=0;
  GenParticle *genPart=0;
  LHEFEvent *event=0;

  Float_t eventWeight;

  Float_t ptMu1, etaMu1, phiMu1, mMu1;
  //Float_t ptGenMu1, etaGenMu1, phiGenMu1;

  Float_t ptMu2, etaMu2, phiMu2, mMu2;
  //Float_t ptGenMu2, etaGenMu2, phiGenMu2;

  Float_t mZ;//, mGenZ;

  TFile *outFile = new TFile(outputfile, "RECREATE");
  TH1D *h = new TH1D("hTot","hTot",1,0,1);
  h->SetBinContent(1,numberOfEntries);

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");

  outTree->Branch("eventWeight", &eventWeight, "eventWeight/f"); // event weight from cross-section and Event->Weight

  outTree->Branch("ptMu1",  &ptMu1,  "ptMu1/f");  // pt(mu1)
  outTree->Branch("etaMu1", &etaMu1, "etaMu1/f"); // eta(mu1)
  outTree->Branch("phiMu1", &phiMu1, "phiMu1/f"); // phi(mu1)
  outTree->Branch("mMu1",   &mMu1,   "mMu1/f");   // m(mu1)

  outTree->Branch("ptMu2",  &ptMu2,  "ptMu2/f");  // pt(mu2)
  outTree->Branch("etaMu2", &etaMu2, "etaMu2/f"); // eta(mu2)
  outTree->Branch("phiMu2", &phiMu2, "phiMu2/f"); // phi(mu2)
  outTree->Branch("mMu2",   &mMu2,   "mMu2/f");   // m(mu2)

  outTree->Branch("mZ",     &mZ,     "mZ/f");     // m(Zmm)

  //outTree->Branch("ptGenMu1",  &ptGenMu1,  "ptGenMu1/f");  // pt(mu1) (gen)
  //outTree->Branch("etaGenMu1", &etaGenMu1, "etaGenMu1/f"); // eta(mu1) (gen)
  //outTree->Branch("phiGenMu1", &phiGenMu1, "phiGenMu1/f"); // phi(mu1) (gen)
  //outTree->Branch("mGenMu1",   &mGenMu1,   "mGenMu1/f");   // m(mu1) (gen)
  //
  //outTree->Branch("ptGenMu2",  &ptGenMu2,  "ptGenMu2/f");  // pt(mu2) (gen)
  //outTree->Branch("etaGenMu2", &etaGenMu2, "etaGenMu2/f"); // eta(mu2) (gen)
  //outTree->Branch("phiGenMu2", &phiGenMu2, "phiGenMu2/f"); // phi(mu2) (gen)
  //outTree->Branch("mGenMu2",   &mGenMu2,   "mGenMu2/f");   // m(mu2) (gen)
  //
  //outTree->Branch("mGenZ",     &mGenZ,     "mGenZ/f");     // m(Zmm) (gen)

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);
    
    ptMu1=-999; etaMu1=-999; phiMu1=-999;
    //ptGenMu1=-999; etaGenMu1=-999; phiGenMu1=-999;
    
    ptMu2=-999; etaMu2=-999; phiMu2=-999;
    //ptGenMu2=-999; etaGenMu2=-999; phiGenMu2=-999;
    
    mZ=-999; //mGenZ=-999;

    // ********************
    // EVENT WEIGHT
    // ********************

    eventWeight = 1;
    //if (branchEvent && sampleNo!=100) {
    //event = (LHEFEvent*) branchEvent->At(0);
    //eventWeight*=event->Weight;
    //}
    eventWeight *= xsec;

    for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
      mu = (Muon*) branchMuon->At(iMuon);
      
      if (fabs(mu->Eta)>4.0) continue;
      if (mu->PT<10) continue;
      //if (mu->IsolationVar>0.4) continue;

      genPart = (GenParticle*) mu->Particle.GetObject();

      if(mu->PT>ptMu1) {
	ptMu2  = ptMu1;
	etaMu2 = etaMu1;
	phiMu2 = phiMu1;
	mMu2   = mMu1;

	ptMu1  = mu->PT;
	etaMu1 = mu->Eta;
	phiMu1 = mu->Phi;
	mMu1   = MUON_MASS;
      }
      else if (mu->PT>ptMu2) {
	ptMu2  = mu->PT;
        etaMu2 = mu->Eta;
	phiMu2 = mu->Phi;
        mMu2   = MUON_MASS;
      }
    }
    
    if (mMu1<0 || mMu2<0) continue;
    
    //std::cout << "hi" << std::endl;

    TLorentzVector m1; m1.SetPtEtaPhiM(ptMu1, etaMu1, phiMu1, mMu1);
    TLorentzVector m2; m2.SetPtEtaPhiM(ptMu2, etaMu2, phiMu2, mMu2);

    mZ = (m1+m2).M();

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
