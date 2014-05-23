//-------------------------------------------------------------------
// Clean up merged files from lxbtch
//
// execute with:
// root -l -q cleanUpMergedFiles(_infile_, _outfile_)
//
// Jay Lawhorn 11/4/13
//-------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
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
#include "Math/LorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void cleanUpMergedFiles(TString infilename="/afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_temp.root",
			  TString outfilename="/afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2/Working/LL-4p-0-100-v1510_14TEV_clean.root") {

  // set up input/output variables and file
  UInt_t nEvents;
  Float_t eventWeight;
  UInt_t tauDecayCat1, tauDecayCat2;
  UInt_t bTag1, bTag2;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;

  TFile* infile = new TFile(infilename); assert(infile);
  TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("eventWeight",    &eventWeight);
  intree->SetBranchAddress("bTag1",          &bTag1);
  intree->SetBranchAddress("bTag2",          &bTag2);
  intree->SetBranchAddress("genB1",          &genB1);
  intree->SetBranchAddress("genB2",          &genB2);
  intree->SetBranchAddress("recoB1",         &recoB1);
  intree->SetBranchAddress("recoB2",         &recoB2);
  intree->SetBranchAddress("tauDecayCat1",   &tauDecayCat1);
  intree->SetBranchAddress("tauDecayCat2",   &tauDecayCat2);
  intree->SetBranchAddress("genTau1",        &genTau1);
  intree->SetBranchAddress("genTau2",        &genTau2);
  intree->SetBranchAddress("genDecayTau1",   &genDecayTau1);
  intree->SetBranchAddress("genDecayTau2",   &genDecayTau2);
  intree->SetBranchAddress("recoTau1",       &recoTau1);
  intree->SetBranchAddress("recoTau2",       &recoTau2);

  TTree* infotree = (TTree*) infile->Get("Info"); assert(infotree);
  infotree->SetBranchAddress("nEvents",      &nEvents);

  Int_t totalEvents=0;

  for (UInt_t iEntry=0; iEntry<infotree->GetEntries(); iEntry++) {
    infotree->GetEntry(iEntry);
    totalEvents+=nEvents;
  }

  TFile *outFile = new TFile(outfilename, "RECREATE");

  // tree to hold information about selected events
  TTree *outTree = new TTree("Events", "Events");
  outTree->Branch("eventWeight",    &eventWeight,    "eventWeight/f");  // event weight from cross-section and Event->Weight
  outTree->Branch("tauDecayCat1",   &tauDecayCat1,   "tauDecayCat1/i"); // leading tau final state - jet, muon, electron
  outTree->Branch("tauDecayCat2",   &tauDecayCat2,   "tauDecayCat2/i"); // second tau final state - jet, muon, electron
  outTree->Branch("bTag1",          &bTag1,          "bTag1/i");        // leading b-jet tag from delphes
  outTree->Branch("bTag2",          &bTag2,          "bTag2/i");        // second b-jet tag from delphes
  outTree->Branch("genTau1",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau1);      // 4-vector for generator leading tau
  outTree->Branch("genTau2",        "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genTau2);      // 4-vector for generator second tau
  outTree->Branch("genDecayTau1",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genDecayTau1); // 4-vector for generator decay product of leading tau
  outTree->Branch("genDecayTau2",   "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genDecayTau2); // 4-vector for generator decay product of second tau
  outTree->Branch("recoTau1",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau1);     // 4-vector for reconstructed leading tau
  outTree->Branch("recoTau2",       "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoTau2);     // 4-vector for reconstructed second tau
  outTree->Branch("genB1",          "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genB1);        // 4-vector for generator leading b-quark
  outTree->Branch("genB2",          "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genB2);        // 4-vector for generator second b-quark
  outTree->Branch("recoB1",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB1);       // 4-vector for reconstructed leading b-jet
  outTree->Branch("recoB2",         "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &recoB2);       // 4-vector for reconstructed second b-jet

  for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // entry loop
    intree->GetEntry(iEntry);

    eventWeight=eventWeight/Float_t(totalEvents);

    outTree->Fill();

  }
  outFile->Write();
  outFile->Close();


}
