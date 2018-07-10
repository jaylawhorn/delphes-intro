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

using namespace std;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

Float_t deltaPhi( const Float_t phi1, const Float_t phi2 );

void analysis(const TString inputfile="/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/ZZZ/ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_10_0.root") {

  Float_t eventWeight;

  Float_t ptMu1, etaMu1, phiMu1, mMu1;

  Float_t ptMu2, etaMu2, phiMu2, mMu2;

  Float_t mZ;

  TFile *inFile = new TFile(inputfile, "READ");
  TTree *inTree = (TTree*) inFile->Get("Events");
  inTree->SetBranchAddress("eventWeight", &eventWeight);
  inTree->SetBranchAddress("ptMu1",  &ptMu1);
  inTree->SetBranchAddress("etaMu1", &etaMu1);
  inTree->SetBranchAddress("phiMu1", &phiMu1);
  inTree->SetBranchAddress("mMu1",   &mMu1);
  inTree->SetBranchAddress("ptMu2",  &ptMu2);
  inTree->SetBranchAddress("etaMu2", &etaMu2);
  inTree->SetBranchAddress("phiMu2", &phiMu2);
  inTree->SetBranchAddress("mMu2",   &mMu2);
  inTree->SetBranchAddress("mZ",     &mZ);


  for (int i=0; i<inTree->GetEntries(); i++) {
    inTree->GetEntry(i);

    // continue as normal
    
  }  

  
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
