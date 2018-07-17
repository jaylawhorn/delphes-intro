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
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TLorentzVector.h>

#include "HttStyles.h"

#include "helperFxns.hh"
#include "drawFxns.hh"

using namespace std;

//not great practice but ok
const Double_t MUON_MASS = 0.105658369;
const Double_t ELE_MASS =  0.000511;
const Double_t W_MASS =   80.385;
const Double_t Z_MASS =   91.1876;
const Float_t MAX_MATCH_DIST = 0.4;
const Float_t LUMI=3000;

Float_t selection(const TString inputfile, HistHolder hists);

Int_t findZ( TClonesArray * array, int &zll_1, int &zll_2, float massRange);

void analysisWWZ() {

  HistHolder signal; initHistHolder(signal, "sig");
  HistHolder ttbar;  initHistHolder(ttbar, "tt");
  HistHolder otherb; initHistHolder(otherb, "bgd");

  float wwz=selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/WWZ_hadd.root", signal);
  float ttz=selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/TTZ_hadd.root", ttbar);
  float tt =selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/TT_hadd.root", ttbar);
  float www=selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/WWW_hadd.root", otherb);
  float wzz=selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/WZZ_hadd.root", otherb);
  float zzz=selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/ZZZ_hadd.root", otherb);
  float ww =selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/WW_hadd.root", otherb);
  float zz =selection("/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-test/ZZTo2L2Q_hadd.root", otherb);

  float bkgds = ttz+tt+www+wzz+zzz+ww+zz;

  cout << "======= Event yields =======" << endl;
  cout << "  WWZ: " << wwz<< endl;
  cout << "Total: " << bkgds<< endl;
  cout << "Backgrounds" << endl;
  cout << "  TTZ: " << ttz<< endl;
  cout << "  TT : " << tt<< endl;
  cout << "  WWW: " << www<< endl;
  cout << "  WZZ: " << wzz<< endl;
  cout << "  ZZZ: " << zzz<< endl;
  cout << "  WW : " << ww<< endl;
  cout << "  ZZ : " << zz<< endl;
  
  cout << endl;

  cout << "S/B = " << wwz / bkgds << endl;
  cout << "S/sqrt(B) = " << wwz / sqrt(bkgds) << endl;

  cout << "============================" << endl;


  TCanvas *c = MakeCanvas("c","c",800, 600);
  gStyle->SetOptStat(0);

  //Drawing plots (look at drawFxns.hh)
  DrawHists(c, "zee.png", signal.zee, ttbar.zee, otherb.zee);
  DrawHists(c, "zmm.png", signal.zmm, ttbar.zmm, otherb.zmm);

  DrawHists(c, "zeeL.png", signal.zeeL, ttbar.zeeL, otherb.zeeL);
  DrawHists(c, "zeeS.png", signal.zeeS, ttbar.zeeS, otherb.zeeS);
  DrawHists(c, "zmmL.png", signal.zmmL, ttbar.zmmL, otherb.zmmL);
  DrawHists(c, "zmmS.png", signal.zmmS, ttbar.zmmS, otherb.zmmS);

  DrawHists(c, "nE.png", signal.nE, ttbar.nE, otherb.nE);
  DrawHists(c, "nM.png", signal.nM, ttbar.nM, otherb.nM);
  DrawHists(c, "nJ.png", signal.nJ, ttbar.nJ, otherb.nJ);

}

Float_t selection(const TString inputfile, HistHolder hists) {

  const Double_t Z_MASS_RANGE = 15;


  Float_t eventWeight;
  uint nEl, nMu, nJet, nBJet;
  Float_t metPt;
  Float_t metPhi;

  TClonesArray *vEl  = new TClonesArray("TLepton");
  TClonesArray *vMu  = new TClonesArray("TLepton");
  TClonesArray *vJet = new TClonesArray("THadJet");

  TFile *inFile = new TFile(inputfile, "READ");

  TH1D *hTot = (TH1D*) inFile->Get("hTot");
  float nTot = hTot->GetBinContent(1);

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

  float nCount=0;
  
  for (uint iEntry=0; iEntry<inTree->GetEntries(); iEntry++) {
    inTree->GetEntry(iEntry);

    float evtWeight = (eventWeight * LUMI) / nTot;

    if (nBJet>0) continue;

    int nEle = vEl->GetEntries();
    int nMu  = vMu->GetEntries();
    int nJet = vJet->GetEntries();

    if (nEle+nMu<4) continue;

    // try to make a Z
    int zee1=-1, zee2=-1;
    int nZee = findZ(vEl, zee1, zee2, Z_MASS_RANGE);
    int zmm1=-1, zmm2=-1;
    int nZmm = findZ(vMu, zmm1, zmm2, Z_MASS_RANGE);

    if (nZee+nZmm==0) continue;

    nCount+=evtWeight;

    hists.nE->Fill(nEle, evtWeight);
    hists.nM->Fill(nMu, evtWeight);
    hists.nJ->Fill(nJet, evtWeight);

    TLepton *lep1=0, *lep2=0;
    if (nZmm==1) {
      lep1 = (TLepton*) (*vMu)[zmm1];
      lep2 = (TLepton*) (*vMu)[zmm2];

      if (lep1->PT > lep2->PT) {
	hists.zmmL->Fill( lep1->PT, evtWeight);
	hists.zmmS->Fill( lep2->PT, evtWeight);
      }
      else {
	hists.zmmL->Fill( lep2->PT, evtWeight);
	hists.zmmS->Fill( lep1->PT, evtWeight);
      }

      TLorentzVector v1; v1. SetPtEtaPhiM(lep1->PT, lep1->Eta, lep1->Phi, MUON_MASS);
      TLorentzVector v2; v2. SetPtEtaPhiM(lep2->PT, lep2->Eta, lep2->Phi, MUON_MASS);

      hists.zmm->Fill( (v1+v2).M(), evtWeight);

    }
    else if (nZee==1) {
      lep1 = (TLepton*) (*vEl)[zee1];
      lep2 = (TLepton*) (*vEl)[zee2];

      if (lep1->PT > lep2->PT) {
	hists.zeeL->Fill( lep1->PT, evtWeight);
	hists.zeeS->Fill( lep2->PT, evtWeight);
      }
      else {
	hists.zeeL->Fill( lep2->PT, evtWeight);
	hists.zeeS->Fill( lep1->PT, evtWeight);
      }
      
      TLorentzVector v1; v1. SetPtEtaPhiM(lep1->PT, lep1->Eta, lep1->Phi, ELE_MASS);
      TLorentzVector v2; v2. SetPtEtaPhiM(lep2->PT, lep2->Eta, lep2->Phi, ELE_MASS);

      hists.zee->Fill( (v1+v2).M(), evtWeight);
    }

  }

  return nCount;

}

Int_t findZ( TClonesArray * vLep, int &zll_1, int &zll_2, float massRange=15) {

  zll_1=-1; zll_2=-1;

  if (vLep->GetEntries()>=2) {

    for (int i=0; i<vLep->GetEntries(); i++) {
      TLepton *lep1 = (TLepton*) (*vLep)[i];
      TLorentzVector v1; v1.SetPtEtaPhiM(lep1->PT, lep1->Eta, lep1->Phi, (lep1->PID == 11) ? ELE_MASS : MUON_MASS);

      for (int j=i+1; j<vLep->GetEntries(); j++) {
	TLepton *lep2 = (TLepton*) (*vLep)[j];
	TLorentzVector v2; v2.SetPtEtaPhiM(lep2->PT, lep2->Eta, lep2->Phi, (lep1->PID == 11) ? ELE_MASS : MUON_MASS);

	if (lep1->Charge == lep2->Charge) continue;

	if ( abs( (v1+v2).M() - Z_MASS ) < massRange ) {
	  zll_1=i; zll_2=j;
	  break;
	}
      }
      if (zll_1 != zll_2) break;
    }
  }

  if (zll_1 == zll_2) return 0;
  else return 1;


}
