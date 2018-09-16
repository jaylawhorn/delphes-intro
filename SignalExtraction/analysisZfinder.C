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

//not great practice but ok
const Double_t MUON_MASS = 0.105658369;
const Double_t ELE_MASS =  0.000511;
const Double_t W_MASS =   80.385;
const Double_t Z_MASS =   91.1876;

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

void thework(const TString isamp, int cutlev);

Int_t findZ( TClonesArray * vLep, int &zll_1, int &zll_2, int &zll_3, int &zll_4, float massRange);
Int_t findZjets( TClonesArray * vLep, int &zll_1, int &zll_2, float massRange);

void analysisZfinder() {

  //TString samples[14] =  {"DYJets_comb", "Higgs_comb", "TT", "TTXX_comb", "TT_350", "TT_950", "VV_comb", "WWW", "WWWTo3L3Nu", "WWZ", "WWZTo4L2Nu", "WZZ", "ZZZ", "ssWW_comb"};
  TString samples[6] =  {"WWW", "WWWTo3L3Nu", "WWZ", "WWZTo4L2Nu", "WZZ", "ZZZ"};
  
  std::cout << " Pre-selection + trigger " << std::endl;
  std::cout << "\t $ 2L2J \t $ 2L4J \t $ 3L0J \t $ 3L2J \t $ 4L0J \t $ 4L2J \t $ 5L0J \t $ 6L0J \\\\" << std::endl;
  for (int i=0; i<6; i++) {
    thework(samples[i],0);
  }

  std::cout << " Pre-selection + trigger +b-jet veto" << std::endl;
  std::cout << "\t $ 2L2J \t $ 2L4J \t $ 3L0J \t $ 3L2J \t $ 4L0J \t $ 4L2J \t $ 5L0J \t $ 6L0J \\\\" << std::endl;
  for (int i=0; i<6; i++) {
    thework(samples[i],1);
  }

  std::cout << " + Z matching" << std::endl;
  std::cout << "\t $ 2L2J \t $ 2L4J \t $ 3L0J \t $ 3L2J \t $ 4L0J \t $ 4L2J \t $ 5L0J \t $ 6L0J \\\\" << std::endl;
  for (int i=0; i<6; i++) {
    thework(samples[i],2);
  }

}

void thework(const TString isamp,int cutlev) {
  const float lumi=3000.;
  //std::cerr << isamp << std::endl;

  const TString ifold="/afs/cern.ch/work/j/jlawhorn/public/delphes-sel-final-comb/";

  Float_t eventWeight;
  uint nEl, nMu, nJet, nBJet;
  Float_t metPt;
  Float_t metPhi;

  TClonesArray *vEl  = new TClonesArray("TLepton");
  TClonesArray *vMu  = new TClonesArray("TLepton");
  TClonesArray *vJet = new TClonesArray("THadJet");

  TFile *inFile = new TFile(ifold+isamp+".root", "READ");
  TTree *inTree = (TTree*) inFile->Get("Events");
  inTree->SetBranchAddress("eventWeightN", &eventWeight);
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

  float c2L2J=0, s2L2J=0;
  float c2L4J=0, s2L4J=0;
  float c3L0J=0, s3L0J=0;
  float c3L2J=0, s3L2J=0;
  float c4L0J=0, s4L0J=0;
  float c4L2J=0, s4L2J=0;
  float c5L0J=0, s5L0J=0;
  float c6L0J=0, s6L0J=0;

  for (int i=0; i<inTree->GetEntries(); i++) {
    inTree->GetEntry(i);

    int trigThresh=30;
    bool triggered=false;

    for (int j=0; j<nMu; j++) {
      if ( ((TLepton*)(*vMu)[j])->PT > trigThresh) {
	triggered=true;
	break;
      }
    }
    if (triggered==false) {
      for (int j=0; j<nEl; j++) {
	if ( ((TLepton*)(*vEl)[j])->PT > trigThresh) {
	  triggered=true;
	  break;
	}
      }
    }

    //if (triggered==false) continue;

    if (cutlev>=1 && nBJet!=0) continue;

    if (cutlev>=2) {
      const Double_t Z_MASS_RANGE = 15;

      int zee1=-1, zee2=-1, zee3=-1, zee4=-1;
      int nZee = findZ(vEl, zee1, zee2, zee3, zee4, Z_MASS_RANGE);
      int zmm1=-1, zmm2=-1, zmm3=-1, zmm4=-1;
      int nZmm = findZ(vMu, zmm1, zmm2, zmm3, zmm4, Z_MASS_RANGE);
      int zjj1=-1, zjj2=-1;
      int nZjj = findZjets(vJet, zjj1, zjj2, Z_MASS_RANGE);

      if (nZee+nZmm>=3) {
	c6L0J+=eventWeight*lumi;
	s6L0J+=eventWeight*eventWeight*lumi*lumi;
      }
      //else if (nZee+nZmm==2 && nEl+nMu>=5) {
      //c5L0J+=eventWeight*lumi;
      //s5L0J+=eventWeight*eventWeight*lumi*lumi;
      //}
      else if (nZee+nZmm==2 && nZjj==1) {
	c4L2J+=eventWeight*lumi;
	s4L2J+=eventWeight*eventWeight*lumi*lumi;
      }

    }
    else {

      if (nEl+nMu>=6) {
	c6L0J+=eventWeight*lumi;
	s6L0J+=eventWeight*eventWeight*lumi*lumi;
      }
      else if (nEl+nMu==5) {
	c5L0J+=eventWeight*lumi;
	s5L0J+=eventWeight*eventWeight*lumi*lumi;
      }
      else if (nEl+nMu==4) {
	if (nJet>=2) {
	  c4L2J+=eventWeight*lumi;
	  s4L2J+=eventWeight*eventWeight*lumi*lumi;
	}    
	else {
	  c4L0J+=eventWeight*lumi;
	  s4L0J+=eventWeight*eventWeight*lumi*lumi;
	}
      }
      else if (nEl+nMu==3) {
	if (nJet>=2) {
	  c3L2J+=eventWeight*lumi;
	  s3L2J+=eventWeight*eventWeight*lumi*lumi;
	}
	else {
	  c3L0J+=eventWeight*lumi;
	  s3L0J+=eventWeight*eventWeight*lumi*lumi;
	}
      }
      else if (nEl+nMu==2) {
	if ( nJet>=4) {
	  c2L4J+=eventWeight*lumi;
	  s2L4J+=eventWeight*eventWeight*lumi*lumi;
	}
	else if ( nJet>=2) {
	  c2L2J+=eventWeight*lumi;
	  s2L2J+=eventWeight*eventWeight*lumi*lumi;
	} 
      }
    }
  }

  s2L2J=sqrt(s2L2J);
  s2L4J=sqrt(s2L4J);
  s3L0J=sqrt(s3L0J);
  s3L2J=sqrt(s3L2J);
  s4L0J=sqrt(s4L0J);
  s4L2J=sqrt(s4L2J);
  s5L0J=sqrt(s5L0J);
  s6L0J=sqrt(s6L0J);  

  std::cout << std::scientific;

  std::cout << isamp << "\t $ ";
  std::cout << std::setprecision(1) << c2L2J << " \\pm " << std::setprecision(1)<< s2L2J <<  " $ " << std::setprecision(1)<< c2L4J << " \\pm " << std::setprecision(1)<< s2L4J << " $ ";
  std::cout << std::setprecision(1) << c3L0J << " \\pm " << std::setprecision(1)<< s3L0J <<  " $ " << std::setprecision(1)<< c3L2J << " \\pm " << std::setprecision(1)<< s3L2J << " $ ";
  std::cout << std::setprecision(1) << c4L0J << " \\pm " << std::setprecision(1)<< s4L0J <<  " $ " << std::setprecision(1)<< c4L2J << " \\pm " << std::setprecision(1)<< s4L2J << " $ ";
  std::cout << std::setprecision(1) << c5L0J << " \\pm " << std::setprecision(1)<< s5L0J <<  " $ " << std::setprecision(1)<< c6L0J << " \\pm " << std::setprecision(1)<< s6L0J << " \\\\ " << std::endl;
  //std::cout  << isamp << " $ " << c2L4J << " $ " << c4L2J << " \\\\ " << endl;

  delete inTree; inTree=0;
  delete inFile; inFile=0;
  
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

Int_t findZ( TClonesArray * vLep, int &zll_1, int &zll_2, int &zll_3, int &zll_4, float massRange) {

  zll_1=-1; zll_2=-1; zll_3=-1; zll_4=-1;

  if (vLep->GetEntries()>=2) {

    for (int i=0; i<vLep->GetEntries(); i++) {
      TLepton *lep1 = (TLepton*) (*vLep)[i];
      TLorentzVector v1; v1.SetPtEtaPhiM(lep1->PT, lep1->Eta, lep1->Phi, (lep1->PID == 11) ? ELE_MASS : MUON_MASS);

      for (int j=i+1; j<vLep->GetEntries(); j++) {
	TLepton *lep2 = (TLepton*) (*vLep)[j];
	TLorentzVector v2; v2.SetPtEtaPhiM(lep2->PT, lep2->Eta, lep2->Phi, (lep2->PID == 11) ? ELE_MASS : MUON_MASS);

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

  if (vLep->GetEntries()>=4) {

    for (int i=0; i<vLep->GetEntries(); i++) {
      if (i==zll_1 || i==zll_2) continue;

      TLepton *lep1 = (TLepton*) (*vLep)[i];
      TLorentzVector v1; v1.SetPtEtaPhiM(lep1->PT, lep1->Eta, lep1->Phi, (lep1->PID == 11) ? ELE_MASS : MUON_MASS);

      for (int j=i+1; j<vLep->GetEntries(); j++) {
	if (j==zll_1 || j==zll_2) continue;

	TLepton *lep2 = (TLepton*) (*vLep)[j];
	TLorentzVector v2; v2.SetPtEtaPhiM(lep2->PT, lep2->Eta, lep2->Phi, (lep2->PID == 11) ? ELE_MASS : MUON_MASS);

	if (lep1->Charge == lep2->Charge) continue;

	if ( abs( (v1+v2).M() - Z_MASS ) < massRange ) {
	  zll_3=i; zll_4=j;
	  break;
	}
      }
      if (zll_3 != zll_4) break;
    }
  }

  if (zll_3 == zll_4) return 1;

  else return 2;
}

Int_t findZjets( TClonesArray * vLep, int &zll_1, int &zll_2, float massRange) {

  zll_1=-1; zll_2=-1;

  if (vLep->GetEntries()>=2) {

    for (int i=0; i<vLep->GetEntries(); i++) {
      THadJet *lep1 = (THadJet*) (*vLep)[i];
      TLorentzVector v1; v1.SetPtEtaPhiM(lep1->PT, lep1->Eta, lep1->Phi, lep1->Mass);

      for (int j=i+1; j<vLep->GetEntries(); j++) {
	THadJet *lep2 = (THadJet*) (*vLep)[j];
	TLorentzVector v2; v2.SetPtEtaPhiM(lep2->PT, lep2->Eta, lep2->Phi, lep2->Mass);

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
