//-------------------------------------------------------------------
// Select bbtautau events from delphes
//
// execute with:
// root -l -q selectDelphes.C+\(\"_inputfile_name_\",_file_cross_section_,\"_outputfile_name\"\)
//
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

#include "modules/Delphes.h"                   // delphes
#include "ExRootAnalysis/ExRootTreeReader.h"   // delphes
#include "classes/DelphesClasses.h"            // delphes

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

void selectDelphes(const TString inputfile="/eos/cms/store/group/upgrade/delphes_output/YR_Delphes/Delphes342pre15/ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_200PU/ZZZ_TuneCUETP8M1_14TeV-amcatnlo-pythia8_10_0.root",
		   const Float_t xsec=1341.36923,
		   const TString outputfile="foo.root") {

  // declare constants
  const Double_t MUON_MASS = 0.105658369;
  const Double_t ELE_MASS  = 0.000511;

  const Int_t B_ID_CODE = 5;
  const Int_t TAU_ID_CODE = 15;

  const Float_t MAX_MATCH_DIST = 0.3;

  // tau decay modes
  enum { hadron=1, electron, muon };

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // set up branches to read in from file
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("MuonLoose");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  // set up loop variables
  GenParticle *genParticle;
  Jet *jet, *genJet;
  Electron *ele;
  Muon *mu;
  // comment out following line for di-higgs samples
  LHEFEvent *event;

  // set up storage variables
  GenParticle *tau1, *tau2;
  Jet *genJetTau1, *genJetTau2;
  Jet *jetTau1, *jetTau2;
  Electron *eleTau1, *eleTau2;
  Muon *muTau1, *muTau2;
  GenParticle *b1, *b2;
  Jet *jetB1, *jetB2;

  // set up output variables and file
  Int_t nEvents;
  Float_t eventWeight;
  UInt_t tauDecayCat1, tauDecayCat2;
  UInt_t bTag1, bTag2;
  LorentzVector *genTau1=0, *genTau2=0, *genDecayTau1=0, *genDecayTau2=0, *recoTau1=0, *recoTau2=0;
  LorentzVector *genB1=0, *genB2=0, *recoB1=0, *recoB2=0;

  TFile *outFile = new TFile(outputfile, "RECREATE");

  // tree to hold the number of events in the file before selection
  TTree *sampTree = new TTree("Info", "Info");
  sampTree->Branch("nEvents",       &nEvents,        "nEvents/i");
  nEvents=numberOfEntries;
  sampTree->Fill();

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

  // set up more storage variables for jet ordering
  Int_t iGenTau1=-1,    iGenTau2=-1;
  Int_t iGenDecay1=-1,  iGenDecay2=-1;
  Int_t iRecoDecay1=-1, iRecoDecay2=-1;
  Int_t iGenB1=-1,      iGenB2=-1;
  Int_t iRecoB1=-1,     iRecoB2=-1;

  // define placeholder vector for things that don't exist
  LorentzVector nothing(999,999,0,999);

  // counters for b-quark and tau candidates
  Int_t bCan=0; Int_t tCan=0;

  //for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
  for (Int_t iEntry=0; iEntry<1000; iEntry++) { // entry loop
    treeReader->ReadEntry(iEntry);

    // comment out following line for di-higgs samples
    event = (LHEFEvent*) branchEvent->At(0);
    eventWeight = 1;
    eventWeight *= xsec;
    // comment out following line for di-higgs samples
    eventWeight *= event->Weight;

    // initialize index holders to -1
    iGenTau1=-1;    iGenTau2=-1;
    iGenDecay1=-1;  iGenDecay2=-1;
    iRecoDecay1=-1; iRecoDecay2=-1;
    iGenB1=-1;      iGenB2=-1;
    iRecoB1=-1;     iRecoB2=-1;

    // initialize tags to 0
    tauDecayCat1=0; tauDecayCat2=0;
    bTag1=0;        bTag2=0;
    
    // initialize counters to 0
    bCan=0; tCan=0;

    // get reconstructed hadronic taus and b jets
    for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) { // reconstructed jet loop
      jet = (Jet*) branchJet->At(iJet);

      if ( jet->TauTag !=0 ) { // tau tag switch
	tCan++; // found a tau candidate

	// logic to order candidates by pT
	if (iRecoDecay1 == -1) {                       // first tau candidate
	  iRecoDecay1 = iJet;                          // save the index in branchJet
	  jetTau1 = (Jet*) branchJet->At(iRecoDecay1); // make jet object from index
	  tauDecayCat1 = hadron;                       // set hadronic decay flag
	}
	else if ( jet->PT > jetTau1->PT ) {
	  iRecoDecay2 = iRecoDecay1;
	  jetTau2 = (Jet*) branchJet->At(iRecoDecay2);
	  tauDecayCat2 = tauDecayCat1;
	  iRecoDecay1 = iJet;
	  jetTau1 = (Jet*) branchJet->At(iRecoDecay1);
	  tauDecayCat1 = hadron;
	}
	else if ( iRecoDecay2 == -1) {
	  iRecoDecay2 = iJet;
	  jetTau2 = (Jet*) branchJet->At(iRecoDecay2);
	  tauDecayCat2 = hadron;
	}
	else if ( jet->PT > jetTau2->PT ) {
	  iRecoDecay2 = iJet;
	  jetTau2 = (Jet*) branchJet->At(iRecoDecay2);	  
	  tauDecayCat2 = hadron;
	}
	  	
      } // end tau tag switch

      if ( jet->BTag !=0) { // b tag switch
	bCan++; // found a b-quark candidate

	// logic to order b-jet candidates by pT
	if (iRecoB1 == -1) {
	  iRecoB1 = iJet;
	  jetB1 = (Jet*) branchJet->At(iRecoB1);
	  bTag1 = jetB1->BTag;                        // store btag flag
	}
	else if ( jet->PT > jetB1->PT ) {
	  iRecoB2 = iRecoB1;
	  jetB2 = (Jet*) branchJet->At(iRecoB2);
	  bTag2 = bTag1;

	  iRecoB1 = iJet;
	  jetB1 = (Jet*) branchJet->At(iRecoB1);
	  bTag1 = jetB1->BTag;
	}
	else if ( iRecoB2 == -1) {
	  iRecoB2 = iJet;
	  jetB2 = (Jet*) branchJet->At(iRecoB2);
	  bTag2 = jetB2->BTag;
	}
	else if ( jet->PT > jetB2->PT ) {
	  iRecoB2 = iJet;
	  jetB2 = (Jet*) branchJet->At(iRecoB2);	  
	  bTag2 = jetB2->BTag;
	}

      } // end b tag switch

    } // end reco jet loop

    // get info for taus that decay to muons
    for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) { // reco muon loop
      mu = (Muon*) branchMuon->At(iMuon);
      tCan++; // found a tau candidate

      // further pT ordering logic
      if (iRecoDecay1 == -1) {
	iRecoDecay1 = iMuon;
	muTau1 = (Muon*) branchMuon->At(iRecoDecay1);
	tauDecayCat1 = muon;
      }
      else if ( ( tauDecayCat1 == hadron) && ( mu->PT > jetTau1->PT ) ) {
	iRecoDecay2 = iRecoDecay1;
	jetTau2 = (Jet*) branchJet->At(iRecoDecay2);
	tauDecayCat2 = tauDecayCat1;
	
	iRecoDecay1 = iMuon;
	muTau1 = (Muon*) branchMuon->At(iRecoDecay1);
	tauDecayCat1 = muon;
      }
      else if ( ( tauDecayCat1 == muon) && ( mu->PT > muTau1->PT ) ) {
	iRecoDecay2 = iRecoDecay1;
	muTau2 = (Muon*) branchMuon->At(iRecoDecay2);
	tauDecayCat2 = tauDecayCat1;
	
	iRecoDecay1 = iMuon;
	muTau1 = (Muon*) branchMuon->At(iRecoDecay1);
	tauDecayCat1 = muon;
      }
      else if ( iRecoDecay2 == -1) {
	iRecoDecay2 = iMuon;
	muTau2 = (Muon*) branchMuon->At(iRecoDecay2);
	tauDecayCat2 = muon;
      }
      else if ( ( tauDecayCat2 == hadron ) && ( mu->PT > jetTau2->PT ) ) {
	iRecoDecay2 = iMuon;
	muTau2 = (Muon*) branchMuon->At(iRecoDecay2);	  
	tauDecayCat2 = muon;
      }
      else if ( ( tauDecayCat2 == muon ) && ( mu->PT > muTau2->PT ) ) {
	iRecoDecay2 = iMuon;
	muTau2 = (Muon*) branchMuon->At(iRecoDecay2);	  
	tauDecayCat2 = muon;
      }
      
    } // end muon loop

    // get info for taus that decay to electrons
    for (Int_t iEle=0; iEle<branchElectron->GetEntries(); iEle++) { // reco ele loop
      ele = (Electron*) branchElectron->At(iEle);
      tCan++; // found a tau candidate

      // and more pT odering logic
      if (iRecoDecay1 == -1) {
	iRecoDecay1 = iEle;
	eleTau1 = (Electron*) branchElectron->At(iRecoDecay1);
	tauDecayCat1 = electron;
      }
      else if ( ( tauDecayCat1 == hadron) && ( ele->PT > jetTau1->PT ) ) {
	iRecoDecay2 = iRecoDecay1;
	jetTau2 = (Jet*) branchJet->At(iRecoDecay2);
	tauDecayCat2 = tauDecayCat1;
	
	iRecoDecay1 = iEle;
	eleTau1 = (Electron*) branchElectron->At(iRecoDecay1);
	tauDecayCat1 = electron;
      }
      else if ( ( tauDecayCat1 == muon) && ( ele->PT > muTau1->PT ) ) {
	iRecoDecay2 = iRecoDecay1;
	muTau2 = (Muon*) branchMuon->At(iRecoDecay2);
	tauDecayCat2 = tauDecayCat1;
	
	iRecoDecay1 = iEle;
	eleTau1 = (Electron*) branchElectron->At(iRecoDecay1);
	tauDecayCat1 = electron;
      }
      else if ( ( tauDecayCat1 == electron) && ( ele->PT > eleTau1->PT ) ) {
	iRecoDecay2 = iRecoDecay1;
	eleTau2 = (Electron*) branchElectron->At(iRecoDecay2);
	tauDecayCat2 = tauDecayCat1;
	
	iRecoDecay1 = iEle;
	eleTau1 = (Electron*) branchElectron->At(iRecoDecay1);
	tauDecayCat1 = electron;
      }
      else if ( iRecoDecay2 == -1) {
	iRecoDecay2 = iEle;
	eleTau2 = (Electron*) branchElectron->At(iRecoDecay2);
	tauDecayCat2 = electron;
      }
      else if ( ( tauDecayCat2 == hadron ) && ( ele->PT > jetTau2->PT ) ) {
	iRecoDecay2 = iEle;
	eleTau2 = (Electron*) branchElectron->At(iRecoDecay2);	  
	tauDecayCat2 = electron;
      }
      else if ( ( tauDecayCat2 == muon ) && ( ele->PT > muTau2->PT ) ) {
	iRecoDecay2 = iEle;
	eleTau2 = (Electron*) branchElectron->At(iRecoDecay2);	  
	tauDecayCat2 = electron;
      }
      else if ( ( tauDecayCat2 == electron ) && ( ele->PT > eleTau2->PT ) ) {
	iRecoDecay2 = iEle;
	eleTau2 = (Electron*) branchElectron->At(iRecoDecay2);	  
	tauDecayCat2 = electron;
      }
      
    } // end electron loop

    // skip event if there weren't two reconstructed b-jets and two tau
    if ( ( tCan<2 ) ) continue;

    // fill 4-vector for leading b-jet
    LorentzVector vRecoB1(0,0,0,0);
    if (bTag1==0) recoB1 = &nothing;
    else {
      vRecoB1.SetPt(jetB1->PT);
      vRecoB1.SetEta(jetB1->Eta);
      vRecoB1.SetPhi(jetB1->Phi);
      vRecoB1.SetM(jetB1->Mass);
      recoB1 = &vRecoB1;
    }

    // fill 4-vector for second b-jet
    LorentzVector vRecoB2(0,0,0,0);
    if (bTag2==0) recoB2 = &nothing;
    else {
      vRecoB2.SetPt(jetB2->PT);
      vRecoB2.SetEta(jetB2->Eta);
      vRecoB2.SetPhi(jetB2->Phi);
      vRecoB2.SetM(jetB2->Mass);
      recoB2 = &vRecoB2;
    }

    // fill 4-vector for leading tau
    LorentzVector vRecoTau1(0,0,0,0);
    if ( tauDecayCat1 == 0 ) recoTau1 = &nothing;
    else if ( tauDecayCat1 == hadron ) {
      vRecoTau1.SetPt(jetTau1->PT);
      vRecoTau1.SetEta(jetTau1->Eta);
      vRecoTau1.SetPhi(jetTau1->Phi);
      vRecoTau1.SetM(jetTau1->Mass);
      recoTau1 = &vRecoTau1;
    }
    else if ( tauDecayCat1 == muon ) {
      vRecoTau1.SetPt(muTau1->PT);
      vRecoTau1.SetEta(muTau1->Eta);
      vRecoTau1.SetPhi(muTau1->Phi);
      vRecoTau1.SetM(MUON_MASS);
      recoTau1 = &vRecoTau1;
    }
    else if ( tauDecayCat1 == electron ) {
      vRecoTau1.SetPt(eleTau1->PT);
      vRecoTau1.SetEta(eleTau1->Eta);
      vRecoTau1.SetPhi(eleTau1->Phi);
      vRecoTau1.SetM(ELE_MASS);
      recoTau1 = &vRecoTau1;
    }

    // fill 4-vector for second tau
    LorentzVector vRecoTau2(0,0,0,0);
    if ( tauDecayCat2 == 0 ) recoTau2 = &nothing;
    else if ( tauDecayCat2 == hadron ) {
      vRecoTau2.SetPt(jetTau2->PT);
      vRecoTau2.SetEta(jetTau2->Eta);
      vRecoTau2.SetPhi(jetTau2->Phi);
      vRecoTau2.SetM(jetTau2->Mass);
      recoTau2 = &vRecoTau2;
    }
    else if ( tauDecayCat2 == muon ) {
      vRecoTau2.SetPt(muTau2->PT);
      vRecoTau2.SetEta(muTau2->Eta);
      vRecoTau2.SetPhi(muTau2->Phi);
      vRecoTau2.SetM(MUON_MASS);
      recoTau2 = &vRecoTau2;
    }
    else if ( tauDecayCat2 == electron ) {
      vRecoTau2.SetPt(eleTau2->PT);
      vRecoTau2.SetEta(eleTau2->Eta);
      vRecoTau2.SetPhi(eleTau2->Phi);
      vRecoTau2.SetM(ELE_MASS);
      recoTau2 = &vRecoTau2;
    }

    // now match generator level particles to reconstructed jets/particles
    for (Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) { // generator particle loop
      genParticle = (GenParticle*) branchParticle->At(iParticle);

      if ( fabs(genParticle->PID) == TAU_ID_CODE ) { // tau switch
	if ( deltaR(genParticle->Eta, vRecoTau1.Eta(), genParticle->Phi, vRecoTau1.Phi()) < MAX_MATCH_DIST ) { // distance check
	  iGenTau1 = iParticle;
	  tau1 = (GenParticle*) branchParticle->At(iGenTau1);
	}
	else if ( deltaR(genParticle->Eta, vRecoTau2.Eta(), genParticle->Phi, vRecoTau2.Phi()) < MAX_MATCH_DIST) {
	  iGenTau2 = iParticle;
	  tau2 = (GenParticle*) branchParticle->At(iGenTau2);
	}
      }

      if ( fabs(genParticle->PID) == B_ID_CODE ) { // b-quark switch
	if ( deltaR(genParticle->Eta, vRecoB1.Eta(), genParticle->Phi, vRecoB1.Phi()) < MAX_MATCH_DIST ) {
	  iGenB1 = iParticle;
	  b1 = (GenParticle*) branchParticle->At(iGenB1);
	}
	if ( deltaR(genParticle->Eta, vRecoB2.Eta(), genParticle->Phi, vRecoB2.Phi()) < MAX_MATCH_DIST ) {
	  iGenB2 = iParticle;
	  b2 = (GenParticle*) branchParticle->At(iGenB2);
	}
      }

    } // end particle loop
    
    // match generator level jets to tau candidates (there will be missing pT from 1-2 neutrinos in tau decay
    for (Int_t iJet=0; iJet<branchGenJet->GetEntries(); iJet++) { // generator level jet loop
      genJet = (Jet*) branchGenJet->At(iJet);
  
      if ( deltaR(genJet->Eta, vRecoTau1.Eta(), genJet->Phi, vRecoTau1.Phi()) < MAX_MATCH_DIST ) {
	iGenDecay1 = iJet;
	genJetTau1 = (Jet*) branchGenJet->At(iGenDecay1);
      }
      else if ( deltaR(genJet->Eta, vRecoTau2.Eta(), genJet->Phi, vRecoTau2.Phi()) < MAX_MATCH_DIST) {
	iGenDecay2 = iJet;
	genJetTau2 = (Jet*) branchGenJet->At(iGenDecay2);
      }
    }
    
    // 
    // OUTPUT 
    //

    // store generator particles
    LorentzVector vGenTau1(0,0,0,0);
    if ( iGenTau1 == -1) genTau1 = &nothing;
    else {
      vGenTau1.SetPt(tau1->PT);
      vGenTau1.SetEta(tau1->Eta);
      vGenTau1.SetPhi(tau1->Phi);
      vGenTau1.SetM(tau1->Mass);
      genTau1 = &vGenTau1;
    }

    LorentzVector vGenTau2(0,0,0,0);
    if ( iGenTau2 == -1) genTau2 = &nothing;
    else {
      vGenTau2.SetPt(tau2->PT);
      vGenTau2.SetEta(tau2->Eta);
      vGenTau2.SetPhi(tau2->Phi);
      vGenTau2.SetM(tau2->Mass);
      genTau2 = &vGenTau2;
    }

    LorentzVector vGenB1(0,0,0,0);
    if ( iGenB1 == -1) genB1 = &nothing;
    else {
      vGenB1.SetPt(b1->PT);
      vGenB1.SetEta(b1->Eta);
      vGenB1.SetPhi(b1->Phi);
      vGenB1.SetM(b1->Mass);
      genB1 = &vGenB1;
    }

    LorentzVector vGenB2(0,0,0,0);
    if ( iGenB2 == -1) genB2 = &nothing;
    else {
      vGenB2.SetPt(b2->PT);
      vGenB2.SetEta(b2->Eta);
      vGenB2.SetPhi(b2->Phi);
      vGenB2.SetM(b2->Mass);
      genB2 = &vGenB2;
    }

    // store generator jets
    LorentzVector vGenJetTau1(0,0,0,0);
    if ( iGenDecay1 == -1 ) genDecayTau1 = &nothing;
    else {
      vGenJetTau1.SetPt(genJetTau1->PT);
      vGenJetTau1.SetEta(genJetTau1->Eta);
      vGenJetTau1.SetPhi(genJetTau1->Phi);
      vGenJetTau1.SetM(genJetTau1->Mass);
      genDecayTau1 = &vGenJetTau1;
    }
    LorentzVector vGenJetTau2(0,0,0,0);
    if ( iGenDecay2 == -1 ) genDecayTau2 = &nothing;
    else {
      vGenJetTau2.SetPt(genJetTau2->PT);
      vGenJetTau2.SetEta(genJetTau2->Eta);
      vGenJetTau2.SetPhi(genJetTau2->Phi);
      vGenJetTau2.SetM(genJetTau2->Mass);
      genDecayTau2 = &vGenJetTau2;
    }

    // save event to tree
    outTree->Fill();

  } // end event loop

  // save file
  outFile->Write();
  // close file
  outFile->Close();

  cout << "----SUMMARY----" << endl;
  cout << " input file " << inputfile << " selection done " << endl;

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
