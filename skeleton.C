#if !defined(__CINT__) || defined(__MAKECINT__)
// include statements for all needed dependencies
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>

// include statements for Delphes
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#endif

void skeleton(const TString inputfile="HHToGGBB_14TeV_0.root")
{

  // read input input file
  TChain chain("Delphes");
  chain.Add(inputfile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // set up loop variable
  Photon *photon;

  // set up storage variables
  Int_t iMaxPtPhoton=-1;
  Photon *maxPtPhoton;

  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  // set up output histogram
  TH1D *hPhotPt = new TH1D("hPhotPt", "hPhotPt", 30, 0, 300);

  for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // event loop
    treeReader->ReadEntry(iEntry);

    // reset counter variable
    iMaxPtPhoton=-1;

    // skip event if it diesn't have at least one photon
    if (branchPhoton->GetEntries()<1) continue;

    for (Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) { // photon loop
      photon = (Photon*) branchPhoton->At(iPhoton);

      // for first photon in event
      if (iMaxPtPhoton==-1) {
	// assume it has the highest pT, store its # from list 
	iMaxPtPhoton=iPhoton;
	// and its object
	maxPtPhoton= (Photon*) branchPhoton->At(iMaxPtPhoton);
      }
      // for all other photons
      else {
	// check if this photon's pT is greater than the current max photon pT
	if (photon->PT > maxPtPhoton->PT ) {
	  // if so, store its # from list
	  iMaxPtPhoton=iPhoton;
	  // and its object
	  maxPtPhoton= (Photon*) branchPhoton->At(iMaxPtPhoton);
	}
      }

    } // end photon loop

    // add pT of photon with highest pT  to histogram
    hPhotPt->Fill(maxPtPhoton->PT);

  } // end event loop

  // initialize a canvas to draw on
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);

  // draw histogram on canvas
  hPhotPt->Draw();

  // add axis labels
  hPhotPt->GetXaxis()->SetTitle("Leading Photon pT");
  hPhotPt->GetYaxis()->SetTitle("Count");
  hPhotPt->SetTitle(""); // title on top

}
