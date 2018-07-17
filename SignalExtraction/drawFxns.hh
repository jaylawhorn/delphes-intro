#ifndef DRAW_FXNS_HH
#define DRAW_FXNS_HH

struct HistHolder {
  TH1D *zee; 
  TH1D *zmm; 

  TH1D *zeeL; 
  TH1D *zeeS;
  TH1D *zmmL; 
  TH1D *zmmS;

  TH1D *nE; 
  TH1D *nM; 
  TH1D *nJ;
};

void initHistHolder(HistHolder &h, TString id) {

  const int nBins_zMass=34, xMin_zMass=75, xMax_zMass=107;
  const int nBins_lPT  =25, xMin_lPT  = 0, xMax_lPT  =150;
  const int nBins_nObj =10, xMin_nObj = 0, xMax_nObj = 10;

  h.zee = new TH1D(Form("%s_zee",id.Data()), Form("%s_zee",id.Data()), 
		   nBins_zMass, xMin_zMass, xMax_zMass);
  h.zee->GetXaxis()->SetTitle("m_{ee} [GeV]");
  h.zee->GetYaxis()->SetTitle("Events");

  h.zmm = new TH1D(Form("%s_zmm",id.Data()), Form("%s_zmm",id.Data()), 
		   nBins_zMass, xMin_zMass, xMax_zMass);

  h.zmm->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  h.zmm->GetYaxis()->SetTitle("Events");

  h.zeeL = new TH1D(Form("%s_zeeL",id.Data()), Form("%s_zeeL",id.Data()),
		    nBins_lPT, xMin_lPT, xMax_lPT);

  h.zeeL->GetXaxis()->SetTitle("lead e p_{T} (Z#arrow ee) [GeV]");
  h.zeeL->GetYaxis()->SetTitle("Events");

  h.zeeS = new TH1D(Form("%s_zeeS",id.Data()), Form("%s_zeeS",id.Data()),
		    nBins_lPT, xMin_lPT, xMax_lPT);

  h.zeeS->GetXaxis()->SetTitle("2nd e p_{T} (Z#arrow ee) [GeV]");
  h.zeeS->GetYaxis()->SetTitle("Events");

  h.zmmL = new TH1D(Form("%s_zmmL",id.Data()), Form("%s_zmmL",id.Data()),
		    nBins_lPT, xMin_lPT, xMax_lPT);

  h.zmmL->GetXaxis()->SetTitle("lead #mu p_{T} (Z#arrow#mu#mu) [GeV]");
  h.zmmL->GetYaxis()->SetTitle("Events");

  h.zmmS = new TH1D(Form("%s_zmmS",id.Data()), Form("%s_zmmS",id.Data()),
		    nBins_lPT, xMin_lPT, xMax_lPT);

  h.zmmS->GetXaxis()->SetTitle("2nd #mu p_{T} (Z#arrow#mu#mu) [GeV]");
  h.zmmS->GetYaxis()->SetTitle("Events");

  h.nE = new TH1D(Form("%s_nE",id.Data()), Form("%s_nE",id.Data()), 
		  nBins_nObj, xMin_nObj, xMax_nObj);

  h.nE->GetXaxis()->SetTitle("n_{e} (p_{T}>10 GeV)");
  h.nE->GetYaxis()->SetTitle("Events");

  h.nM = new TH1D(Form("%s_nM",id.Data()), Form("%s_nM",id.Data()), 
		  nBins_nObj, xMin_nObj, xMax_nObj);

  h.nM->GetXaxis()->SetTitle("n_{#mu} (p_{T}>10 GeV)");
  h.nM->GetYaxis()->SetTitle("Events");

  h.nJ = new TH1D(Form("%s_nJ",id.Data()), Form("%s_nJ",id.Data()), 
		  nBins_nObj, xMin_nObj, xMax_nObj);

  h.nJ->GetXaxis()->SetTitle("n_{j} (p_{T}>15 GeV)");
  h.nJ->GetYaxis()->SetTitle("Events");

}

void DrawHists(TCanvas *c, TString name, TH1D* signal, TH1D* ttbar, TH1D* otherb) {

  signal->SetLineColor(kBlue); signal->SetFillColor(kBlue); signal->SetFillStyle(1001);
  ttbar->SetLineColor(kRed); ttbar->SetFillColor(kRed); ttbar->SetFillStyle(1001);
  otherb->SetLineColor(kGray); otherb->SetFillColor(kGray); otherb->SetFillStyle(1001);

  otherb->Add(ttbar);
  signal->Add(otherb);

  signal->SetTitle("");

  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);

  leg->AddEntry(signal, "WWZ", "lf");
  leg->AddEntry(ttbar, "TT(Z)", "lf");
  leg->AddEntry(otherb, "Other", "lf");

  signal->Draw("hist");
  otherb->Draw("hist same");
  ttbar->Draw("hist same");

  leg->Draw();

  c->SaveAs(name);
  
}

#endif
