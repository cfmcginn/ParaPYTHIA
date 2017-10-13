#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "include/getLogBins.h"
#include "include/kirchnerPalette.h"

int plotQGFraction(const std::string inFileName)
{
  const int nRParam = 5;
  const double rParams[nRParam] = {0.2, 0.3, 0.4, 0.5, 0.6};

  kirchnerPalette colors;

  const Int_t nMaxJets = 10000;
  Int_t ngenjt_[nRParam];
  Float_t genjtpt_[nRParam][nMaxJets];
  Int_t genjtpart_[nRParam][nMaxJets];

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p[nRParam];

  for(int i = 0; i < nRParam; ++i){
    std::string treeName = "ak" + std::to_string(int(rParams[i]*10)) + "GenJetTree";
    inTree_p[i] = (TTree*)inFile_p->Get(treeName.c_str());

    inTree_p[i]->SetBranchStatus("ngen", 1);
    inTree_p[i]->SetBranchStatus("genpt", 1);
    inTree_p[i]->SetBranchStatus("genpart", 1);

    inTree_p[i]->SetBranchAddress("ngen", &ngenjt_[i]);
    inTree_p[i]->SetBranchAddress("genpt", genjtpt_[i]);
    inTree_p[i]->SetBranchAddress("genpart", genjtpart_[i]);
  }

  const Int_t nJtPtBins = 20;
  const Float_t jtPtBinsLow = 80;
  const Float_t jtPtBinsHi = 350;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtBinsLow, jtPtBinsHi, nJtPtBins, jtPtBins);

  TCanvas* qgFrac_p[nRParam];
  TH1F* incSpect_h[nRParam];
  TH1F* qSpect_h[nRParam];
  TH1F* gSpect_h[nRParam];
  TH1F* uSpect_h[nRParam];
  for(int i = 0; i < nRParam; ++i){
    const std::string canvName = "qgFracR" + std::to_string(int(rParams[i]*10)) + "_c";
    const std::string incName = "incSpectR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string qName = "qSpectR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string gName = "gSpectR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string uName = "uSpectR" + std::to_string(int(rParams[i]*10)) + "_h";

    qgFrac_p[i] = new TCanvas(canvName.c_str(), canvName.c_str(), 500, 500);
    incSpect_h[i] = new TH1F(incName.c_str(), ";Jet p_{T} (Inc.);Events", nJtPtBins, jtPtBins);
    qSpect_h[i] = new TH1F(qName.c_str(), ";Jet p_{T} (Quark);Events", nJtPtBins, jtPtBins);
    gSpect_h[i] = new TH1F(gName.c_str(), ";Jet p_{T} (Gluon);Events", nJtPtBins, jtPtBins);
    uSpect_h[i] = new TH1F(uName.c_str(), ";Jet p_{T} (Untagged);Events", nJtPtBins, jtPtBins);

    incSpect_h[i]->Sumw2();
    qSpect_h[i]->Sumw2();
    gSpect_h[i]->Sumw2();
    uSpect_h[i]->Sumw2();

    incSpect_h[i]->GetXaxis()->CenterTitle();
    qSpect_h[i]->GetXaxis()->CenterTitle();
    gSpect_h[i]->GetXaxis()->CenterTitle();
    uSpect_h[i]->GetXaxis()->CenterTitle();

    incSpect_h[i]->GetYaxis()->CenterTitle();
    qSpect_h[i]->GetYaxis()->CenterTitle();
    gSpect_h[i]->GetYaxis()->CenterTitle();
    uSpect_h[i]->GetYaxis()->CenterTitle();
  }

  const Int_t nEntries = inTree_p[0]->GetEntries();

  for(Int_t entry = 0; entry < nEntries; ++entry){
    for(int i = 0; i < nRParam; ++i){
      inTree_p[i]->GetEntry(entry);

      for(Int_t j = 0; j < ngenjt_[i]; ++j){
	if(genjtpt_[i][j] > jtPtBinsHi) continue;
	if(genjtpt_[i][j] < jtPtBinsLow) continue;

	incSpect_h[i]->Fill(genjtpt_[i][j]);
	if(TMath::Abs(genjtpart_[i][j]) < 9) qSpect_h[i]->Fill(genjtpt_[i][j]);
	else if(TMath::Abs(genjtpart_[i][j]) == 21) gSpect_h[i]->Fill(genjtpt_[i][j]);
	else uSpect_h[i]->Fill(genjtpt_[i][j]);
      }
    }
  }

  

  for(int i = 0; i < nRParam; ++i){
    qSpect_h[i]->Divide(incSpect_h[i]);
    gSpect_h[i]->Divide(incSpect_h[i]);
    uSpect_h[i]->Divide(incSpect_h[i]);

    qgFrac_p[i]->cd();
    gStyle->SetOptStat(0);
    qgFrac_p[i]->SetTopMargin(0.01);
    qgFrac_p[i]->SetRightMargin(0.01);
    qSpect_h[i]->SetFillColor(kBlue);
    gSpect_h[i]->SetFillColor(kRed);
    uSpect_h[i]->SetFillColor(kYellow+2);

    gSpect_h[i]->SetMinimum(0.);
    gSpect_h[i]->SetMaximum(1.1);

    gSpect_h[i]->GetXaxis()->SetTitle("Gen. Jet p_{T} (GeV/c)");
    gSpect_h[i]->GetYaxis()->SetTitle(("Fraction (#color[2]{G},#color[4]{Q},#color[" + std::to_string(kYellow+2) + "]{Untagged})").c_str());
    gSpect_h[i]->DrawCopy("HIST");
    qSpect_h[i]->Add(gSpect_h[i]);
    qSpect_h[i]->DrawCopy("HIST SAME");

    uSpect_h[i]->Add(qSpect_h[i]);
    uSpect_h[i]->DrawCopy("HIST SAME");

    qSpect_h[i]->DrawCopy("HIST SAME");
    gSpect_h[i]->DrawCopy("HIST SAME");

    gPad->SetLogx();
    gPad->RedrawAxis();

    std::string saveName = "pdfDir/qgFracR" + std::to_string(int(rParams[i]*10)) + ".pdf";
    qgFrac_p[i]->SaveAs(saveName.c_str());

    delete qgFrac_p[i];
    delete incSpect_h[i];
    delete qSpect_h[i];
    delete gSpect_h[i];
    delete uSpect_h[i];
  }  

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./plotQGFraction.exe <inFileName>" << std::endl;
    return 1;
  }
 
  int retVal = 0;
  if(argc == 2) retVal += plotQGFraction(argv[1]);
  return retVal;
}
