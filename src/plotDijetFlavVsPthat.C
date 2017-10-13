#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TLine.h"

#include "include/getLogBins.h"
#include "include/getLinBins.h"
#include "include/etaPhiFunc.h"
#include "include/kirchnerPalette.h"
#include "include/gausJetEnergyLoss.h"
#include "include/expJetEnergyLoss.h"
#include "include/poissonJetEnergyLoss.h"

void addWithoutErr(TH1F* inHist_p, TH1F* addPlot_p, Double_t multFact = 1)
{
  for(Int_t bIter = 0; bIter < inHist_p->GetNbinsX(); ++bIter){
    inHist_p->SetBinContent(bIter+1, inHist_p->GetBinContent(bIter+1) + addPlot_p->GetBinContent(bIter+1)*multFact);
  }
  return;
}


std::string getSimpleString(Double_t inNum)
{
  std::string retStr = std::to_string(inNum);
  if(retStr.find(".") != std::string::npos){
    retStr.replace(retStr.find("."), 1, "p");
    while(retStr.substr(retStr.size()-1, 1).find("0") != std::string::npos){retStr.replace(retStr.size()-1,1,"");}
  }

  return retStr;
}

std::string getSimpleString2(Double_t inNum)
{
  std::string retStr = std::to_string(inNum);
  if(retStr.find(".") != std::string::npos){
    while(retStr.substr(retStr.size()-1, 1).find("0") != std::string::npos){retStr.replace(retStr.size()-1,1,"");}
  }

  return retStr;
}


void drawSyst(TCanvas* canv_p, TH1F* hist_p, Double_t ySys[], bool doDelta, Int_t col)
{
  TLine* line_p = new TLine();
  line_p->SetLineWidth(2);
  line_p->SetLineColor(col);
  canv_p->cd();

  for(Int_t i = 0; i < hist_p->GetNbinsX(); ++i){
    Double_t low = hist_p->GetBinLowEdge(i+1);
    Double_t w = hist_p->GetBinWidth(i+1);

    Double_t err = ySys[i];
    if(doDelta) err = TMath::Abs(err-hist_p->GetBinContent(i+1));
    Double_t errLow = TMath::Max(0., hist_p->GetBinContent(i+1)-err);
    Double_t errHi = hist_p->GetBinContent(i+1)+err;
    Double_t delErr = errHi-errLow;
    
    line_p->DrawLine(low + w/10., errHi, low+w - w/10., errHi);
    line_p->DrawLine(low + w/10., errHi, low+w/10., errHi-delErr/10.);
    line_p->DrawLine(low +w - w/10., errHi, low+w-w/10., errHi-delErr/10.);


    line_p->DrawLine(low + w/10., errLow, low+w -w/10., errLow);
    line_p->DrawLine(low + w/10., errLow, low+w/10., errLow+delErr/10.);
    line_p->DrawLine(low +w - w/10., errLow, low+w-w/10., errLow+delErr/10.);
  }

  delete line_p;
  return;
}

int plotDijetFlavVsPthat(const std::string inFileName, const std::string pdfName, std::vector<Double_t> dParam, std::vector<Int_t> collParam)
{
  TDatime* date = new TDatime();
  std::string histAppendStr = std::to_string(date->GetDate());
  std::string labelStr;
  if(pdfName.find("Gaus") != std::string::npos && pdfName.size() == 4){
    histAppendStr = "_Gaus_QM" + getSimpleString(dParam.at(0)) + "QS" + getSimpleString(dParam.at(1)) + "_GM" + getSimpleString(dParam.at(2)) + "GS" + getSimpleString(dParam.at(3)) + "_" + histAppendStr;
    labelStr = "Gaus. (#mu_{Q}=" + getSimpleString2(dParam.at(0)) + ", #sigma_{Q}=" + getSimpleString2(dParam.at(1)) + "; #mu_{G}=" + getSimpleString2(dParam.at(2)) +", #sigma_{G}=" + getSimpleString2(dParam.at(3)) + ")";
  }
  else if(pdfName.find("Exp") != std::string::npos && pdfName.size() == 3){
    histAppendStr = "_Exp_QM" + getSimpleString(dParam.at(0)) + "_GM" + getSimpleString(dParam.at(2)) + "_" + histAppendStr;
    labelStr = "Exp. (#tau_{Q}=" + getSimpleString2(dParam.at(0)) + "; #tau_{G}=" + getSimpleString2(dParam.at(2)) + ")";
  }
  else if(pdfName.find("Pois") != std::string::npos && pdfName.size() == 4){
    histAppendStr = "_Pois_QColl" + std::to_string(collParam.at(0)) + "_GColl" + std::to_string(collParam.at(1)) + "_" + "EPer" + getSimpleString(dParam.at(0)) + "_" + histAppendStr;
    labelStr = "Pois. (#lambda_{Q}=" + std::to_string(collParam.at(0)) + "; #lambda_{G}=" + std::to_string(collParam.at(1)) + "; E/Coll=" + getSimpleString2(dParam.at(0)) + ")";
  }
  else{
    histAppendStr = "_NoLoss_" + histAppendStr;
    labelStr = "No Energy Loss";
  }

  delete date;

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

  kirchnerPalette colors;

  const int nRParam = 1;
  //  const double rParams[nRParam] = {0.2, 0.3, 0.4, 0.5, 0.6};
  const double rParams[nRParam] = {0.4};
  const double xjMax[nRParam] = {4.0};

  const double fracMax = 1.2;
  Double_t fracYBins[20+1];
  getLinBins(0.0, fracMax, 20, fracYBins);

  gausJetEnergyLoss lGaus(dParam.at(0), dParam.at(1), dParam.at(2), dParam.at(3));
  expJetEnergyLoss lExp(dParam.at(0), dParam.at(2));
  poissonJetEnergyLoss lPois(collParam.at(0), collParam.at(1), dParam.at(0));
  
  const Int_t nMaxJets = 10000;
  Float_t pthat_[nRParam];
  Int_t ngenjt_[nRParam];
  Float_t genjtpt_[nRParam][nMaxJets];
  Float_t genjteta_[nRParam][nMaxJets];
  Float_t genjtphi_[nRParam][nMaxJets];
  Int_t genjtpart_[nRParam][nMaxJets];
  //  Int_t genjtpartF_[nRParam][nMaxJets];

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p[nRParam];

  for(int i = 0; i < nRParam; ++i){
    std::string treeName = "ak" + std::to_string(int(rParams[i]*10)) + "GenJetTree";
    inTree_p[i] = (TTree*)inFile_p->Get(treeName.c_str());

    inTree_p[i]->SetBranchStatus("pthat", 1);
    inTree_p[i]->SetBranchStatus("ngen", 1);
    inTree_p[i]->SetBranchStatus("genpt", 1);
    inTree_p[i]->SetBranchStatus("geneta", 1);
    inTree_p[i]->SetBranchStatus("genphi", 1);
    inTree_p[i]->SetBranchStatus("genpart", 1);
    //    inTree_p[i]->SetBranchStatus("genpartF", 1);

    inTree_p[i]->SetBranchAddress("pthat", &pthat_[i]);
    inTree_p[i]->SetBranchAddress("ngen", &ngenjt_[i]);
    inTree_p[i]->SetBranchAddress("genpt", genjtpt_[i]);
    inTree_p[i]->SetBranchAddress("geneta", genjteta_[i]);
    inTree_p[i]->SetBranchAddress("genphi", genjtphi_[i]);
    inTree_p[i]->SetBranchAddress("genpart", genjtpart_[i]);
    //    inTree_p[i]->SetBranchAddress("genpartF", genjtpartF_[i]);
  }

  const Int_t nJtPtBins = 10;
  const Float_t jtPtBinsLow = 50;
  const Float_t jtPtBinsHi = 350;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtBinsLow, jtPtBinsHi, nJtPtBins, jtPtBins);

  const Int_t nAJBins = 20;
  const Float_t ajLow = 0.;
  const Float_t ajHi = 1.;
  Double_t ajBins[nAJBins+1];
  getLinBins(ajLow, ajHi, nAJBins, ajBins);

  const Int_t nXJBins = 10;
  //  const Float_t xjLow = 0.3;
  //  const Float_t xjHi = 1.;
  Double_t xjBins[nXJBins+1] = {.3, 0.35647008, 0.39893624, 0.4471955, 0.5002996, 0.56299543, 0.63131267, 0.7078205, 0.79498863, 0.8919211, 1.};
  Double_t xjYVals[nXJBins+1] = {0.05951643, 0.28270304, 0.80099195, 1.7582145, 2.5096095, 2.3409796, 1.7755735, 1.4283943, 1.3763174, 1.3465592};
  Double_t xjYValErr[nXJBins+1] = {0.119180635, 0.48665425, 1.1446307, 2.2122905, 2.9174426, 2.6021104, 1.9540658, 1.5716946, 1.5270019, 1.5046555};
  //  getLinBins(xjLow, xjHi, nXJBins, xjBins);

  Double_t xjYVals_pp[nXJBins+1] = {0.15736067, 0.33045575, 0.57538545, 0.8599136, 1.1220732, 1.3123116, 1.4875939, 1.7272302, 2.1749492, 2.5110443};
  Double_t xjYValErr_pp[nXJBins+1] = {0.20206584, 0.38759342, 0.6324549, 0.92693263, 1.1816528, 1.3620235, 1.5372466, 1.7742089, 2.216944, 2.5406566};



  TH1F* xjATLASR0p4_h = new TH1F("xjATLASR0p4_h", ";;", nXJBins, xjBins);
  TH1F* xjATLASR0p4_pp_h = new TH1F("xjATLASR0p4_pp_h", ";;", nXJBins, xjBins);
  for(Int_t bIter = 0; bIter < xjATLASR0p4_h->GetNbinsX(); ++bIter){
    xjATLASR0p4_h->SetBinContent(bIter+1, xjYVals[bIter]);
    xjATLASR0p4_pp_h->SetBinContent(bIter+1, xjYVals_pp[bIter]);
  }

  xjATLASR0p4_h->SetMarkerStyle(34);
  xjATLASR0p4_h->SetMarkerSize(1.);
  xjATLASR0p4_h->SetMarkerColor(colors.getColor(5));
  xjATLASR0p4_h->SetLineColor(1);

  xjATLASR0p4_pp_h->SetMarkerStyle(34);
  xjATLASR0p4_pp_h->SetMarkerSize(1.);
  xjATLASR0p4_pp_h->SetMarkerColor(colors.getColor(6));
  xjATLASR0p4_pp_h->SetLineColor(1);

  TCanvas* dijetFlav_p[nRParam];
  TCanvas* pthatPastCut_p[nRParam];
  TCanvas* ajFlav_p[nRParam];
  TCanvas* xjFlav_p[nRParam];

  TCanvas* ajFlavFrac_p[nRParam];
  TCanvas* xjFlavFrac_p[nRParam];

  TH1F* incAJ_h[nRParam];
  TH1F* qqAJ_h[nRParam];
  TH1F* ggAJ_h[nRParam];
  TH1F* qgAJ_h[nRParam];
  TH1F* gqAJ_h[nRParam];
  TH1F* uAJ_h[nRParam];

  TH1F* incXJ_h[nRParam];
  TH1F* qqXJ_h[nRParam];
  TH1F* ggXJ_h[nRParam];
  TH1F* qgXJ_h[nRParam];
  TH1F* gqXJ_h[nRParam];
  TH1F* uXJ_h[nRParam];

  TH1F* incPthat_h[nRParam];
  TH1F* qqPthat_h[nRParam];
  TH1F* ggPthat_h[nRParam];
  TH1F* qgPthat_h[nRParam];
  TH1F* gqPthat_h[nRParam];
  TH1F* uPthat_h[nRParam];

  for(int i = 0; i < nRParam; ++i){
    const std::string canvName = "dijetFlavR" + std::to_string(int(rParams[i]*10)) + "_c";
    const std::string canvName2 = "pthatPastCutR" + std::to_string(int(rParams[i]*10)) + "_c";

    const std::string ajCanvName = "ajFlavR" + std::to_string(int(rParams[i]*10)) + "_c";
    const std::string xjCanvName = "xjFlavR" + std::to_string(int(rParams[i]*10)) + "_c";

    const std::string ajFracCanvName = "ajFlavFracR" + std::to_string(int(rParams[i]*10)) + "_c";
    const std::string xjFracCanvName = "xjFlavFracR" + std::to_string(int(rParams[i]*10)) + "_c";

    const std::string incName = "incPthatR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string qqName = "qqPthatR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string ggName = "ggPthatR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string qgName = "qgPthatR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string gqName = "gqPthatR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string uName = "uPthatR" + std::to_string(int(rParams[i]*10)) + "_h";

    const std::string incAJName = "incAJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string qqAJName = "qqAJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string ggAJName = "ggAJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string qgAJName = "qgAJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string gqAJName = "gqAJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string uAJName = "uAJR" + std::to_string(int(rParams[i]*10)) + "_h";

    const std::string incXJName = "incXJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string qqXJName = "qqXJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string ggXJName = "ggXJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string qgXJName = "qgXJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string gqXJName = "gqXJR" + std::to_string(int(rParams[i]*10)) + "_h";
    const std::string uXJName = "uXJR" + std::to_string(int(rParams[i]*10)) + "_h";

    dijetFlav_p[i] = new TCanvas(canvName.c_str(), canvName.c_str(), 500, 500);
    pthatPastCut_p[i] = new TCanvas(canvName2.c_str(), canvName2.c_str(), 500, 500);
    incPthat_h[i] = new TH1F(incName.c_str(), ";pthat (Inc.);Events", nJtPtBins, jtPtBins);
    qqPthat_h[i] = new TH1F(qqName.c_str(), ";pthat (QQ);Events", nJtPtBins, jtPtBins);
    ggPthat_h[i] = new TH1F(ggName.c_str(), ";pthat (GG);Events", nJtPtBins, jtPtBins);
    qgPthat_h[i] = new TH1F(qgName.c_str(), ";pthat (QG);Events", nJtPtBins, jtPtBins);
    gqPthat_h[i] = new TH1F(gqName.c_str(), ";pthat (GQ);Events", nJtPtBins, jtPtBins);
    uPthat_h[i] = new TH1F(uName.c_str(), ";pthat (U);Events", nJtPtBins, jtPtBins);

    ajFlav_p[i] = new TCanvas(ajCanvName.c_str(), ajCanvName.c_str(), 500, 500);
    ajFlavFrac_p[i] = new TCanvas(ajFracCanvName.c_str(), ajCanvName.c_str(), 500, 500);
    incAJ_h[i] = new TH1F(incAJName.c_str(), ";A_{J} (Inc.);Events", nAJBins, ajBins);
    qqAJ_h[i] = new TH1F(qqAJName.c_str(), ";A_{J} (QQ);Events", nAJBins, ajBins);
    ggAJ_h[i] = new TH1F(ggAJName.c_str(), ";A_{J} (GG);Events", nAJBins, ajBins);
    qgAJ_h[i] = new TH1F(qgAJName.c_str(), ";A_{J} (QG);Events", nAJBins, ajBins);
    gqAJ_h[i] = new TH1F(gqAJName.c_str(), ";A_{J} (GQ);Events", nAJBins, ajBins);
    uAJ_h[i] = new TH1F(uAJName.c_str(), ";A_{J} (U);Events", nAJBins, ajBins);

    xjFlav_p[i] = new TCanvas(xjCanvName.c_str(), xjCanvName.c_str(), 500, 500);
    xjFlavFrac_p[i] = new TCanvas(xjFracCanvName.c_str(), xjCanvName.c_str(), 500, 500);
    incXJ_h[i] = new TH1F(incXJName.c_str(), ";x_{J} (Inc.);Events", nXJBins, xjBins);
    qqXJ_h[i] = new TH1F(qqXJName.c_str(), ";x_{J} (QQ);Events", nXJBins, xjBins);
    ggXJ_h[i] = new TH1F(ggXJName.c_str(), ";x_{J} (GG);Events", nXJBins, xjBins);
    qgXJ_h[i] = new TH1F(qgXJName.c_str(), ";x_{J} (QG);Events", nXJBins, xjBins);
    gqXJ_h[i] = new TH1F(gqXJName.c_str(), ";x_{J} (GQ);Events", nXJBins, xjBins);
    uXJ_h[i] = new TH1F(uXJName.c_str(), ";x_{J} (U);Events", nXJBins, xjBins);

    incPthat_h[i]->Sumw2();
    qqPthat_h[i]->Sumw2();
    ggPthat_h[i]->Sumw2();
    qgPthat_h[i]->Sumw2();
    gqPthat_h[i]->Sumw2();
    uPthat_h[i]->Sumw2();

    incPthat_h[i]->GetXaxis()->CenterTitle();
    qqPthat_h[i]->GetXaxis()->CenterTitle();
    ggPthat_h[i]->GetXaxis()->CenterTitle();
    qgPthat_h[i]->GetXaxis()->CenterTitle();
    gqPthat_h[i]->GetXaxis()->CenterTitle();
    uPthat_h[i]->GetXaxis()->CenterTitle();

    incPthat_h[i]->GetYaxis()->CenterTitle();
    qqPthat_h[i]->GetYaxis()->CenterTitle();
    ggPthat_h[i]->GetYaxis()->CenterTitle();
    qgPthat_h[i]->GetYaxis()->CenterTitle();
    gqPthat_h[i]->GetYaxis()->CenterTitle();
    uPthat_h[i]->GetYaxis()->CenterTitle();

    incAJ_h[i]->Sumw2();
    qqAJ_h[i]->Sumw2();
    ggAJ_h[i]->Sumw2();
    qgAJ_h[i]->Sumw2();
    gqAJ_h[i]->Sumw2();
    uAJ_h[i]->Sumw2();

    incAJ_h[i]->GetXaxis()->CenterTitle();
    qqAJ_h[i]->GetXaxis()->CenterTitle();
    ggAJ_h[i]->GetXaxis()->CenterTitle();
    qgAJ_h[i]->GetXaxis()->CenterTitle();
    gqAJ_h[i]->GetXaxis()->CenterTitle();
    uAJ_h[i]->GetXaxis()->CenterTitle();

    incAJ_h[i]->GetYaxis()->CenterTitle();
    qqAJ_h[i]->GetYaxis()->CenterTitle();
    ggAJ_h[i]->GetYaxis()->CenterTitle();
    qgAJ_h[i]->GetYaxis()->CenterTitle();
    gqAJ_h[i]->GetYaxis()->CenterTitle();
    uAJ_h[i]->GetYaxis()->CenterTitle();

    incXJ_h[i]->Sumw2();
    qqXJ_h[i]->Sumw2();
    ggXJ_h[i]->Sumw2();
    qgXJ_h[i]->Sumw2();
    gqXJ_h[i]->Sumw2();
    uXJ_h[i]->Sumw2();

    incXJ_h[i]->GetXaxis()->CenterTitle();
    qqXJ_h[i]->GetXaxis()->CenterTitle();
    ggXJ_h[i]->GetXaxis()->CenterTitle();
    qgXJ_h[i]->GetXaxis()->CenterTitle();
    gqXJ_h[i]->GetXaxis()->CenterTitle();
    uXJ_h[i]->GetXaxis()->CenterTitle();

    incXJ_h[i]->GetYaxis()->CenterTitle();
    qqXJ_h[i]->GetYaxis()->CenterTitle();
    ggXJ_h[i]->GetYaxis()->CenterTitle();
    qgXJ_h[i]->GetYaxis()->CenterTitle();
    gqXJ_h[i]->GetYaxis()->CenterTitle();
    uXJ_h[i]->GetYaxis()->CenterTitle();
  }

  const Int_t nEntries = inTree_p[0]->GetEntries();

  for(Int_t entry = 0; entry < nEntries; ++entry){
    for(int i = 0; i < nRParam; ++i){
      inTree_p[i]->GetEntry(entry);

      Float_t leadingJtPt_ = -999.;
      Float_t subleadingJtPt_ = -999.;
      Float_t leadingJtPhi_ = -999.;
      Float_t subleadingJtPhi_ = -999.;
      Int_t leadingJtFlav_ = -999;
      Int_t subleadingJtFlav_ = -999;

      for(Int_t j = 0; j < ngenjt_[i]; ++j){
	if(pdfName.find("Gaus") != std::string::npos && pdfName.size() == 4) genjtpt_[i][j] = lGaus.getNewEnergy(genjtpt_[i][j], genjtpart_[i][j]);
	else if(pdfName.find("Exp") != std::string::npos && pdfName.size() == 3) genjtpt_[i][j] = lExp.getNewEnergy(genjtpt_[i][j], genjtpart_[i][j]);
	else if(pdfName.find("Pois") != std::string::npos && pdfName.size() == 4) genjtpt_[i][j] = lPois.getNewEnergy(genjtpt_[i][j], genjtpart_[i][j]);
      }

      for(Int_t j = 0; j < ngenjt_[i]; ++j){
	if(TMath::Abs(genjteta_[i][j]) > 2.) continue;

	if(genjtpt_[i][j] > leadingJtPt_){
	  subleadingJtPt_ = leadingJtPt_;
	  subleadingJtPhi_ = leadingJtPhi_;
	  subleadingJtFlav_ = leadingJtFlav_;

	  leadingJtPt_ = genjtpt_[i][j];
	  leadingJtPhi_ = genjtphi_[i][j];
	  leadingJtFlav_ = genjtpart_[i][j];
	  //	  if(leadingJtFlav_ == -999) leadingJtFlav_ = genjtpartF_[i][j];
	}
	else if(genjtpt_[i][j] > subleadingJtPt_){
	  subleadingJtPt_ = genjtpt_[i][j];
	  subleadingJtPhi_ = genjtphi_[i][j];
	  subleadingJtFlav_ = genjtpart_[i][j];
	  //	  if(subleadingJtFlav_ == -999) subleadingJtFlav_ = genjtpartF_[i][j];
	}
      }

      if(leadingJtPt_ < 100) continue;
      if(leadingJtPt_ > 126) continue;
      if(subleadingJtPt_ < 30) continue;
      if(TMath::Abs(getDPHI(leadingJtPhi_, subleadingJtPhi_)) < 7.*TMath::Pi()/8.) continue;

      //      if(leadingJtFlav_ == -999 || subleadingJtFlav_ == -999) continue;

      incPthat_h[i]->Fill(pthat_[i]);
      incAJ_h[i]->Fill((leadingJtPt_ - subleadingJtPt_)/(leadingJtPt_ + subleadingJtPt_));
      incXJ_h[i]->Fill(subleadingJtPt_/leadingJtPt_);

      if(TMath::Abs(leadingJtFlav_) < 7 && TMath::Abs(subleadingJtFlav_) < 7){
	qqPthat_h[i]->Fill(pthat_[i]);
	qqAJ_h[i]->Fill((leadingJtPt_ - subleadingJtPt_)/(leadingJtPt_ + subleadingJtPt_));
	qqXJ_h[i]->Fill(subleadingJtPt_/leadingJtPt_);
      }
      else if(TMath::Abs(leadingJtFlav_) < 7 && subleadingJtFlav_ == 21){
	qgPthat_h[i]->Fill(pthat_[i]);
	qgAJ_h[i]->Fill((leadingJtPt_ - subleadingJtPt_)/(leadingJtPt_ + subleadingJtPt_));
	qgXJ_h[i]->Fill(subleadingJtPt_/leadingJtPt_);
      }
      else if(leadingJtFlav_ == 21 && TMath::Abs(subleadingJtFlav_) < 7){
	gqPthat_h[i]->Fill(pthat_[i]);
	gqAJ_h[i]->Fill((leadingJtPt_ - subleadingJtPt_)/(leadingJtPt_ + subleadingJtPt_));
	gqXJ_h[i]->Fill(subleadingJtPt_/leadingJtPt_);
      }
      else if(leadingJtFlav_ == 21 && subleadingJtFlav_ == 21){
	ggPthat_h[i]->Fill(pthat_[i]);
	ggAJ_h[i]->Fill((leadingJtPt_ - subleadingJtPt_)/(leadingJtPt_ + subleadingJtPt_));
	ggXJ_h[i]->Fill(subleadingJtPt_/leadingJtPt_);
      }
      else{
	uPthat_h[i]->Fill(pthat_[i]);
	uAJ_h[i]->Fill((leadingJtPt_ - subleadingJtPt_)/(leadingJtPt_ + subleadingJtPt_));
	uXJ_h[i]->Fill(subleadingJtPt_/leadingJtPt_);
      }
    }
  }

  for(int i = 0; i < nRParam; ++i){
    qqPthat_h[i]->Divide(incPthat_h[i]);
    ggPthat_h[i]->Divide(incPthat_h[i]);
    qgPthat_h[i]->Divide(incPthat_h[i]);
    gqPthat_h[i]->Divide(incPthat_h[i]);
    uPthat_h[i]->Divide(incPthat_h[i]);

    dijetFlav_p[i]->cd();
    gStyle->SetOptStat(0);
    dijetFlav_p[i]->SetTopMargin(0.01);
    dijetFlav_p[i]->SetRightMargin(0.01);
    dijetFlav_p[i]->SetLeftMargin(dijetFlav_p[i]->GetLeftMargin()*1.5);
    dijetFlav_p[i]->SetBottomMargin(dijetFlav_p[i]->GetLeftMargin());
    qqPthat_h[i]->SetFillColor(colors.getColor(0));
    qgPthat_h[i]->SetFillColor(colors.getColor(1));
    gqPthat_h[i]->SetFillColor(colors.getColor(2));
    ggPthat_h[i]->SetFillColor(colors.getColor(3));
    uPthat_h[i]->SetFillColor(colors.getColor(4));

    incPthat_h[i]->SetLineColor(1);
    qqPthat_h[i]->SetLineColor(1);
    qgPthat_h[i]->SetLineColor(1);
    gqPthat_h[i]->SetLineColor(1);
    ggPthat_h[i]->SetLineColor(1);
    uPthat_h[i]->SetLineColor(1);

    ggPthat_h[i]->SetMinimum(0.);
    ggPthat_h[i]->SetMaximum(fracMax);

    ggPthat_h[i]->GetXaxis()->SetTitleOffset(ggPthat_h[i]->GetXaxis()->GetTitleOffset()*1.8);
    ggPthat_h[i]->GetYaxis()->SetTitleOffset(ggPthat_h[i]->GetYaxis()->GetTitleOffset()*1.8);
    ggPthat_h[i]->GetXaxis()->SetTitle("pthat");
    ggPthat_h[i]->GetYaxis()->SetTitle(("Fraction (#color[" + std::to_string(colors.getColor(3)) + "]{GG},#color[" + std::to_string(colors.getColor(2)) + "]{GQ},#color[" + std::to_string(colors.getColor(1)) + "]{QG},#color[" + std::to_string(colors.getColor(0)) + "]{QQ},#color[" + std::to_string(colors.getColor(4)) + "]{Untagged})").c_str());
    ggPthat_h[i]->DrawCopy("HIST E1");

    addWithoutErr(gqPthat_h[i], ggPthat_h[i], 1);
    gqPthat_h[i]->DrawCopy("HIST E1 SAME");

    addWithoutErr(qgPthat_h[i], gqPthat_h[i], 1);
    qgPthat_h[i]->DrawCopy("HIST E1 SAME");


    addWithoutErr(qqPthat_h[i], qgPthat_h[i], 1);
    qqPthat_h[i]->DrawCopy("HIST E1 SAME");

    addWithoutErr(uPthat_h[i], qqPthat_h[i], 1);
    uPthat_h[i]->DrawCopy("HIST E1 SAME");

    qqPthat_h[i]->DrawCopy("HIST E1 SAME");
    qgPthat_h[i]->DrawCopy("HIST E1 SAME");
    gqPthat_h[i]->DrawCopy("HIST E1 SAME");
    ggPthat_h[i]->DrawCopy("HIST E1 SAME");

    gPad->SetLogx();
    gPad->RedrawAxis();

    label_p->DrawLatex(jtPtBins[1], fracYBins[19], ("100<p_{T,1}<126; 30<p_{T,2}; |#eta_{1,2}| < 2.; #Delta#phi > 7#pi/8 (R=0." + std::to_string(int(rParams[i]*10)) + ")").c_str());
    label_p->DrawLatex(jtPtBins[1], fracYBins[18], labelStr.c_str());

    std::string saveName = "pdfDir/dijetFlavR" + std::to_string(int(rParams[i]*10)) + histAppendStr + ".pdf";
    dijetFlav_p[i]->SaveAs(saveName.c_str());

    pthatPastCut_p[i]->cd();
    pthatPastCut_p[i]->SetTopMargin(0.01);
    pthatPastCut_p[i]->SetRightMargin(0.01);
    pthatPastCut_p[i]->SetLeftMargin(pthatPastCut_p[i]->GetLeftMargin()*1.5);
    pthatPastCut_p[i]->SetBottomMargin(pthatPastCut_p[i]->GetLeftMargin());

    incPthat_h[i]->Scale(1./incPthat_h[i]->Integral());

    for(Int_t bIter = 0; bIter < incPthat_h[i]->GetNbinsX(); ++bIter){
      incPthat_h[i]->SetBinContent(bIter+1, incPthat_h[i]->GetBinContent(bIter+1)/incPthat_h[i]->GetBinWidth(bIter+1));
      incPthat_h[i]->SetBinError(bIter+1, incPthat_h[i]->GetBinError(bIter+1)/incPthat_h[i]->GetBinWidth(bIter+1));
    }

    incPthat_h[i]->GetYaxis()->SetTitle("#frac{1}{N_{Dijet}}#times#frac{dN_{dijet}}{dpthat}");
    incPthat_h[i]->SetMaximum(incPthat_h[i]->GetMaximum()*5);
    incPthat_h[i]->GetXaxis()->SetTitleOffset(incPthat_h[i]->GetXaxis()->GetTitleOffset()*1.8);
    incPthat_h[i]->GetYaxis()->SetTitleOffset(incPthat_h[i]->GetYaxis()->GetTitleOffset()*1.8);
    incPthat_h[i]->DrawCopy("HIST E1");

    gPad->SetLogx();
    gPad->SetLogy();
    label_p->DrawLatex(jtPtBins[1], incPthat_h[i]->GetMaximum()*1./2., ("100<p_{T,1}<126; 30<p_{T,2}; |#eta_{1,2}| < 2.; #Delta#phi > 7#pi/8 (R=0." + std::to_string(int(rParams[i]*10)) + ")").c_str());
    label_p->DrawLatex(jtPtBins[1], incPthat_h[i]->GetMaximum()*1./10., labelStr.c_str());

    saveName = "pdfDir/pthatPastCutR" + std::to_string(int(rParams[i]*10)) + histAppendStr + ".pdf";
    pthatPastCut_p[i]->SaveAs(saveName.c_str());

    delete dijetFlav_p[i];
    delete pthatPastCut_p[i];
    delete incPthat_h[i];
    delete qqPthat_h[i];
    delete ggPthat_h[i];
    delete qgPthat_h[i];
    delete gqPthat_h[i];
    delete uPthat_h[i];
  }  


  for(int i = 0; i < nRParam; ++i){
    ajFlav_p[i]->cd();
    gStyle->SetOptStat(0);
    ajFlav_p[i]->SetTopMargin(0.01);
    ajFlav_p[i]->SetRightMargin(0.01);
    ajFlav_p[i]->SetLeftMargin(ajFlav_p[i]->GetLeftMargin()*1.5);
    ajFlav_p[i]->SetBottomMargin(ajFlav_p[i]->GetLeftMargin());

    qqAJ_h[i]->Scale(1./incAJ_h[i]->Integral());
    qgAJ_h[i]->Scale(1./incAJ_h[i]->Integral());
    gqAJ_h[i]->Scale(1./incAJ_h[i]->Integral());
    ggAJ_h[i]->Scale(1./incAJ_h[i]->Integral());
    uAJ_h[i]->Scale(1./incAJ_h[i]->Integral());
    incAJ_h[i]->Scale(1./incAJ_h[i]->Integral());

    for(Int_t bIter = 0; bIter < incAJ_h[i]->GetNbinsX(); ++bIter){
      qqAJ_h[i]->SetBinContent(bIter+1, qqAJ_h[i]->GetBinContent(bIter+1)/qqAJ_h[i]->GetBinWidth(bIter+1)); 
      qqAJ_h[i]->SetBinError(bIter+1, qqAJ_h[i]->GetBinError(bIter+1)/qqAJ_h[i]->GetBinWidth(bIter+1));

      qgAJ_h[i]->SetBinContent(bIter+1, qgAJ_h[i]->GetBinContent(bIter+1)/qgAJ_h[i]->GetBinWidth(bIter+1)); 
      qgAJ_h[i]->SetBinError(bIter+1, qgAJ_h[i]->GetBinError(bIter+1)/qgAJ_h[i]->GetBinWidth(bIter+1));

      gqAJ_h[i]->SetBinContent(bIter+1, gqAJ_h[i]->GetBinContent(bIter+1)/gqAJ_h[i]->GetBinWidth(bIter+1)); 
      gqAJ_h[i]->SetBinError(bIter+1, gqAJ_h[i]->GetBinError(bIter+1)/gqAJ_h[i]->GetBinWidth(bIter+1));

      ggAJ_h[i]->SetBinContent(bIter+1, ggAJ_h[i]->GetBinContent(bIter+1)/ggAJ_h[i]->GetBinWidth(bIter+1)); 
      ggAJ_h[i]->SetBinError(bIter+1, ggAJ_h[i]->GetBinError(bIter+1)/ggAJ_h[i]->GetBinWidth(bIter+1));

      uAJ_h[i]->SetBinContent(bIter+1, uAJ_h[i]->GetBinContent(bIter+1)/uAJ_h[i]->GetBinWidth(bIter+1)); 
      uAJ_h[i]->SetBinError(bIter+1, uAJ_h[i]->GetBinError(bIter+1)/uAJ_h[i]->GetBinWidth(bIter+1)); 

      incAJ_h[i]->SetBinContent(bIter+1, incAJ_h[i]->GetBinContent(bIter+1)/incAJ_h[i]->GetBinWidth(bIter+1)); 
      incAJ_h[i]->SetBinError(bIter+1, incAJ_h[i]->GetBinError(bIter+1)/incAJ_h[i]->GetBinWidth(bIter+1));
    }

    incAJ_h[i]->SetMarkerColor(1);
    incAJ_h[i]->SetMarkerStyle(20);
    incAJ_h[i]->SetMarkerSize(.8);
    qqAJ_h[i]->SetFillColor(colors.getColor(0));
    qgAJ_h[i]->SetFillColor(colors.getColor(1));
    gqAJ_h[i]->SetFillColor(colors.getColor(2));
    ggAJ_h[i]->SetFillColor(colors.getColor(3));
    uAJ_h[i]->SetFillColor(colors.getColor(4));

    incAJ_h[i]->SetLineColor(1);
    qqAJ_h[i]->SetLineColor(1);
    qgAJ_h[i]->SetLineColor(1);
    gqAJ_h[i]->SetLineColor(1);
    ggAJ_h[i]->SetLineColor(1);
    uAJ_h[i]->SetLineColor(1);

    incAJ_h[i]->GetXaxis()->SetTitleOffset(incAJ_h[i]->GetXaxis()->GetTitleOffset()*1.8);
    incAJ_h[i]->GetYaxis()->SetTitleOffset(incAJ_h[i]->GetYaxis()->GetTitleOffset()*1.8);
    incAJ_h[i]->GetXaxis()->SetTitle("A_{J}");
    incAJ_h[i]->GetYaxis()->SetTitle(("#frac{1}{N_{Dijet}} #times #frac{dN_{Dijet}}{dA_{J}} (#color[" + std::to_string(colors.getColor(3)) + "]{GG},#color[" + std::to_string(colors.getColor(2)) + "]{GQ},#color[" + std::to_string(colors.getColor(1)) + "]{QG},#color[" + std::to_string(colors.getColor(0)) + "]{QQ},#color[" + std::to_string(colors.getColor(4)) + "]{Untagged})").c_str());
    incAJ_h[i]->SetMinimum(0.0);
    incAJ_h[i]->SetMaximum(8.0);
    incAJ_h[i]->DrawCopy("E1 P");

    ggAJ_h[i]->DrawCopy("HIST E1 SAME");

    addWithoutErr(gqAJ_h[i], ggAJ_h[i], 1);
    gqAJ_h[i]->DrawCopy("HIST E1 SAME");

    addWithoutErr(qgAJ_h[i], gqAJ_h[i], 1);
    qgAJ_h[i]->DrawCopy("HIST E1 SAME");

    addWithoutErr(qqAJ_h[i], qgAJ_h[i], 1);
    qqAJ_h[i]->DrawCopy("HIST E1 SAME");

    addWithoutErr(uAJ_h[i], qqAJ_h[i], 1);
    uAJ_h[i]->DrawCopy("HIST E1 SAME");

    incAJ_h[i]->DrawCopy("E1 P SAME");
    qqAJ_h[i]->DrawCopy("HIST E1 SAME");
    qgAJ_h[i]->DrawCopy("HIST E1 SAME");
    gqAJ_h[i]->DrawCopy("HIST E1 SAME");
    ggAJ_h[i]->DrawCopy("HIST E1 SAME");

    gPad->RedrawAxis();

    label_p->DrawLatex(ajBins[1], incAJ_h[i]->GetMaximum()*9./10., ("100<p_{T,1}<126; 30<p_{T,2}; |#eta_{1,2}| < 2.; #Delta#phi > 7#pi/8 (R=0." + std::to_string(int(rParams[i]*10)) + ")").c_str());
    label_p->DrawLatex(ajBins[1], incAJ_h[i]->GetMaximum()*8./10., labelStr.c_str());

    std::string saveName = "pdfDir/ajFlavR" + std::to_string(int(rParams[i]*10)) + histAppendStr + ".pdf";
    ajFlav_p[i]->SaveAs(saveName.c_str());

    ajFlavFrac_p[i]->cd();
    ajFlavFrac_p[i]->SetTopMargin(0.01);
    ajFlavFrac_p[i]->SetRightMargin(0.01);
    ajFlavFrac_p[i]->SetLeftMargin(ajFlavFrac_p[i]->GetLeftMargin()*1.5);
    ajFlavFrac_p[i]->SetBottomMargin(ajFlavFrac_p[i]->GetLeftMargin());
    
    addWithoutErr(uAJ_h[i], qqAJ_h[i], -1);
    addWithoutErr(qqAJ_h[i], qgAJ_h[i], -1);
    addWithoutErr(qgAJ_h[i], gqAJ_h[i], -1);
    addWithoutErr(gqAJ_h[i], ggAJ_h[i], -1);



    qqAJ_h[i]->Divide(incAJ_h[i]);
    qgAJ_h[i]->Divide(incAJ_h[i]);
    gqAJ_h[i]->Divide(incAJ_h[i]);
    ggAJ_h[i]->Divide(incAJ_h[i]);
    uAJ_h[i]->Divide(incAJ_h[i]);
    incAJ_h[i]->Divide(incAJ_h[i]);

    incAJ_h[i]->SetMaximum(fracMax);
    incAJ_h[i]->SetMinimum(0.);

    incAJ_h[i]->GetYaxis()->SetTitle(("Fraction #frac{1}{N_{Dijet}} #times #frac{dN_{Dijet}}{dA_{J}} (#color[" + std::to_string(colors.getColor(3)) + "]{GG},#color[" + std::to_string(colors.getColor(2)) + "]{GQ},#color[" + std::to_string(colors.getColor(1)) + "]{QG},#color[" + std::to_string(colors.getColor(0)) + "]{QQ},#color[" + std::to_string(colors.getColor(4)) + "]{Untagged})").c_str());

    incAJ_h[i]->DrawCopy("E1 P");

    addWithoutErr(gqAJ_h[i], ggAJ_h[i], 1);
    addWithoutErr(qgAJ_h[i], gqAJ_h[i], 1);
    addWithoutErr(qqAJ_h[i], qgAJ_h[i], 1);
    addWithoutErr(uAJ_h[i], qqAJ_h[i], 1);

    uAJ_h[i]->DrawCopy("HIST SAME E1");
    qqAJ_h[i]->DrawCopy("HIST SAME E1");
    qgAJ_h[i]->DrawCopy("HIST SAME E1");
    gqAJ_h[i]->DrawCopy("HIST SAME E1");
    ggAJ_h[i]->DrawCopy("HIST SAME E1");

    gPad->RedrawAxis();

    label_p->DrawLatex(ajBins[1], fracYBins[19], ("100<p_{T,1}<126; 30<p_{T,2}; |#eta_{1,2}| < 2.; #Delta#phi > 7#pi/8 (R=0." + std::to_string(int(rParams[i]*10)) + ")").c_str());
    label_p->DrawLatex(ajBins[1], fracYBins[18], labelStr.c_str());

    saveName = "pdfDir/ajFlavFracR" + std::to_string(int(rParams[i]*10)) + histAppendStr + ".pdf";
    ajFlavFrac_p[i]->SaveAs(saveName.c_str());

    delete ajFlav_p[i];
    delete ajFlavFrac_p[i];
    delete incAJ_h[i];
    delete qqAJ_h[i];
    delete ggAJ_h[i];
    delete qgAJ_h[i];
    delete gqAJ_h[i];
    delete uAJ_h[i];
  }  

  for(int i = 0; i < nRParam; ++i){
    xjFlav_p[i]->cd();
    gStyle->SetOptStat(0);
    xjFlav_p[i]->SetTopMargin(0.01);
    xjFlav_p[i]->SetRightMargin(0.01);
    xjFlav_p[i]->SetLeftMargin(xjFlav_p[i]->GetLeftMargin()*1.5);
    xjFlav_p[i]->SetBottomMargin(xjFlav_p[i]->GetLeftMargin());

    qqXJ_h[i]->Scale(1./incXJ_h[i]->Integral());
    qgXJ_h[i]->Scale(1./incXJ_h[i]->Integral());
    gqXJ_h[i]->Scale(1./incXJ_h[i]->Integral());
    ggXJ_h[i]->Scale(1./incXJ_h[i]->Integral());
    uXJ_h[i]->Scale(1./incXJ_h[i]->Integral());
    incXJ_h[i]->Scale(1./incXJ_h[i]->Integral());

    for(Int_t bIter = 0; bIter < incXJ_h[i]->GetNbinsX(); ++bIter){
      incXJ_h[i]->SetBinContent(bIter+1, incXJ_h[i]->GetBinContent(bIter+1)/incXJ_h[i]->GetBinWidth(bIter+1)); 
      incXJ_h[i]->SetBinError(bIter+1, incXJ_h[i]->GetBinError(bIter+1)/incXJ_h[i]->GetBinWidth(bIter+1));

      qqXJ_h[i]->SetBinContent(bIter+1, qqXJ_h[i]->GetBinContent(bIter+1)/qqXJ_h[i]->GetBinWidth(bIter+1)); 
      qqXJ_h[i]->SetBinError(bIter+1, qqXJ_h[i]->GetBinError(bIter+1)/qqXJ_h[i]->GetBinWidth(bIter+1));

      qgXJ_h[i]->SetBinContent(bIter+1, qgXJ_h[i]->GetBinContent(bIter+1)/qgXJ_h[i]->GetBinWidth(bIter+1)); 
      qgXJ_h[i]->SetBinError(bIter+1, qgXJ_h[i]->GetBinError(bIter+1)/qgXJ_h[i]->GetBinWidth(bIter+1));

      gqXJ_h[i]->SetBinContent(bIter+1, gqXJ_h[i]->GetBinContent(bIter+1)/gqXJ_h[i]->GetBinWidth(bIter+1)); 
      gqXJ_h[i]->SetBinError(bIter+1, gqXJ_h[i]->GetBinError(bIter+1)/gqXJ_h[i]->GetBinWidth(bIter+1));

      ggXJ_h[i]->SetBinContent(bIter+1, ggXJ_h[i]->GetBinContent(bIter+1)/ggXJ_h[i]->GetBinWidth(bIter+1)); 
      ggXJ_h[i]->SetBinError(bIter+1, ggXJ_h[i]->GetBinError(bIter+1)/ggXJ_h[i]->GetBinWidth(bIter+1));

      uXJ_h[i]->SetBinContent(bIter+1, uXJ_h[i]->GetBinContent(bIter+1)/uXJ_h[i]->GetBinWidth(bIter+1)); 
      uXJ_h[i]->SetBinError(bIter+1, uXJ_h[i]->GetBinError(bIter+1)/uXJ_h[i]->GetBinWidth(bIter+1)); 
    }

    incXJ_h[i]->SetMarkerColor(1);
    incXJ_h[i]->SetMarkerStyle(20);
    incXJ_h[i]->SetMarkerSize(.8);
    qqXJ_h[i]->SetFillColor(colors.getColor(0));
    qgXJ_h[i]->SetFillColor(colors.getColor(1));
    gqXJ_h[i]->SetFillColor(colors.getColor(2));
    ggXJ_h[i]->SetFillColor(colors.getColor(3));
    uXJ_h[i]->SetFillColor(colors.getColor(4));

    incXJ_h[i]->SetLineColor(1);
    qqXJ_h[i]->SetLineColor(1);
    qgXJ_h[i]->SetLineColor(1);
    gqXJ_h[i]->SetLineColor(1);
    ggXJ_h[i]->SetLineColor(1);
    uXJ_h[i]->SetLineColor(1);
  
    incXJ_h[i]->GetXaxis()->SetTitleOffset(incXJ_h[i]->GetXaxis()->GetTitleOffset()*1.8);
    incXJ_h[i]->GetYaxis()->SetTitleOffset(incXJ_h[i]->GetYaxis()->GetTitleOffset()*1.8);
    incXJ_h[i]->GetXaxis()->SetTitle("x_{J}");
    incXJ_h[i]->GetYaxis()->SetTitle(("#frac{1}{N_{Dijet}} #times #frac{dN_{Dijet}}{dx_{J}} (#color[" + std::to_string(colors.getColor(3)) + "]{GG},#color[" + std::to_string(colors.getColor(2)) + "]{GQ},#color[" + std::to_string(colors.getColor(1)) + "]{QG},#color[" + std::to_string(colors.getColor(0)) + "]{QQ},#color[" + std::to_string(colors.getColor(4)) + "]{Untagged})").c_str());
    incXJ_h[i]->SetMinimum(0.0);
    incXJ_h[i]->SetMaximum(xjMax[i]);
    incXJ_h[i]->DrawCopy("E1 P");

    ggXJ_h[i]->DrawCopy("HIST E1 SAME");

    addWithoutErr(gqXJ_h[i], ggXJ_h[i], 1);
    addWithoutErr(qgXJ_h[i], gqXJ_h[i], 1);
    addWithoutErr(qqXJ_h[i], qgXJ_h[i], 1);
    addWithoutErr(uXJ_h[i], qqXJ_h[i], 1);

    incXJ_h[i]->DrawCopy("E1 P SAME");
    uXJ_h[i]->DrawCopy("HIST E1 SAME");
    qqXJ_h[i]->DrawCopy("HIST E1 SAME");
    qgXJ_h[i]->DrawCopy("HIST E1 SAME");
    gqXJ_h[i]->DrawCopy("HIST E1 SAME");
    ggXJ_h[i]->DrawCopy("HIST E1 SAME");

    if(int(10*rParams[i]) == 4){
      xjATLASR0p4_h->DrawCopy("P SAME");
      drawSyst(xjFlav_p[i], xjATLASR0p4_h, xjYValErr, true, colors.getColor(5));

      xjATLASR0p4_pp_h->DrawCopy("P SAME");
      drawSyst(xjFlav_p[i], xjATLASR0p4_pp_h, xjYValErr_pp, true, colors.getColor(6));
    }

    gPad->RedrawAxis();

    Double_t yLabelBins[20+1];
    getLinBins(0.0, xjMax[i], 20, yLabelBins);

    label_p->DrawLatex(xjBins[1], yLabelBins[19], ("100<p_{T,1}<126; 30<p_{T,2}; |#eta_{1,2}| < 2.; #Delta#phi > 7#pi/8 (R=0." + std::to_string(int(rParams[i]*10)) + ")").c_str());
    label_p->DrawLatex(xjBins[1], yLabelBins[18], labelStr.c_str());
    label_p->DrawLatex(xjBins[1], yLabelBins[17], "PYTHIA 6 (MSEL=1)");

    std::string saveName = "pdfDir/xjFlavR" + std::to_string(int(rParams[i]*10)) + histAppendStr + ".pdf";
    xjFlav_p[i]->SaveAs(saveName.c_str());

    xjFlavFrac_p[i]->cd();
    xjFlavFrac_p[i]->SetTopMargin(0.01);
    xjFlavFrac_p[i]->SetRightMargin(0.01);
    xjFlavFrac_p[i]->SetLeftMargin(xjFlavFrac_p[i]->GetLeftMargin()*1.5);
    xjFlavFrac_p[i]->SetBottomMargin(xjFlavFrac_p[i]->GetLeftMargin());

    addWithoutErr(uXJ_h[i], qqXJ_h[i], -1);
    addWithoutErr(qqXJ_h[i], qgXJ_h[i], -1); 
    addWithoutErr(qgXJ_h[i], gqXJ_h[i], -1);
    addWithoutErr(gqXJ_h[i], ggXJ_h[i], -1);

    qqXJ_h[i]->Divide(incXJ_h[i]);
    qgXJ_h[i]->Divide(incXJ_h[i]);
    gqXJ_h[i]->Divide(incXJ_h[i]);
    ggXJ_h[i]->Divide(incXJ_h[i]);
    uXJ_h[i]->Divide(incXJ_h[i]);
    incXJ_h[i]->Divide(incXJ_h[i]);

    incXJ_h[i]->SetMaximum(fracMax);
    incXJ_h[i]->SetMinimum(0.);

    incXJ_h[i]->GetYaxis()->SetTitle(("Fraction #frac{1}{N_{Dijet}} #times #frac{dN_{Dijet}}{dA_{J}} (#color[" + std::to_string(colors.getColor(3)) + "]{GG},#color[" + std::to_string(colors.getColor(2)) + "]{GQ},#color[" + std::to_string(colors.getColor(1)) + "]{QG},#color[" + std::to_string(colors.getColor(0)) + "]{QQ},#color[" + std::to_string(colors.getColor(4)) + "]{Untagged})").c_str());


    addWithoutErr(gqXJ_h[i], ggXJ_h[i], 1);
    addWithoutErr(qgXJ_h[i], gqXJ_h[i], 1);
    addWithoutErr(qqXJ_h[i], qgXJ_h[i], 1); 
    addWithoutErr(uXJ_h[i], qqXJ_h[i], 1);

    incXJ_h[i]->DrawCopy("E1 P");
    uXJ_h[i]->DrawCopy("HIST SAME E1");
    qqXJ_h[i]->DrawCopy("HIST SAME E1");
    qgXJ_h[i]->DrawCopy("HIST SAME E1");
    gqXJ_h[i]->DrawCopy("HIST SAME E1");
    ggXJ_h[i]->DrawCopy("HIST SAME E1");

    gPad->RedrawAxis();

    label_p->DrawLatex(xjBins[1], fracYBins[19], ("100<p_{T,1}<126; 30<p_{T,2}; |#eta_{1,2}| < 2.; #Delta#phi > 7#pi/8 (R=0." + std::to_string(int(rParams[i]*10)) + ")").c_str());
    label_p->DrawLatex(xjBins[1], fracYBins[18], labelStr.c_str());

    saveName = "pdfDir/xjFlavFracR" + std::to_string(int(rParams[i]*10)) + histAppendStr + ".pdf";
    xjFlavFrac_p[i]->SaveAs(saveName.c_str());

    delete xjFlav_p[i];
    delete xjFlavFrac_p[i];
    delete incXJ_h[i];
    delete qqXJ_h[i];
    delete ggXJ_h[i];
    delete qgXJ_h[i];
    delete gqXJ_h[i];
    delete uXJ_h[i];
  }  

  delete xjATLASR0p4_h;

  inFile_p->Close();
  delete inFile_p;

  delete label_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./plotDijetFlavVsPthat.exe <inFileName>" << std::endl;
    return 1;
  }

  const Int_t nMean = 8;
  const Int_t nSigma = 6;
 
  const Double_t means[nMean] = {0.5, 0.45, 0.4, 0.35, 0.3, 0.2, 0.1, 0.};
  const Double_t sigmas[nSigma] = {0.25,0.2,0.15,0.1,0.05, 0.01};

  //  const Int_t nPoisson = 11;
  //  const Int_t nPoisson2 = 11;
  //  const Int_t counts[nPoisson] = {1,3,5,7,9,11,13,15,17,19,21};
  //  const Double_t ePerCount[nPoisson] = {.001,.002,.004,.006,.008,.01,.02, .04, .06, .08, .1};
  
  int retVal = 0;
  retVal += plotDijetFlavVsPthat(argv[1], "", {0, 0, 0, 0}, {0, 0});
  
  for(Int_t i = 0; i < nMean; ++i){
    if(i < 4) continue;

    for(Int_t j = 0; j < nMean; ++j){

      for(Int_t k = 0; k < nSigma; ++k){
	for(Int_t l = 0; l < nSigma; ++l){	  

	  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", {means[i], sigmas[k], means[j], sigmas[l]}, {0, 0});
	}
      }

    }
  }

  /*
  for(Int_t i = 0; i < nPoisson2; ++i){
    for(Int_t j = 0; j < nPoisson; ++j){
      retVal += plotDijetFlavVsPthat(argv[1], "Pois", {ePerCount[j], 0., 0., 0.}, {counts[i], 9*counts[i]/4});
    }
  }
  */
  /*
  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", .1, .1, .2, .2);
  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", .2, .2, .4, .4);
  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", .1, .01, .2, .01);
  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", .2, .01, .4, .01);
  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", .05, .05, .2, .01);
  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", .05, .05, .4, .01);
  retVal += plotDijetFlavVsPthat(argv[1], "Gaus", .05, .05, .6, .01);
  retVal += plotDijetFlavVsPthat(argv[1], "Exp", .1, 0., .2, 0.);
  retVal += plotDijetFlavVsPthat(argv[1], "Exp", .2, 0., .4, 0.);
  */
  return retVal;
}
