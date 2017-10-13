#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TNamed.h"

int extractWeightAndErr(const std::string inFileName1, const std::string inFileName2, const int thresh)
{
  Float_t pthat1_;
  Float_t pthat2_;

  TFile* inFile1_p = new TFile(inFileName1.c_str(), "READ");
  TTree* inTree1_p = (TTree*)inFile1_p->Get("ak3GenJetTree");
  inTree1_p->SetBranchStatus("*", 0);
  inTree1_p->SetBranchStatus("pthat", 1);

  inTree1_p->SetBranchAddress("pthat", &pthat1_);

  TFile* inFile2_p = new TFile(inFileName2.c_str(), "READ");
  TTree* inTree2_p = (TTree*)inFile2_p->Get("ak3GenJetTree");
  inTree2_p->SetBranchStatus("*", 0);
  inTree2_p->SetBranchStatus("pthat", 1);

  inTree2_p->SetBranchAddress("pthat", &pthat2_);

  
  const Double_t binLow = thresh-.00001;
  const Double_t binHigh = thresh+100;
  const Int_t nBins = (binHigh - binLow)/5;

  TH1F* noWeightSpectra_h = new TH1F("noWeightSpectra_h", ";p_{T} Hat;Events", nBins, binLow, binHigh);
  TH1F* weightSpectra_h = new TH1F("weightSpectra_h", ";p_{T} Hat;Events", nBins, binLow, binHigh);
  TH1F* dividedSpectra_h = new TH1F("dividedSpectra_h", ";p_{T} Hat;Events", nBins, binLow, binHigh);

  const int nEntries1 = inTree1_p->GetEntries();
  const int nEntries2 = inTree2_p->GetEntries();

  for(int entry = 0; entry < nEntries1; ++entry){
    inTree1_p->GetEntry(entry);

    if(pthat1_ < binLow) continue;
    if(pthat1_ >= binHigh) continue;

    noWeightSpectra_h->Fill(pthat1_);
    weightSpectra_h->Fill(pthat1_);
    dividedSpectra_h->Fill(pthat1_);
  }

  for(int entry = 0; entry < nEntries2; ++entry){
    inTree2_p->GetEntry(entry);

    if(pthat2_ < binLow) continue;
    if(pthat2_ >= binHigh) continue;

    weightSpectra_h->Fill(pthat2_);
  }

  noWeightSpectra_h->Sumw2();
  weightSpectra_h->Sumw2();
  dividedSpectra_h->Sumw2();
  
  TFile* outFile_p = new TFile("extractFile.root", "RECREATE");

  noWeightSpectra_h->Write("", TObject::kOverwrite);
  weightSpectra_h->Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  outFile_p = new TFile("extractFile.root", "UPDATE");

  dividedSpectra_h->Divide(weightSpectra_h);
  
  TF1* f1_p = new TF1("f1_p", "[0]", thresh, thresh+100);

  dividedSpectra_h->Fit("f1_p", "M N", "", thresh, thresh+100);
  dividedSpectra_h->Write("", TObject::kOverwrite);

  TNamed fitPar0("fitPar0", std::to_string(f1_p->GetParameter(0)));
  TNamed fitPar0Err("fitPar0Err", std::to_string(f1_p->GetParError(0)));
  
  fitPar0.Write("", TObject::kOverwrite);
  fitPar0Err.Write("", TObject::kOverwrite);

  delete f1_p;

  outFile_p->Close();
  delete outFile_p;

  delete noWeightSpectra_h;
  delete weightSpectra_h;
  delete dividedSpectra_h;

  inFile2_p->Close();
  delete inFile2_p;

  inFile1_p->Close();
  delete inFile1_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "./Usage: ./extractWeightAndErr.exe <inFileName1> <inFileName2> <thresh>" << std::endl;
    return 1;
  }


  int retVal = 0;
  retVal += extractWeightAndErr(argv[1], argv[2], std::stoi(argv[3]));
  return retVal;
}
