#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "fastjet/ClusterSequence.hh"

#include "include/etaPhiFunc.h"

int clusterGenToJet(const std::string inGenFile, const std::string outFileName = "outFile_Clustered.root")
{
  const Int_t nMaxPart = 10000;
  Float_t pthat_;
  Int_t nPart_;
  Float_t pt_[nMaxPart];
  Float_t phi_[nMaxPart];
  Float_t eta_[nMaxPart];

  const Int_t nMaxQG_ = 500;
  Int_t nQG_;
  Float_t qgPt_[nMaxQG_];
  Float_t qgPhi_[nMaxQG_];
  Float_t qgEta_[nMaxQG_];
  Int_t qgPDGID_[nMaxQG_];

  Int_t nQGF_;
  Float_t qgFPt_[nMaxQG_];
  Float_t qgFPhi_[nMaxQG_];
  Float_t qgFEta_[nMaxQG_];
  Int_t qgFPDGID_[nMaxQG_];

  /*
  const int nRParam = 5;
  const double rParams[nRParam] = {0.2, 0.3, 0.4, 0.5, 0.6};
  */

  const int nRParam = 2;
  const double rParams[nRParam] = {0.3, 0.4};

  const Int_t nMaxJets = 10000;
  Int_t ngenjt_[nRParam];
  Float_t genjtpt_[nRParam][nMaxJets];
  Float_t genjtphi_[nRParam][nMaxJets];
  Float_t genjteta_[nRParam][nMaxJets];
  Int_t genjtpart_[nRParam][nMaxJets];
  Float_t genjtpartpt_[nRParam][nMaxJets];
  Int_t genjtpartH_[nRParam][nMaxJets];
  Float_t genjtpartHpt_[nRParam][nMaxJets];
  Int_t genjtpartF_[nRParam][nMaxJets];
  Float_t genjtpartFpt_[nRParam][nMaxJets];

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p[nRParam];

  for(int i = 0; i < nRParam; ++i){
    std::string treeName = "ak" + std::to_string(int(rParams[i]*10)) + "GenJetTree";
    outTree_p[i] = new TTree(treeName.c_str(), "");

    outTree_p[i]->Branch("pthat", &pthat_, "pthat/F");
    outTree_p[i]->Branch("ngen", &ngenjt_[i], "ngen/I");
    outTree_p[i]->Branch("genpt", genjtpt_[i], "genpt[ngen]/F");
    outTree_p[i]->Branch("genphi", genjtphi_[i], "genphi[ngen]/F");
    outTree_p[i]->Branch("geneta", genjteta_[i], "geneta[ngen]/F");
    outTree_p[i]->Branch("genpart", genjtpart_[i], "genpart[ngen]/I");
    outTree_p[i]->Branch("genpartpt", genjtpartpt_[i], "genpartpt[ngen]/F");
    outTree_p[i]->Branch("genpartH", genjtpartH_[i], "genpartH[ngen]/I");
    outTree_p[i]->Branch("genpartHpt", genjtpartHpt_[i], "genpartHpt[ngen]/F");
    outTree_p[i]->Branch("genpartF", genjtpartF_[i], "genpartF[ngen]/I");
    outTree_p[i]->Branch("genpartFpt", genjtpartFpt_[i], "genpartFpt[ngen]/F");
  }

  TFile* inFile_p = new TFile(inGenFile.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("genTree");

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("pthat", 1);
  inTree_p->SetBranchStatus("nPart", 1);
  inTree_p->SetBranchStatus("pt", 1);
  inTree_p->SetBranchStatus("phi", 1);
  inTree_p->SetBranchStatus("eta", 1);
  inTree_p->SetBranchStatus("nQG", 1);
  inTree_p->SetBranchStatus("qgPt", 1);
  inTree_p->SetBranchStatus("qgPhi", 1);
  inTree_p->SetBranchStatus("qgEta", 1);
  inTree_p->SetBranchStatus("qgPDGID", 1);
  inTree_p->SetBranchStatus("nQGF", 1);
  inTree_p->SetBranchStatus("qgFPt", 1);
  inTree_p->SetBranchStatus("qgFPhi", 1);
  inTree_p->SetBranchStatus("qgFEta", 1);
  inTree_p->SetBranchStatus("qgFPDGID", 1);

  inTree_p->SetBranchAddress("pthat", &pthat_);
  inTree_p->SetBranchAddress("nPart", &nPart_);
  inTree_p->SetBranchAddress("pt", pt_);
  inTree_p->SetBranchAddress("phi", phi_);
  inTree_p->SetBranchAddress("eta", eta_);
  inTree_p->SetBranchAddress("nQG", &nQG_);
  inTree_p->SetBranchAddress("qgPt", qgPt_);
  inTree_p->SetBranchAddress("qgPhi", qgPhi_);
  inTree_p->SetBranchAddress("qgEta", qgEta_);
  inTree_p->SetBranchAddress("qgPDGID", qgPDGID_);
  inTree_p->SetBranchAddress("nQGF", &nQGF_);
  inTree_p->SetBranchAddress("qgFPt", qgFPt_);
  inTree_p->SetBranchAddress("qgFPhi", qgFPhi_);
  inTree_p->SetBranchAddress("qgFEta", qgFEta_);
  inTree_p->SetBranchAddress("qgFPDGID", qgFPDGID_);

  const Int_t nEntries = inTree_p->GetEntries();

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%1000 == 0) std::cout << "Entry " << entry << "/" << nEntries << std::endl;
    
    inTree_p->GetEntry(entry);

    std::vector<fastjet::PseudoJet> particles;
    for(int i = 0; i < nPart_; ++i){
      TLorentzVector temp;
      temp.SetPtEtaPhiM(pt_[i], eta_[i], phi_[i], 0);
      particles.push_back(fastjet::PseudoJet(temp.Px(), temp.Py(), temp.Pz(), temp.E()));
    }

    for(Int_t j = 0; j < nRParam; ++j){
      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParams[j]);
      fastjet::ClusterSequence cs(particles, jetDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      
      //      if(jets.at(0).pt() < 50) continue;
      ngenjt_[j] = 0;

      const Int_t nUsedQG = nQG_;
      Bool_t isUsedQG[nUsedQG];
      for(Int_t i = 0; i < nQG_; ++i){isUsedQG[i] = false;}

      for(Int_t l = 0; l < nQG_; ++l){
	for(Int_t k = l+1; k < nQG_; ++k){
	  if(qgPt_[l] < qgPt_[k]){
	    Float_t tempPt = qgPt_[k];
	    Float_t tempPhi = qgPhi_[k];
	    Float_t tempEta = qgEta_[k];
	    Int_t tempPDGID = qgPDGID_[k];

	    qgPt_[k] = qgPt_[l];
	    qgPhi_[k] = qgPhi_[l];
	    qgEta_[k] = qgEta_[l];
	    qgPDGID_[k] = qgPDGID_[l];

	    qgPt_[l] = tempPt;
            qgPhi_[l] = tempPhi;
            qgEta_[l] = tempEta;
            qgPDGID_[l] = tempPDGID;
	  }
	}
      }

      const Int_t nUsedQGF = nQGF_;
      Bool_t isUsedQGF[nUsedQGF];
      for(Int_t i = 0; i < nQGF_; ++i){isUsedQGF[i] = false;}

      for(Int_t l = 0; l < nQGF_; ++l){
	for(Int_t k = l+1; k < nQGF_; ++k){
	  if(qgFPt_[l] < qgFPt_[k]){
	    Float_t tempPt = qgFPt_[k];
	    Float_t tempPhi = qgFPhi_[k];
	    Float_t tempEta = qgFEta_[k];
	    Int_t tempPDGID = qgFPDGID_[k];

	    qgFPt_[k] = qgFPt_[l];
	    qgFPhi_[k] = qgFPhi_[l];
	    qgFEta_[k] = qgFEta_[l];
	    qgFPDGID_[k] = qgFPDGID_[l];

	    qgFPt_[l] = tempPt;
            qgFPhi_[l] = tempPhi;
            qgFEta_[l] = tempEta;
            qgFPDGID_[l] = tempPDGID;
	  }
	}
      }

      for(unsigned int i = 0; i < jets.size(); ++i){
	if(jets.at(i).pt() < 15) break;
	//      if(TMath::Abs(jets.at(i).eta()) > 3.) continue;
	
	genjtpt_[j][ngenjt_[j]] = jets.at(i).pt();
	genjtphi_[j][ngenjt_[j]] = jets.at(i).phi_std();
	genjteta_[j][ngenjt_[j]] = jets.at(i).eta();

	genjtpart_[j][ngenjt_[j]] = -999;
	genjtpartpt_[j][ngenjt_[j]] = 0;

	genjtpartH_[j][ngenjt_[j]] = -999;
	genjtpartHpt_[j][ngenjt_[j]] = 0;

	genjtpartF_[j][ngenjt_[j]] = -999;
	genjtpartFpt_[j][ngenjt_[j]] = 0;

	for(Int_t l = 0; l < nQG_; ++l){
	  if(isUsedQG[l]) continue;

	  if(getDR(jets.at(i).eta(), jets.at(i).phi(), qgEta_[l], qgPhi_[l]) < rParams[j]){
	    genjtpart_[j][ngenjt_[j]] = qgPDGID_[l];
	    genjtpartpt_[j][ngenjt_[j]] = qgPt_[l];

	    genjtpartH_[j][ngenjt_[j]] = qgPDGID_[l];
	    genjtpartHpt_[j][ngenjt_[j]] = qgPt_[l];

	    isUsedQG[l] = true;
	    break;
	  }
	}
      
	for(Int_t l = 0; l < nQGF_; ++l){
	  if(isUsedQGF[l]) continue;
	  if(qgFPt_[l] < 10.) continue;

	  if(getDR(jets.at(i).eta(), jets.at(i).phi(), qgFEta_[l], qgFPhi_[l]) < rParams[j]){
	    if(genjtpart_[j][ngenjt_[j]] == -999){
	      genjtpart_[j][ngenjt_[j]] = qgFPDGID_[l];
	      genjtpartpt_[j][ngenjt_[j]] = qgFPt_[l];
	    }

	    genjtpartF_[j][ngenjt_[j]] = qgFPDGID_[l];
	    genjtpartFpt_[j][ngenjt_[j]] = qgFPt_[l];
	    isUsedQGF[l] = true;
	    break;
	  }
	}
      

	++ngenjt_[j];
      }

      outTree_p[j]->Fill();
    }

    particles.clear();
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(int i = 0; i < nRParam; ++i){
    outTree_p[i]->Write("", TObject::kOverwrite);
    delete outTree_p[i];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./clusterGenToJet.exe <inGenFile> <outFileName>" << std::endl;
    return 1;
  }
 
  int retVal = 0;
  if(argc == 2) retVal += clusterGenToJet(argv[1]);
  else if(argc == 3) retVal += clusterGenToJet(argv[1], argv[2]);
  return retVal;
}
