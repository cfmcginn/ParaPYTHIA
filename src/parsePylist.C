#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

int parsePylist(const std::string inPyFileName, const std::string outFileName = "outFile.root")
{
  const std::string numStr = "0123456789.";
  const Int_t nMaxPart = 10000;

  std::string pthatStr = "";
  Float_t pthat_=-1;
  Int_t nPart_;
  Float_t pt_[nMaxPart];
  Float_t phi_[nMaxPart];
  Float_t eta_[nMaxPart];
  Int_t pdgid_[nMaxPart];

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
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("genTree", "genTree");

  outTree_p->Branch("pthat", &pthat_, "pthat/F");
  outTree_p->Branch("nPart", &nPart_, "nPart/I");
  outTree_p->Branch("pt", pt_, "pt[nPart]/F");
  outTree_p->Branch("phi", phi_, "phi[nPart]/F");
  outTree_p->Branch("eta", eta_, "eta[nPart]/F");
  outTree_p->Branch("pdgid", pdgid_, "pdgid[nPart]/I");

  outTree_p->Branch("nQG", &nQG_, "nQG/I");
  outTree_p->Branch("qgPt", qgPt_, "qgPt[nQG]/F");
  outTree_p->Branch("qgPhi", qgPhi_, "qgPhi[nQG]/F");
  outTree_p->Branch("qgEta", qgEta_, "qgEta[nQG]/F");
  outTree_p->Branch("qgPDGID", qgPDGID_, "qgPDGID[nQG]/I");

  outTree_p->Branch("nQGF", &nQGF_, "nQGF/I");
  outTree_p->Branch("qgFPt", qgFPt_, "qgFPt[nQGF]/F");
  outTree_p->Branch("qgFPhi", qgFPhi_, "qgFPhi[nQGF]/F");
  outTree_p->Branch("qgFEta", qgFEta_, "qgFEta[nQGF]/F");
  outTree_p->Branch("qgFPDGID", qgFPDGID_, "qgFPDGID[nQGF]/I");
  
  std::ifstream file(inPyFileName);
  std::string getStr;

  int evtCnt = 0;

  while(std::getline(file, getStr)){
    //    std::cout << getStr << std::endl;
    if(getStr.size() == 0) continue;
    if(getStr.find("==") != std::string::npos) continue;
    if(getStr.find("Event listing") != std::string::npos) continue;
    if(getStr.find("sum") != std::string::npos) continue;
    if(getStr.find("particle") != std::string::npos) continue;
    if(getStr.find("PYSTAT:  Statistics on Number of Events and Cross-sections") != std::string::npos) break;

    if(getStr.find("pre-event") != std::string::npos){
      if(evtCnt%500 == 0) std::cout << "Event count: " << evtCnt << std::endl;
      ++evtCnt;
      if(pthat_ > 0) outTree_p->Fill();
      nPart_=0;
      nQG_=0;
      nQGF_=0;

      std::getline(file, getStr);
      pthatStr = getStr;
      pthat_ = std::stof(getStr);
      continue;
    }

    if(pthat_ > 0){
      for(int i = 0; i < (int)numStr.size(); ++i){
	std::string searchStr = numStr.substr(i,1) + "-";
	std::string repStr = numStr.substr(i,1) + " -";

	while(getStr.find(searchStr) != std::string::npos){
	  getStr.replace(getStr.find(searchStr),searchStr.size(),repStr);
	}
      }
      
      while(getStr.substr(0,1).find(" ") != std::string::npos){getStr.replace(0,1,"");}
      
      std::vector<std::string> strVect;
      while(getStr.find(" ") != std::string::npos){
	std::string tempString = getStr.substr(0, getStr.find(" "));
	strVect.push_back(tempString);
	getStr.replace(0, getStr.find(" ")+1, "");
	while(getStr.substr(0,1).find(" ") != std::string::npos){getStr.replace(0,1,"");}
      }
      if(getStr.size() != 0) strVect.push_back(getStr);
      if(strVect.size() == 11) strVect.erase(strVect.begin()+2);

      if(std::stoi(strVect.at(2)) != 1 && std::stoi(strVect.at(2)) != 21 && std::stoi(strVect.at(2)) != 12){strVect.clear(); continue;}

      if(std::stoi(strVect.at(2)) == 21 && std::stoi(strVect.at(3)) == 2212) continue;

      //      std::cout << " check " << strVect.at(0) << ", " <<  strVect.at(1) << ", " <<  strVect.at(2) << ", " <<  strVect.at(3) << std::endl;

      Int_t id = std::stoi(strVect.at(3));
      if(TMath::Abs(id) == 12) continue;
      if(TMath::Abs(id) == 14) continue;
      if(TMath::Abs(id) == 16) continue;
      
      Float_t px = std::stof(strVect.at(5));
      Float_t py = std::stof(strVect.at(6));
      Float_t pz = std::stof(strVect.at(7));
      Float_t E = std::stof(strVect.at(8));
      TLorentzVector t(px,py,pz,E);
      
      if(TMath::Abs(px) < .005 && TMath::Abs(py) < .005) continue;
      if(TMath::Abs(t.Eta()) > 5.) continue;
      if(t.Pt() < .1) continue;
      
      if(std::stoi(strVect.at(2)) == 1){
	if(t.Eta() > 100){
	  std::cout << "Possible error: " << pthat_ << ", " << pthatStr << ", " << strVect.at(5) << ", " << strVect.at(6) << ", " << strVect.at(7) << ", " << strVect.at(8) << std::endl;
	}
	
	pt_[nPart_] = t.Pt();
	phi_[nPart_] = t.Phi();
	eta_[nPart_] = t.Eta();
	pdgid_[nPart_] = std::stoi(strVect.at(3));
	++nPart_;
	strVect.clear();
      }
      else if(std::stoi(strVect.at(2)) == 21){
        qgPt_[nQG_] = t.Pt();
        qgPhi_[nQG_] = t.Phi();
        qgEta_[nQG_] = t.Eta();
        qgPDGID_[nQG_] = std::stoi(strVect.at(3));
        ++nQG_;
        strVect.clear();
      }
      else{
        qgFPt_[nQGF_] = t.Pt();
        qgFPhi_[nQGF_] = t.Phi();
        qgFEta_[nQGF_] = t.Eta();
        qgFPDGID_[nQGF_] = std::stoi(strVect.at(3));
        ++nQGF_;
        strVect.clear();
      }
    }
  }
  outTree_p->Fill();
  file.close();

  outFile_p->cd();
  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;
  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./parsePylist.exe <inPyFileName> <outFileName>" << std::endl;
    return 1;
  }
 
  int retVal = 0;
  if(argc == 2) retVal += parsePylist(argv[1]);
  else if(argc == 3) retVal += parsePylist(argv[1], argv[2]);
  return retVal;
}
