#ifndef POISSONJETENERGYLOSS_H
#define POISSONJETENERGYLOSS_H

#include "TRandom3.h"
#include "TMath.h"

class poissonJetEnergyLoss{
 public:
  poissonJetEnergyLoss(Int_t meanCollQ, Int_t meanCollG, Double_t ePerColl);
  ~poissonJetEnergyLoss();
  Double_t getNewEnergy(Double_t jtPt, Int_t jtFlavor);

  Int_t meanQ;
  Int_t meanG;
  Double_t ePer;
  TRandom3* randGen_p;
};

poissonJetEnergyLoss::poissonJetEnergyLoss(Int_t meanCollQ, Int_t meanCollG, Double_t ePerColl)
{
  meanQ = meanCollQ;
  meanG = meanCollG;
  ePer = ePerColl;

  randGen_p = new TRandom3(0);

  return;
}

poissonJetEnergyLoss::~poissonJetEnergyLoss(){delete randGen_p; return;}

Double_t poissonJetEnergyLoss::getNewEnergy(Double_t jtPt, Int_t jtFlavor)
{
  Double_t lostE = jtPt;
  Int_t counts = 0;

  if(TMath::Abs(jtFlavor) < 7) counts = randGen_p->Poisson(meanQ);
  else counts = randGen_p->Poisson(meanG);

  while(counts > 0){
    lostE -= lostE*ePer;
    --counts;
  }

  if(lostE < 0) lostE = 0;
  return lostE;
}

#endif
