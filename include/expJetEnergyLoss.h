#ifndef EXPJETENERGYLOSS_H
#define EXPJETENERGYLOSS_H

#include "TRandom3.h"
#include "TMath.h"

class expJetEnergyLoss{
 public:
  expJetEnergyLoss(Double_t qMean, Double_t gMean);
  ~expJetEnergyLoss();
  Double_t getNewEnergy(Double_t jtPt, Int_t jtFlavor);

  Double_t qM,gM;
  TRandom3* randGen_p;
};

expJetEnergyLoss::expJetEnergyLoss(Double_t qMean, Double_t gMean)
{
  qM = qMean;
  gM = gMean;
  randGen_p = new TRandom3(0);
  return;
}

expJetEnergyLoss::~expJetEnergyLoss(){delete randGen_p; return;}

Double_t expJetEnergyLoss::getNewEnergy(Double_t jtPt, Int_t jtFlavor)
{
  Double_t lostE = jtPt;
  if(TMath::Abs(jtFlavor) < 7) lostE *= randGen_p->Exp(qM);
  else lostE *= randGen_p->Exp(gM);

  if(lostE < 0) lostE = 0;
  jtPt < lostE ? jtPt = 0 : jtPt -= lostE;
  return jtPt;
}

#endif
