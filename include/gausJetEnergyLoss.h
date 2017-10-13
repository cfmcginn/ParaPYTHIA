#ifndef GAUSJETENERGYLOSS_H
#define GAUSJETENERGYLOSS_H

#include "TRandom3.h"
#include "TMath.h"

class gausJetEnergyLoss{
 public:
  gausJetEnergyLoss(Double_t qMean, Double_t qSigma, Double_t gMean, Double_t gSigma);
  ~gausJetEnergyLoss();
  Double_t getNewEnergy(Double_t jtPt, Int_t jtFlavor);

  Double_t qM,qS,gM,gS;
  TRandom3* randGen_p;
};

gausJetEnergyLoss::gausJetEnergyLoss(Double_t qMean, Double_t qSigma, Double_t gMean, Double_t gSigma)
{
  qM = qMean;
  qS = qSigma;
  gM = gMean;
  gS = gSigma;

  randGen_p = new TRandom3(0);

  return;
}

gausJetEnergyLoss::~gausJetEnergyLoss(){delete randGen_p; return;}

Double_t gausJetEnergyLoss::getNewEnergy(Double_t jtPt, Int_t jtFlavor)
{
  Double_t lostE = jtPt;
  if(TMath::Abs(jtFlavor) < 7) lostE *= randGen_p->Gaus(qM,qS);
  else lostE *= randGen_p->Gaus(gM,gS);

  if(lostE < 0) lostE = 0;
  jtPt < lostE ? jtPt = 0 : jtPt -= lostE;
  return jtPt;
}

#endif
