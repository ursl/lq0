#ifndef LEPTON_H
#define LEPTON_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"

struct lepton {
  TLorentzVector p4; 
  int q; 
  Muon *pMuon;
  Electron *pElectron;
};

#endif
