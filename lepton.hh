#ifndef LEPTON_H
#define LEPTON_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"

struct lepton {
  TLorentzVector fP4; 
  int q; 
  Muon *fpMuon;
  Electron *fpElectron;
};

#endif
