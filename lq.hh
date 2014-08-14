#ifndef LQ_H
#define LQ_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"

#include "lepton.hh"
#include "jet.hh"

struct lq {
  TLorentzVector fP4; 
  int q; 
  jet *fpJet; 
  lepton *fpLepton; 
};

#endif
