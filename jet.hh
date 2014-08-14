#ifndef JET_H
#define JET_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"

struct jet {
  TLorentzVector fP4; 
  Jet *fpJet;
};



#endif
