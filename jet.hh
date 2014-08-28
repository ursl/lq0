#ifndef JET_H
#define JET_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"

struct jet {
  TLorentzVector p4; 
  Jet *pJet;
};



#endif
