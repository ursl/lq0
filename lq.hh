#ifndef LQ_H
#define LQ_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"

struct lq {
  TLorentzVector p4; 
  int q; 
  int idxJ, idxL, idxK;
};

#endif
