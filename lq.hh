#ifndef LQ_H
#define LQ_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"

struct lq {
  TLorentzVector p4; 
  int q, tm; 
  int idxJ, idxL, idxK;
  double st, mljmin, mll;
};

#endif
