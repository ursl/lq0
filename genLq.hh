#ifndef GENLQ_H
#define GENLQ_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"


struct genLq {
  int q; 
  TLorentzVector p4LQ, p4L, p4Q, p4J, p4K, p4LJ;
  GenParticle *pLQ, *pL, *pQ, *pK; 
  Jet *pJ;
};

#endif
