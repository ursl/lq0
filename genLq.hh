#ifndef GENLQ_H
#define GENLQ_H

#include <TLorentzVector.h>
#include "classes/DelphesClasses.h"


struct genLq {
  int q, tm; 
  TLorentzVector p4LQ, p4L, p4Q, p4J, p4K, p4I, p4LJ;
  GenParticle *pLQ, *pL, *pQ, *pK; 
  Jet *pJ, *pI;
};

#endif
