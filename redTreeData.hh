#ifndef REDTREEDATA_H
#define REDTREEDATA_H

struct redTreeData {
  int type, channel; 
  double w8; 

  int ngen;
  // LQ 
  double gm[2], gpt[2], geta[2], gphi[2];
  // gen-level lepton charge (q)
  int glq[2], gtm[2]; 
  // lepton + gen jet
  double gmlj[2]; 
  // gen-level particles
  double glpt[2], gleta[2], glphi[2], 
    gqpt[2], gqeta[2], gqphi[2], 
    gjpt[2], gjeta[2], gjphi[2], 
    gkpt[2], gketa[2], gkphi[2],
    gipt[2], gieta[2], giphi[2];

  // -- reco 
  int nrec;
  // LQ
  double m[2], pt[2], eta[2], phi[2]; 
  // rec-level lepton charge (q)
  int lq[2], tm[2];  
  // reco-level particles
  double lpt[2], leta[2], lphi[2];  // lepton
  double jpt[2], jeta[2], jphi[2];  // jet
  double kpt[2], keta[2], kphi[2];  // bachelor lepton k

  // other information
  double st; 
  double mll;
  double mljmin;

};

#endif
