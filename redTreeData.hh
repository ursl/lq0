#ifndef REDTREEDATA_H
#define REDTREEDATA_H

struct redTreeData {
  int type, channel; 
  double w8; 

  int ngen;
  // LQ 
  double gm[4], gpt[4], geta[4], gphi[4];
  // gen-level lepton charge (q)
  int glq[4], gtm[4]; 
  // lepton + gen jet
  double gmlj[4]; 
  // gen-level particles
  double glpt[4], gleta[4], glphi[4], 
    gqpt[4], gqeta[4], gqphi[4], 
    gjpt[4], gjeta[4], gjphi[4], 
    gkpt[4], gketa[4], gkphi[4],
    gipt[4], gieta[4], giphi[4];
  int gkq[4]; 

  // -- reco 
  int nrec;
  // LQ
  double m[4], pt[4], eta[4], phi[4]; 
  // rec-level lepton charge (q)
  int lq[4], tm[4];  
  // reco-level particles
  double lpt[4], leta[4], lphi[4];  // lepton
  double jpt[4], jeta[4], jphi[4];  // jet
  double kpt[4], keta[4], kphi[4];  // bachelor lepton k
  int kq[4]; 

  // other information
  double st; 
  double mll;
  double mljmin;

};

#endif
