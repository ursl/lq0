#ifndef REDTREEDATA_H
#define REDTREEDATA_H

#define NLQ 10

struct redTreeData {
  int channel; 
  double w8; 

  int ngen;
  // LQ 
  double gm[NLQ], gpt[NLQ], geta[NLQ], gphi[NLQ];
  // gen-level lepton charge (q)
  int glq[NLQ], gtm[NLQ]; 
  // lepton + gen jet
  double gmlj[NLQ]; 
  // gen-level particles
  double glpt[NLQ], gleta[NLQ], glphi[NLQ], 
    gqpt[NLQ], gqeta[NLQ], gqphi[NLQ], 
    gjpt[NLQ], gjeta[NLQ], gjphi[NLQ], 
    gkpt[NLQ], gketa[NLQ], gkphi[NLQ],
    gipt[NLQ], gieta[NLQ], giphi[NLQ];
  int gkq[NLQ]; 

  // -- reco 
  int nrec;
  // LQ
  double m[NLQ], pt[NLQ], eta[NLQ], phi[NLQ]; 
  // rec-level lepton charge (q)
  int lq[NLQ], tm[NLQ];  
  // reco-level particles
  double lpt[NLQ], leta[NLQ], lphi[NLQ];  // lepton
  double jpt[NLQ], jeta[NLQ], jphi[NLQ];  // jet
  double kpt[NLQ], keta[NLQ], kphi[NLQ];  // bachelor lepton k
  int kq[NLQ]; 

  // other information since this information can depend on the lq combination, it is per LQ 
  // (this is certainly redundant for pair production ...)
  double st[NLQ]; 
  double mll[NLQ];
  double mljmin[NLQ];
  int    type[NLQ];
};

#endif
