#ifndef PLOTLQ_h
#define PLOTLQ_h

#include "plotClass.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotLq: public plotClass {

public :
                 plotLq(std::string dir = "hpt0", std::string files = "plotLq.files", std::string setup = "m");
  virtual        ~plotLq();

  virtual void   loadFiles(std::string afiles);
  void           splitType(std::string stype, std::string &name, double &mass, double &lambda);

  // -- Main analysis methods 
  virtual void   makeAll(int bitmask = 0);
  virtual void   treeAnalysis(); 
  virtual void   normOverlay(std::string f1, std::string f2); 

  void           overlayAll();
  
  virtual void   bookHist(std::string name); 
  virtual void   setupTree(TTree *t); 
  virtual void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  virtual void   loopFunction1(); 
  virtual void   loopFunction2(); 
  
  void           optimizePairCuts(std::string sg, std::string bg, double lumi = 20.); 


private: 

  bool fPair;
  bool fGoodEvent, fGoodCandLQp, fGoodCandLQn; 

  struct redTreeData fRtd; 

  // -- optimization
  int                   fOptMode; 
  std::vector<selpoint> fSelPoints;
  
  // ----------------------------------------------------------------------
  ClassDef(plotLq,1) 

};

#endif
