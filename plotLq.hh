#ifndef PLOTLQ_h
#define PLOTLQ_h

#include "plotClass.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotLq: public plotClass {

public :
                 plotLq(std::string dir = "results", std::string files = "lq.files", std::string setup = "m");
  virtual        ~plotLq();

  virtual void   loadFiles(std::string afiles);
  void           splitType(std::string stype, std::string &name, double &mass, double &lambda);

  // -- Main analysis methods 
  virtual void   makeAll(int bitmask = 0);
  virtual void   treeAnalysis(std::string cds1, std::string cds2, std::string cds3); 
  virtual void   normOverlay(std::string f1, std::string f2); 

  void           overlayAll();
  
  virtual void   bookHist(std::string name); 
  virtual void   setupTree(TTree *t); 
  virtual void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  virtual void   loopFunction1(); 
  virtual void   loopFunction2(); 
  virtual void   loopFunction3(); 
  virtual void   loopFunction4(); 
  
  void           optimizeCuts(std::string dir, std::string sg, std::string bg, 
			      double lumi = 20., int nevts = -1); 
  void           displayOptimization(std::string file, std::string tree); 


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
