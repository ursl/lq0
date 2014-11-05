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
  void   makeAll(int bitmask = 0);
  void   treeAnalysis(std::string cds1, std::string cds2, std::string cds3); 
  void   normOverlay(std::string f1, std::string f2); 
  void   genMass(std::string type = "lq_pair", int offset = 29, int nplot = 4);    
  void   overlayAll();

  void   genMassPlots(std::string dir = "single"); 
  
  void   bookHist(std::string name); 
  void   setupTree(TTree *t); 
  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  void   loopFunction1(); 
  void   loopFunction2(); 
  void   loopFunction3(); 
  void   loopFunction4(); 
  
  void  optimizeCuts(std::vector<std::string> samples, std::string dir, double lumi = 20., int nevts = -1); 
  void  optAnalysis(int mode = 1, std::string filename = "opt.root", std::string treename = "opt_single");
  void  displayOptimization(std::string file, std::string tree); 


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
