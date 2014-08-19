#ifndef PLOTLQ_h
#define PLOTLQ_h


#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include "dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotLq: public TObject {

public :
  plotLq(std::string dir = "hpt0", std::string files = "plotLq.files", std::string setup = "m");
  ~plotLq();

  void   loadFiles(std::string afiles);
  TFile* loadFile(std::string afiles);
  void   splitType(std::string stype, std::string &name, double &mass, double &lambda);


  // -- Main analysis methods 
  void makeAll(int bitmask = 0);
  void treeAnalysis(); 

  void overlayAll();
  void overlay(TH1* h1, std::string f1, TH1 *h2, std::string f2, bool legend = true);
  void overlay(std::string f1, std::string h1name, std::string f2, std::string h2name, bool legend = true);

  void bookHist(std::string name); 
  TTree* getTree(std::string ds); 
  void setupTree(TTree *t); 
  void loopOverTree(TTree *t, int nevts = -1, int nstart = 0); 
  void candAnalysis(); 
  void loopFunction(); 

  void cd(std::string dataset) {fDS[dataset]->cd("");}
  void replaceAll(std::string &sInput, const std::string &oldString, const std::string &newString);
  void newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  void makeCanvas(int i = 3);
  void normHist(TH1 *, std::string ds="", int method = NONORM); 

private: 

  enum HistNorm {NONORM,     // do not touch the normalization
		 SOMETHING,  // normalize to what is given in fNorm
		 UNITY,      // normalize all to 1
		 XSECTION,   // the resulting histograms will be cross sections
		 LUMI        // according to the number provided in fLumi
  };
  double fNorm, fLumi;

  int    fVerbose; 
  double fEpsilon; 

  std::string fDirectory, fSetup, fSuffix;   

  int    NBINS; 

  bool fPair;
  bool fGoodEvent, fGoodCandLQp, fGoodCandLQn; 

  struct redTreeData fRtd; 

  int             fOptMode; 
  std::map<std::string, TH1*> fHists;

  // -- datasets (files and associated information)
  std::map<std::string, dataset*> fDS; 
  // -- current dataset for analysis
  std::string fCds; 

  // -- Display utilities
  int fFont; 
  double fSize; 
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *legg;
  TLegendEntry *legge;


  // ----------------------------------------------------------------------
  ClassDef(plotLq,1) 

};


#endif
