#include "plotLq.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"

#include "util.hh"

ClassImp(plotLq)

using namespace std; 

// ----------------------------------------------------------------------
plotLq::plotLq(string dir,  string files, string setup): plotClass(dir, files, setup) {

  loadFiles(files);

  fLumi = 100.;

  fHistFile = TFile::Open(Form("%s/plotLq.root", dir.c_str()), "RECREATE"); 
}


// ----------------------------------------------------------------------
plotLq::~plotLq() {
  fHistFile->Close();
}

// ----------------------------------------------------------------------
void plotLq::bookHist(string name) {

  char cds[200];
  sprintf(cds, "%s", fCds.c_str());
  
  // -- start gen level
  char hist[200];
  sprintf(hist, "%s_%s", "ngen", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 5, 0, 5.))); 
  setHistTitles(fHists[hist], fDS[name], "n_{Gen}(LQ)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptl0", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(l0)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptl1", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 120, 0, 1200.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(l1)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptq0", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(q0)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptq1", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 120, 0, 1200.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(q1)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptlq", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 120, 0, 1200.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(lq)", "Entries/bin");

  sprintf(hist, "%s_%s", "etaq0", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -4., 4.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta_{T}(q0)", "Entries/bin");

  sprintf(hist, "%s_%s", "etaq1", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -4., 4.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta_{T}(q1)", "Entries/bin");

  sprintf(hist, "%s_%s", "etalq", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -4., 4.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta_{T}(lq)", "Entries/bin");


  sprintf(hist, "%s_%s", "dRll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, 0, 6.))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta R(l,l)", "Entries/bin");


  sprintf(hist, "%s_%s", "dRql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, 0, 6.))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta R(q,l)", "Entries/bin");

  sprintf(hist, "%s_%s", "dRxql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, 0, 6.))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta RX(q,l)", "Entries/bin");

  sprintf(hist, "%s_%s", "dFql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3.15, 3.15))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta #phi(q,l)", "Entries/bin");

  sprintf(hist, "%s_%s", "dFxql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3.15, 3.15))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta #phiX(q,l)", "Entries/bin");




  sprintf(hist, "%s_%s", "dPtll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, -100., 900.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(l0) - pT(l1) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "sPtll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0., 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(l0) + pT(l1) [GeV]", "Entries/bin");


  sprintf(hist, "%s_%s", "dPtqq", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, -100., 900.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(q0) - pT(q1) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "sPtqq", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0., 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(q0) + pT(q1) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "mxql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "mX(q,l) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "mqXlY", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "m(qX,lY) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "mqXlX", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "m(qX,lX) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "mll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "m(l,k) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "ptll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "pt(l0 + l1) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "ptqq", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0., 2000.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(q0 + q1) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "st", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 200, 0., 5000.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(l0) + pT(l1) + pT(q0) + pT(q1) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "dst", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, -100., 1400.))); 
  setHistTitles(fHists[hist], fDS[name], "(pT(l0) + pT(q0)) -  (pT(l1)  + pT(q1)) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "dpt", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0., 500.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(l0+q0) - (pT(l1_q1) [GeV]", "Entries/bin");

  // -- end gen level


  // -- m
  fHists.insert(make_pair(Form("m_%s", name.c_str()), 
			  new TH1D(Form("m_%s", name.c_str()), Form("m_%s", name.c_str()), 100, 0, 2000.))); 
  setTitles(fHists[Form("m_%s", name.c_str())], "m [GeV]", "Entries/bin");
  setHist(fHists[Form("m_%s", name.c_str())], fDS[name]);

  // -- st
//   fHists.insert(make_pair(Form("st_%s", name.c_str()), 
// 			  new TH1D(Form("st_%s", name.c_str()), Form("st_%s", name.c_str()), 35, 0, 3500.))); 
//   setTitles(fHists[Form("st_%s", name.c_str())], "S_{T} [GeV]", "Entries/bin");
//   setHist(fHists[Form("st_%s", name.c_str())], fDS[name]);

//   // -- mll
//   fHists.insert(make_pair(Form("mll_%s", name.c_str()), 
// 			  new TH1D(Form("mll_%s", name.c_str()), Form("mll_%s", name.c_str()), 50, 0, 1500.))); 
//   setTitles(fHists[Form("mll_%s", name.c_str())], "m_{l l} [GeV]", "Entries/bin");
//   setHist(fHists[Form("mll_%s", name.c_str())], fDS[name]);

  // -- mljetmin
  fHists.insert(make_pair(Form("mljmin_%s", name.c_str()), 
			  new TH1D(Form("mljmin_%s", name.c_str()), Form("mljmin_%s", name.c_str()), 60, 0, 1500.))); 
  setTitles(fHists[Form("mljmin_%s", name.c_str())], "m_{l jet}^{min} [GeV]", "Entries/bin");
  setHist(fHists[Form("mljmin_%s", name.c_str())], fDS[name]);
  
  fHists.insert(make_pair(Form("pt_%s", name.c_str()), 
			  new TH1D(Form("pt_%s", name.c_str()), Form("pt_%s", name.c_str()), 100, 0, 1000.))); 
  setHist(fHists[Form("pt_%s", name.c_str())], fDS[name]);
  setTitles(fHists[Form("pt_%s", name.c_str())], "p_{T} [GeV]", "Entries/bin");
}


// ----------------------------------------------------------------------
void plotLq::makeAll(int bitmask) {
  if (bitmask & 0x1) {
    genMass();
    c1->SaveAs("genMass-lq_pair.pdf");
    genMass("lq_pair_down", 37, 3);
    c1->SaveAs("genMass-lq_pair_down.pdf");
    genMass("lq_pair_up", 37, 3);
    c1->SaveAs("genMass-lq_pair_up.pdf");
    genMass("lq_single_up", 37, 3);
    c1->SaveAs("genMass-lq_single_up.pdf");
    genMass("lq_single_down", 37, 3);
    c1->SaveAs("genMass-lq_single_down.pdf");
  }

  if (bitmask & 0x2) {
    vector<string> samples;
    samples.push_back("lq_single_up_37");
    samples.push_back("dy_el");
    optimizeCuts(samples, "opt.root", 100.);
  }
      
    
}

// ----------------------------------------------------------------------
void plotLq::treeAnalysis(string cds1, string cds2, string cds3) {
  if (0) {
    // -- pair
    fPair = false;
    fCds = cds1; 
    bookHist(fCds); 
    TTree *tp = getTree(fCds, "pair"); 
    setupTree(tp); 
    loopOverTree(tp, 3, 20000); 
    
    // -- single
    fPair = false;
    fCds = cds2;
    bookHist(fCds); 
    TTree *ts = getTree(fCds, "pair"); 
    setupTree(ts); 
    loopOverTree(ts, 3, 20000); 
  }

  if (1) {
    // -- pair
    fPair = false;
    fCds = cds1; 
    bookHist(fCds); 
    TTree *tp = getTree(fCds, "lq"); 
    setupTree(tp); 
    loopOverTree(tp, 4, 20000); 
    
    // -- single
    fPair = false;
    fCds = cds2;
    bookHist(fCds); 
    TTree *ts = getTree(fCds, "lq"); 
    setupTree(ts); 
    loopOverTree(ts, 4, 20000); 

    // -- single
    fPair = false;
    fCds = cds3;
    bookHist(fCds); 
    TTree *tb = getTree(fCds, "lq"); 
    setupTree(tb); 
    loopOverTree(tb, 4, 20000); 
  }

  char hist[200];
  int ipad(1); 
  zone(5, 5);
  sprintf(hist, "ngen");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptl0");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptl1");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dPtll");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "sPtll");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptll");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptq0");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptq1");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dPtqq");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "sPtqq");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "st");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "mqXlY");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "mqXlX");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "mll");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dst");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptlq");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dpt");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "etaq0");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "etaq1");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "etalq");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dFql");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dFxql");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);


  c0->cd(++ipad); 
  sprintf(hist, "dRll");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);


  c0->cd(++ipad); 
  sprintf(hist, "dRql");
  overlay(fHists[Form("%s_%s", hist, cds1.c_str())], cds1, fHists[Form("%s_%s", hist, cds2.c_str())], cds2, UNITY, false, false);



}

// ----------------------------------------------------------------------
void plotLq::normOverlay(std::string ds0, std::string ds1) {
  //  string ds0 = "dy_pair";
  fPair = true;
  fCds = ds0; 
  bookHist(ds0); 
  TTree *t = getTree(ds0); 
  setupTree(t); 
  loopOverTree(t, 1); 

  //  string ds1 = "lq_pair_01"; 
  fPair = true;
  fCds = ds1; 
  bookHist(ds1); 
  t = getTree(ds1); 
  setupTree(t);
  loopOverTree(t, 1); 

  string hist("m");
  overlay(fHists[Form("%s_%s", hist.c_str(), ds0.c_str())], ds0, fHists[Form("%s_%s", hist.c_str(), ds1.c_str())], ds1); 
  c0->SaveAs(Form("%s-%s-%s.pdf", hist.c_str(), ds0.c_str(), ds1.c_str()));

  hist = "pt";
  overlay(fHists[Form("%s_%s", hist.c_str(), ds0.c_str())], ds0, fHists[Form("%s_%s", hist.c_str(), ds1.c_str())], ds1); 
  c0->SaveAs(Form("%s-%s-%s.pdf", hist.c_str(), ds0.c_str(), ds1.c_str()));

  hist = "st";
  overlay(fHists[Form("%s_%s", hist.c_str(), ds0.c_str())], ds0, fHists[Form("%s_%s", hist.c_str(), ds1.c_str())], ds1); 
  c0->SaveAs(Form("%s-%s-%s.pdf", hist.c_str(), ds0.c_str(), ds1.c_str()));

  hist = "mljmin";
  overlay(fHists[Form("%s_%s", hist.c_str(), ds0.c_str())], ds0, fHists[Form("%s_%s", hist.c_str(), ds1.c_str())], ds1); 
  c0->SaveAs(Form("%s-%s-%s.pdf", hist.c_str(), ds0.c_str(), ds1.c_str()));
  
  hist = "mll";
  overlay(fHists[Form("%s_%s", hist.c_str(), ds0.c_str())], ds0, fHists[Form("%s_%s", hist.c_str(), ds1.c_str())], ds1); 
  c0->SaveAs(Form("%s-%s-%s.pdf", hist.c_str(), ds0.c_str(), ds1.c_str()));

}

// ----------------------------------------------------------------------
void plotLq::overlayAll() {
  // -- simple overlays
  c0->cd(1); 
}

// ----------------------------------------------------------------------
void plotLq::displayOptimization(string file, string tree) {


}

// ----------------------------------------------------------------------
void plotLq::optimizeCuts(vector<string> samples, string fname, double lumi, int nevts) {

  TFile *f(0); 
  cout << "open file " << Form("%s/%s", 
			       fDirectory.c_str(), fname.c_str()) << " RECREATE" << endl;
  f = TFile::Open(Form("%s/%s", 
		       fDirectory.c_str(), fname.c_str()), "RECREATE");     


  // -- ST
  static const double stArr[] = { 500., 700., 1000., 1200., 1400., 1700., 2200., 2500., 3000.
  };
  vector<double> stCuts(stArr, stArr + sizeof(stArr)/sizeof(stArr[0]));
  
  // -- m(l,l)
  static const double mllArr[] = {100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1200., 1500., 2000.
  };
  vector<double> mllCuts(mllArr, mllArr + sizeof(mllArr)/sizeof(mllArr[0]));
  
  // -- m(l,j)
  static const double mljArr[] = {100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200.
  };
  vector<double> mljCuts(mljArr, mljArr + sizeof(mljArr)/sizeof(mljArr[0]));
  
  for (unsigned int i = 0; i < stCuts.size(); ++i) {
    for (unsigned int j = 0; j < mljCuts.size(); ++j) {
      for (unsigned int k = 0; k < mllCuts.size(); ++k) {
	selpoint s; 
	s.fLargerThan.push_back(make_pair(&fRtd.st[0], stCuts[i])); 
	s.fLargerThan.push_back(make_pair(&fRtd.mljmin[0], mljCuts[j])); 
	s.fLargerThan.push_back(make_pair(&fRtd.mll[0], mllCuts[k])); 
	fSelPoints.push_back(s); 
      }
    }
  }

  string dir("pair"); 
  if (string::npos != samples[0].find("single")) {
    dir = "single";
  }

  if (!dir.compare("pair")) {
    fPair = true;
  } else {
    fPair = false;
  }
  
  // -- signal
  fCds = samples[0];
  fOptMode = 0; 
  bookHist(fCds); 
  double nevtsScale(1.); 
  TTree *ts = getTree(fCds, dir); 
  if (nevts > -1) {
    nevtsScale = fDS[fCds]->fLumi*nevts/ts->GetEntries();
  }
  double sgScale = lumi/fDS[fCds]->fLumi;
  setupTree(ts); 
  loopOverTree(ts, 2, nevts); 

  // -- background
  vector<double> bgScale; 
  for (unsigned int i = 1; i < samples.size(); ++i) {
    fCds = samples[i];
    fOptMode = i; 
    bookHist(fCds); 
    bgScale.push_back(lumi/fDS[fCds]->fLumi);
    TTree *tb = getTree(fCds, dir); 
    setupTree(tb); 
    loopOverTree(tb, 2, nevts); 
  }

  f->cd();
  TTree *t = new TTree(Form("opt_%s", dir.c_str()), Form("opt_%s", dir.c_str()));
  double s, b, s0, b0; 
  double st, mll, mlj;
  double ssb, sb;
  t->Branch("s",    &s,     "s/D");
  t->Branch("b",    &b,     "b/D");
  t->Branch("s0",   &s0,    "s0/D");
  t->Branch("b0",   &b0,    "b0/D");
  t->Branch("ssb",  &ssb,   "ssb/D");
  t->Branch("sb",   &sb,    "sb/D");

  t->Branch("st",   &st,    "st/D");
  t->Branch("mll",  &mll,   "mll/D");
  t->Branch("mlj",  &mlj,   "mlj/D");
  
  // FIXME this still has only ONE background!
  for (unsigned int i = 0; i < fSelPoints.size(); ++i) {
    s0  = fSelPoints[i].fCnt[0];
    b0  = fSelPoints[i].fCnt[1];
    s   = sgScale * fSelPoints[i].fCnt[0];
    b   = bgScale[0] * fSelPoints[i].fCnt[1];
    ssb = (s+b>0? s/TMath::Sqrt(s+b):0.);
    sb  = (b>0?   s/TMath::Sqrt(b)  :0.);
    st  = fSelPoints[i].fLargerThan[0].second;
    mll = fSelPoints[i].fLargerThan[1].second;
    mlj = fSelPoints[i].fLargerThan[2].second;
    t->Fill();
  }

  t->Write();
  f->Close();
}


// ----------------------------------------------------------------------
void plotLq::optAnalysis(int mode, string filename, string treename) {
  TFile *f(0); 
  if (filename.compare("")) {
    f = TFile::Open(Form("%s/%s", fDirectory.c_str(), filename.c_str())); 
  }
  TTree *t = (TTree*)f->Get(treename.c_str()); 
  if (0 == t) {
    cout << "no tree with name " << treename << " found in file " << filename << endl;
    return;
  }

  // -- setup tree
  double ssb, sb, s, b, b1; 
  double st, mll, mlj;
  t->SetBranchAddress("ssb",   &ssb);
  t->SetBranchAddress("s",     &s);
  t->SetBranchAddress("b",     &b);
  t->SetBranchAddress("st",    &st);
  t->SetBranchAddress("mll",   &mll);
  t->SetBranchAddress("mlj",   &mlj);

  TH1D *hst = new TH1D("hst", "", 300, 0., 3001.);
  TH1D *hmll = new TH1D("hmll", "", 400, 0., 2001.);
  TH1D *hmlj = new TH1D("hmlj", "", 400, 0., 2000.);
  
  // -- loop over tree
  int ibin(0);
  int nentries = Int_t(t->GetEntries());
  double fom(0.); 
  for (int jentry = 0; jentry < nentries; jentry++) {
    t->GetEntry(jentry);
    if (1 == mode) fom = ssb;
    if (fom > hst->GetBinContent(hst->FindBin(st))) {
      hst->SetBinContent(hst->FindBin(st), fom); 
    }

    if (fom > hmll->GetBinContent(hst->FindBin(mll))) {
      hmll->SetBinContent(hmll->FindBin(mll), fom); 
    }

    if (fom > hmlj->GetBinContent(hst->FindBin(mlj))) {
      hmlj->SetBinContent(hmlj->FindBin(mlj), fom); 
    }
  }

  zone(2,2);
  hst->Draw();
  c0->cd(2);
  hmll->Draw();
  c0->cd(3);
  hmlj->Draw();
  c0->SaveAs(Form("%s/opt-%d.pdf", fDirectory.c_str(), mode)); 
    

}


// ----------------------------------------------------------------------
void plotLq::loopFunction1() {

  // -- cuts
  fGoodEvent   = false; 
  
  if (fPair) {
    if (fRtd.m[0] > 0. && fRtd.m[1] > 0.) fGoodEvent = true;
    if (fRtd.st[0] < 685.)       fGoodEvent = false;
    if (fRtd.mll[0] < 150.)      fGoodEvent = false;
    if (fRtd.mljmin[0] < 155.) fGoodEvent = false;
  } else {
    if (fRtd.m[0] > 0.) fGoodEvent = true;
  }

  char cds[200];
  sprintf(cds, "%s", fCds.c_str());

  if (fGoodEvent) { 
    //    fHists[Form("st_%s", cds)]->Fill(fRtd.st); 
    //    fHists[Form("mll_%s", cds)]->Fill(fRtd.mll); 
    fHists[Form("mljmin_%s", cds)]->Fill(fRtd.mljmin[0]); 
      
    for (int i = 0; i < fRtd.nrec; ++i) {
      fHists[Form("m_%s", cds)]->Fill(fRtd.m[i]); 
      fHists[Form("pt_%s", cds)]->Fill(fRtd.pt[i]); 
    }
  }

}


// ----------------------------------------------------------------------
// gen-level analysis 
void plotLq::loopFunction4() {

  char cds[200];
  sprintf(cds, "%s", fCds.c_str());

  for (int i = 0; i < fRtd.ngen; ++i) {

    

  }
  

}



// ----------------------------------------------------------------------
// gen-level analysis for single vs pair production
void plotLq::loopFunction3() {

  char cds[200];
  sprintf(cds, "%s", fCds.c_str());

  TLorentzVector pL[2], pQ[2], pX[2], pXL[2], pXQ[2]; 
  bool single2(false);

  if (fPair) {
    single2 = true;
    // -- first the correct combination for the LQ, ordered according to lepton pT!
    if (fRtd.glpt[0] > fRtd.glpt[1]) {
      pXL[0].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pXL[1].SetPtEtaPhiM(fRtd.glpt[1], fRtd.gleta[1], fRtd.glphi[1], 0);
      pXQ[0].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      pXQ[1].SetPtEtaPhiM(fRtd.gjpt[1], fRtd.gjeta[1], fRtd.gjphi[1], 0);
    } else {
      pXL[1].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pXL[0].SetPtEtaPhiM(fRtd.glpt[1], fRtd.gleta[1], fRtd.glphi[1], 0);
      pXQ[1].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      pXQ[0].SetPtEtaPhiM(fRtd.gjpt[1], fRtd.gjeta[1], fRtd.gjphi[1], 0);
    }
    pX[0] = pXL[0] + pXQ[0];
    pX[1] = pXL[1] + pXQ[1];

    // -- now sorted into leading and subleading 
    if (fRtd.glpt[0] > fRtd.glpt[1]) {
      pL[0].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pL[1].SetPtEtaPhiM(fRtd.glpt[1], fRtd.gleta[1], fRtd.glphi[1], 0);
      //       pQ[0].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      //       pQ[1].SetPtEtaPhiM(fRtd.gjpt[1], fRtd.gjeta[1], fRtd.gjphi[1], 0);
    } else {
      pL[1].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pL[0].SetPtEtaPhiM(fRtd.glpt[1], fRtd.gleta[1], fRtd.glphi[1], 0);
      //       pQ[1].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      //       pQ[0].SetPtEtaPhiM(fRtd.gjpt[1], fRtd.gjeta[1], fRtd.gjphi[1], 0);
    }
    if (fRtd.gjpt[0] > fRtd.gjpt[1]) {
      pQ[0].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      pQ[1].SetPtEtaPhiM(fRtd.gjpt[1], fRtd.gjeta[1], fRtd.gjphi[1], 0);
    } else {
      pQ[1].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      pQ[0].SetPtEtaPhiM(fRtd.gjpt[1], fRtd.gjeta[1], fRtd.gjphi[1], 0);
    }
  } else {
    // -- first the correct combination for the LQ
    pXL[0].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
    pXL[1].SetPtEtaPhiM(fRtd.gkpt[0], fRtd.gketa[0], fRtd.gkphi[0], 0);
    pXQ[0].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
    if (fRtd.gipt[0] > 0.) {
      pXQ[1].SetPtEtaPhiM(fRtd.gipt[0], fRtd.gieta[0], fRtd.giphi[0], 0);
      single2 = true;
    }

    pX[0] = pXL[0] + pXQ[0];
    if (single2) pX[1] = pXL[1] + pXQ[1];

    // -- now the ordering for leading and subleading
    if (fRtd.glpt[0] > fRtd.gkpt[0]) {
      pL[0].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      //       pQ[0].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      pL[1].SetPtEtaPhiM(fRtd.gkpt[0], fRtd.gketa[0], fRtd.gkphi[0], 0);
    } else {
      cout << "we should not get here" << endl;
      pL[1].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pL[0].SetPtEtaPhiM(fRtd.gkpt[0], fRtd.gketa[0], fRtd.gkphi[0], 0);
      //       pQ[0].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
    }

    if (fRtd.gjpt[0] > fRtd.gipt[0]) {
      pQ[0].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      if (single2) pQ[1].SetPtEtaPhiM(fRtd.gipt[0], fRtd.gieta[0], fRtd.giphi[0], 0);
    } else {
      pQ[1].SetPtEtaPhiM(fRtd.gjpt[0], fRtd.gjeta[0], fRtd.gjphi[0], 0);
      if (single2) pQ[0].SetPtEtaPhiM(fRtd.gipt[0], fRtd.gieta[0], fRtd.giphi[0], 0);
    }      
  }
    
  fHists[Form("ngen_%s", cds)]->Fill(fRtd.ngen); 
  
  fHists[Form("ptl0_%s", cds)]->Fill(pL[0].Pt()); 
  fHists[Form("ptl1_%s", cds)]->Fill(pL[1].Pt()); 
  fHists[Form("ptq0_%s", cds)]->Fill(pQ[0].Pt()); 
  if (single2)   fHists[Form("ptq1_%s", cds)]->Fill(pQ[1].Pt()); 
  fHists[Form("ptlq_%s", cds)]->Fill(pX[0].Pt()); 
  if (single2) fHists[Form("ptlq_%s", cds)]->Fill(pX[1].Pt()); 
  
  fHists[Form("etaq0_%s", cds)]->Fill(pQ[0].Eta()); 
  if (single2) fHists[Form("etaq1_%s", cds)]->Fill(pQ[1].Eta()); 
  fHists[Form("etalq_%s", cds)]->Fill(pX[0].Eta()); 
  if (single2) fHists[Form("etalq_%s", cds)]->Fill(pXQ[1].Eta()); 
  
  fHists[Form("dRll_%s", cds)]->Fill(pL[0].DeltaR(pL[1])); 
  
  fHists[Form("dRxql_%s", cds)]->Fill(pXQ[0].DeltaR(pXL[0])); 
  if (single2) fHists[Form("dRxql_%s", cds)]->Fill(pXQ[1].DeltaR(pXL[1])); 
  fHists[Form("dRql_%s", cds)]->Fill(pQ[0].DeltaR(pL[0])); 
  if (single2) fHists[Form("dRql_%s", cds)]->Fill(pQ[1].DeltaR(pL[1])); 
  fHists[Form("dFxql_%s", cds)]->Fill(pXQ[0].DeltaPhi(pXL[0])); 
  if (single2) fHists[Form("dFxql_%s", cds)]->Fill(pXQ[1].DeltaPhi(pXL[1])); 
  fHists[Form("dFql_%s", cds)]->Fill(pQ[0].DeltaPhi(pL[0])); 
  if (single2) fHists[Form("dFql_%s", cds)]->Fill(pQ[1].DeltaPhi(pL[1])); 
  
  fHists[Form("mxql_%s", cds)]->Fill(pX[0].M()); 
  if (single2) fHists[Form("mxql_%s", cds)]->Fill(pX[1].M()); 

  fHists[Form("mqXlX_%s", cds)]->Fill((pQ[0]+pL[0]).M()); 
  if (single2) fHists[Form("mqXlX_%s", cds)]->Fill((pQ[1]+pL[1]).M()); 
  
  fHists[Form("mqXlY_%s", cds)]->Fill((pQ[0]+pL[1]).M()); 
  if (single2) fHists[Form("mqXlY_%s", cds)]->Fill((pQ[1]+pL[0]).M()); 
  
  fHists[Form("mll_%s", cds)]->Fill((pL[0]+pL[1]).M()); 
  
  fHists[Form("dPtll_%s", cds)]->Fill(pL[0].Pt() - pL[1].Pt()); 

  fHists[Form("sPtll_%s", cds)]->Fill(pL[0].Pt() + pL[1].Pt()); 
  fHists[Form("ptll_%s", cds)]->Fill((pL[0]+pL[1]).Pt()); 
  
  if (single2) {
    fHists[Form("dPtqq_%s", cds)]->Fill(pQ[0].Pt() - pQ[1].Pt()); 
    fHists[Form("sPtqq_%s", cds)]->Fill(pQ[0].Pt() + pQ[1].Pt()); 
    fHists[Form("ptqq_%s", cds)]->Fill((pQ[0]+pQ[1]).Pt()); 
    fHists[Form("st_%s", cds)]->Fill(pL[0].Pt() + pL[1].Pt() + pQ[0].Pt() + pQ[1].Pt()); 
    fHists[Form("dst_%s", cds)]->Fill((pXL[0].Pt() + pXQ[0].Pt()) - (pXL[1].Pt() + pXQ[1].Pt())); 
    fHists[Form("dpt_%s", cds)]->Fill((pXL[0] + pXQ[0]).Pt() - (pXL[1] + pXQ[1]).Pt()); 
  }
}


// ----------------------------------------------------------------------
void plotLq::loopFunction2() {

  // -- add possible cand selection cuts that are not optimized
  // .. HERE ..

  for (unsigned int i = 0; i < fSelPoints.size(); ++i) {
    fSelPoints[i].eval(fOptMode, 1.);
  }

}

// ----------------------------------------------------------------------
void plotLq::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries); 
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0; 
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart; 
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts; 
    } else {
      nend = t->GetEntries();
    }
  }
  
  nentries = nend - nstart; 
  
  int step(1000000); 
  if (nentries < 5000000)  step = 500000; 
  if (nentries < 1000000)  step = 100000; 
  if (nentries < 100000)   step = 10000; 
  if (nentries < 10000)    step = 1000; 
  if (nentries < 1000)     step = 100; 
  if (2 == ifunc)          step = 10000; 
  cout << "==> plotLq::loopOverTree> loop over dataset " << fCds << " in file " 
       << t->GetDirectory()->GetName() 
       << " with " << nentries << " entries"  << " looping from  " << nbegin << " .. " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotLq::*pF)(void);
  if (ifunc == 1) pF = &plotLq::loopFunction1;
  if (ifunc == 2) pF = &plotLq::loopFunction2;
  if (ifunc == 3) pF = &plotLq::loopFunction3;
  if (ifunc == 4) pF = &plotLq::loopFunction4;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    (this->*pF)();
  }

}

// ----------------------------------------------------------------------
void plotLq::setupTree(TTree *t) {

  t->SetBranchAddress("type",    &fRtd.type);
  t->SetBranchAddress("w8",      &fRtd.w8);
  
  t->SetBranchAddress("ngen",    &fRtd.ngen);
  t->SetBranchAddress("gm",       fRtd.gm);
  t->SetBranchAddress("gpt",      fRtd.gpt);
  t->SetBranchAddress("geta",     fRtd.geta);
  t->SetBranchAddress("gphi",     fRtd.gphi);
  t->SetBranchAddress("gmlj",     fRtd.gmlj);
  t->SetBranchAddress("glq",      fRtd.glq);
  t->SetBranchAddress("gtm",      fRtd.gtm);

  t->SetBranchAddress("glpt",     fRtd.glpt);
  t->SetBranchAddress("gleta",    fRtd.gleta);
  t->SetBranchAddress("glphi",    fRtd.glphi);

  t->SetBranchAddress("gqpt",     fRtd.gqpt);
  t->SetBranchAddress("gqeta",    fRtd.gqeta);
  t->SetBranchAddress("gqphi",    fRtd.gqphi);

  t->SetBranchAddress("gjpt",     fRtd.gjpt);
  t->SetBranchAddress("gjeta",    fRtd.gjeta);
  t->SetBranchAddress("gjphi",    fRtd.gjphi);

  t->SetBranchAddress("gipt",     fRtd.gipt);
  t->SetBranchAddress("gieta",    fRtd.gieta);
  t->SetBranchAddress("giphi",    fRtd.giphi);

  t->SetBranchAddress("gkpt",     fRtd.gkpt);
  t->SetBranchAddress("gketa",    fRtd.gketa);
  t->SetBranchAddress("gkphi",    fRtd.gkphi);

  t->SetBranchAddress("nrec",    &fRtd.nrec);
  t->SetBranchAddress("m",        fRtd.m);
  t->SetBranchAddress("pt",       fRtd.pt);
  t->SetBranchAddress("eta" ,     fRtd.eta);
  t->SetBranchAddress("phi",      fRtd.phi);
  t->SetBranchAddress("lq",        fRtd.lq);

  t->SetBranchAddress("lpt",      fRtd.lpt);
  t->SetBranchAddress("leta" ,    fRtd.leta);
  t->SetBranchAddress("lphi",     fRtd.lphi);

  t->SetBranchAddress("jpt",      fRtd.jpt);
  t->SetBranchAddress("jeta" ,    fRtd.jeta);
  t->SetBranchAddress("jphi",     fRtd.jphi);

  t->SetBranchAddress("kpt",      fRtd.kpt);
  t->SetBranchAddress("keta" ,    fRtd.keta);
  t->SetBranchAddress("kphi",     fRtd.kphi);

  t->SetBranchAddress("st",       fRtd.st);
  t->SetBranchAddress("mll" ,     fRtd.mll);
  t->SetBranchAddress("mljmin",   fRtd.mljmin);



}

// ----------------------------------------------------------------------
void plotLq::splitType(string stype, string &name, double &mass, double &lambda) {    
  istringstream ss(stype);
  string token;
  
  getline(ss, token, ',');
  name = token;
  cout << "name: " << name << endl;

  getline(ss, token, ','); 
  mass = atof(token.c_str());
  cout << "mass: " << mass << endl;

  getline(ss, token, ','); 
  lambda = atof(token.c_str());
  cout << "lambda: " << lambda << endl;
}

// ----------------------------------------------------------------------
void plotLq::loadFiles(string afiles) {
  
  string files = fDirectory + "/" + afiles;
  cout << "==> Loading files listed in " << files << endl;

  char buffer[1000];
  ifstream is(files.c_str());
  string sname, sdecay; 
  double mass, lambda; 

  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    if (buffer[0] == '\n') {continue;}
    
    string sbuffer = string(buffer); 
    replaceAll(sbuffer, "\t", " ");
    replaceAll(sbuffer, "  ", " ");

    string::size_type m1 = sbuffer.find("xsec="); 
    string stype = sbuffer.substr(5, m1-6); 
    splitType(stype, sname, mass, lambda); 
    
    cout << "sname: " << sname << endl;

    string::size_type m2 = sbuffer.find("file="); 
    string sxsec = sbuffer.substr(m1+5, m2-m1-6); 
    string sfile = sbuffer.substr(m2+5); 

    TFile *pF(0); 
    // -- MC
    pF = loadFile(sfile); 
    sdecay = "single";
    if (string::npos != sname.find("pair")) sdecay = "pair";
    TTree *t = (TTree*)pF->Get(Form("%s/events", sdecay.c_str())); 
    int nevt = t->GetEntries();
    if (string::npos != sname.find("dy")) {
      dataset *ds = new dataset(); 
      sdecay = "Drell-Yan";
      ds->fColor = kGreen+2; 
      ds->fLcolor = kGreen+2; 
      ds->fFcolor = kGreen+2; 
      ds->fSymbol = 25; 

      ds->fF      = pF; 
      ds->fXsec   = atof(sxsec.c_str());          // [xsec] = pb
      ds->fBf     = 1.;
      ds->fMass   = -1.;
      ds->fLambda = -1.;
      ds->fLumi   = nevt/ds->fXsec/ds->fBf/1000.; // [lumi] = 1/fb
      //      ds->fName   = "MadGraph " + sdecay; 
      ds->fName   = sdecay; 
      ds->fFillStyle = 3657; 
      ds->fSize = 0.5; 
      ds->fWidth = 1.0; 
      fDS.insert(make_pair(sname, ds)); 
    }

    if (string::npos != sname.find("ttbarjj")) {
      dataset *ds = new dataset(); 
      sdecay = "ttbar#rightarrow jj";
      ds->fColor = kMagenta+2; 
      ds->fLcolor = kMagenta+2; 
      ds->fFcolor = kMagenta+2; 
      ds->fSymbol = 26; 

      ds->fF      = pF; 
      ds->fXsec   = atof(sxsec.c_str());          // [xsec] = pb
      ds->fBf     = 1.;
      ds->fMass   = -1.;
      ds->fLambda = -1.;
      ds->fLumi   = nevt/ds->fXsec/ds->fBf/1000.; // [lumi] = 1/fb
      //      ds->fName   = "MadGraph " + sdecay; 
      ds->fName   = sdecay; 
      ds->fFillStyle = 3375; 
      ds->fSize = 0.5; 
      ds->fWidth = 1.0; 
      fDS.insert(make_pair(sname, ds)); 
    }


    if (string::npos != sname.find("lq")) {
      dataset *ds = new dataset(); 
      sdecay = "LQ";
      if (string::npos != sname.find("pair")) sdecay = "LQ #bar{LQ}";
      sdecay = Form("%s (%.0fGeV, #Lambda=%2.1f)", sdecay.c_str(), mass, lambda);
      if (string::npos != sdecay.find("bar")) {
	ds->fColor = kBlue; 
	ds->fLcolor = kBlue; 
	ds->fFcolor = kBlue; 
	ds->fSymbol = 24; 
	ds->fFillStyle = 3653; 
      } else {
	ds->fColor = kRed; 
	ds->fLcolor = kRed; 
	ds->fFcolor = kRed; 
	ds->fSymbol = 25; 
	ds->fFillStyle = 3635; 
      }
      ds->fF      = pF; 
      ds->fXsec   = atof(sxsec.c_str());          // [xsec] = pb
      ds->fBf     = 1.;
      ds->fMass   = mass;
      ds->fLambda = lambda;
      ds->fLumi   = nevt/ds->fXsec/ds->fBf/1000.; // [lumi] = 1/fb
      ds->fName   = sdecay; 
      ds->fSize = 0.5; 
      ds->fWidth = 1.0; 
      fDS.insert(make_pair(sname, ds)); 
    }

    // mb ub nb pb fb 
    cout << "opened MC file "  << sfile  << " as " << sname << " (" << stype << ") with xsec = " << sxsec
	 << Form(" = %8.5f", fDS[sname]->fXsec)
	 << Form(", equivalent lumi = %5.1f/fb", fDS[sname]->fLumi)
	 << endl;

  }

}


// ----------------------------------------------------------------------
void plotLq::genMass(string type, int offset, int nplot) {
  makeCanvas(1); 
  c1->Clear();

  newLegend(0.2, 0.3, 0.6, 0.8);

  string dir = "pair";
  if (string::npos != type.find("single")) {
    dir = "single"; 
  }
  string hname = Form("%s/m", dir.c_str());

  vector<string> ds; 
  for (int i = 0; i < nplot; ++i) {
    ds.push_back(Form("%s_%d", type.c_str(), offset+i)); 
  }
//   ds.push_back(Form("lq_%s_%d", type.c_str(), offset+i)); 
//   ds.push_back(Form("lq_%s_%d", type.c_str(), offset)); 
//   ds.push_back(Form("lq_%s_%d", type.c_str(), offset)); 

  
  TH1D *h(0); 
  for (unsigned int i = 0; i < ds.size(); ++i) {
    cout << ds[i] << " mass = " << fDS[ds[i]]->fMass << "GeV, lambda = " << fDS[ds[i]]->fLambda << endl;
    h = fDS[ds[i]]->getHist(hname);
    if (0 == i) setHist(h, kBlack); 
    if (1 == i) setHist(h, kBlue); 
    if (2 == i) setHist(h, kMagenta); 
    if (3 == i) setHist(h, kRed); 
    if (0 == i) {
      h->SetMaximum(2.*h->GetMaximum());
      h->DrawCopy();
    } else {
      h->DrawCopy("same");
    }
    legg->AddEntry(h, Form("m = %4.0f, #lambda = %3.2f", fDS[ds[i]]->fMass, fDS[ds[i]]->fLambda), "l");
  }

  legg->Draw();
	      
//   TH1D *h1 = fDS[Form("lq_%s_30", type.c_str())]->getHist(hname);
//   setHist(h1, kBlue); 
//   h1->Draw("same");
//   legg->Add(h0, Form("m = %4.0f, #Lambda = %3.2f", fDS[Form("lq_%s_30", type.c_str())]->fMass, fDS[Form("lq_%s_30", type.c_str())]->fLambda)

//   TH1D *h2 = fDS[Form("lq_%s_31", type.c_str())]->getHist(hname);
//   setHist(h2, kMagenta); 
//   h2->Draw("same");

//   TH1D *h3 = fDS[Form("lq_%s_32", type.c_str())]->getHist(hname);
//   setHist(h3, kRed); 
//   h3->DrawCopy("same");

	    
}


// ----------------------------------------------------------------------
void plotLq::genMassPlots(string dir) {

  map<string, dataset*>::iterator ids = fDS.begin();
  map<string, dataset*>::iterator eds = fDS.end();
  TH1D *ha = new TH1D("ha", "", 200, 0., 2000.); setHist(ha, kBlack); 
  setTitles(ha, "m [GeV]", "Entries/bin", 0.05, 1.1, 2.0); 
  TH1D *hr = new TH1D("hr", "", 200, 0., 2000.); setFilledHist(hr, kBlue, kBlue, 3653); 
  TH1D *hn = new TH1D("hn", "", 200, 0., 2000.); setFilledHist(hn, kRed, kRed, 3635); 
  double na(0.), nr(0.), nn(0.); 
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.03);
  for (; ids != eds; ++ids) {
    if (string::npos == ids->first.find("lq")) continue;
    if (string::npos == ids->first.find(dir)) continue;
    if (string::npos == ids->first.find("down") && string::npos == ids->first.find("up")) continue;
    ha->Reset(); 
    ha->SetTitle(Form("%s (m = %4.0f GeV,    #lambda = %3.2f)", 
		      ids->first.c_str(), ids->second->fMass, ids->second->fLambda)); 
    hr->Reset(); 
    hn->Reset(); 
    cout << ids->first << endl;
    TTree *t = getTree(ids->first, dir);
    t->Draw("gm>>ha", "");
    t->Draw("gm>>hr", "gtm==0");
    int oflow = t->Draw("gm>>hn", "gtm>0&&gm>2000");
    t->Draw("gm>>hn", "gtm>0");
    na = ha->Integral();
    nr = hr->Integral();
    nn = hn->Integral();
    cout << "all: " << na
	 <<  " resonant: " << nr
	 <<  " non-resonant: " << nn
	 << " overflow: " << oflow
	 << endl;
    zone();
    shrinkPad(0.15, 0.2); 
    gPad->SetLogy(1);
    ha->SetMinimum(0.5);
    double hamax = hr->GetMaximum(); 
    if (hn->GetMaximum() > hamax) hamax = hn->GetMaximum();
    ha->SetMaximum(10.*hamax); 
    ha->Draw("hist");
    hr->Draw("histsame");
    hn->Draw("histsame");

    tl->DrawLatex(0.25, 0.86, Form("%5.0f (all)", na)); 
    tl->DrawLatex(0.25, 0.83, Form("%5.0f (resonant)", nr)); 
    tl->DrawLatex(0.25, 0.80, Form("%5.0f (non-resonant)", nn)); 

    c0->SaveAs(Form("%s/genmass-%s-%s.pdf", fDirectory.c_str(), dir.c_str(), ids->first.c_str())); 
  }


}


// ----------------------------------------------------------------------
void plotLq::showVar(string var, string cuts, string signal, string dir) {

  gStyle->SetOptStat(0); 

  int NBINS(100); 
  double xlo(0.), xhi(2000.);
  string xtitle("m(LQ) [GeV]");
  if (var == "pt") {
    NBINS = 100;
    xlo = 0.;
    xhi = 1000.;
  }

  // -- signal
  TH1D *hS = new TH1D("hS", "", NBINS, xlo, xhi); setFilledHist(hS, kBlue, kBlue, 3653); 
  // -- cross feed
  TH1D *hX = new TH1D("hX", "", NBINS, xlo, xhi); setFilledHist(hX, kRed, kRed, 3635); 
  // -- DY
  TH1D *hD = new TH1D("hD", "", NBINS, xlo, xhi); setFilledHist(hD, kRed, kRed, 3635); 
  setHist(hD, fDS["dy"]); 
  // -- ttbar
  TH1D *hT = new TH1D("hT", "", NBINS, xlo, xhi); setFilledHist(hT, kRed, kRed, 3635); 
  setHist(hT, fDS["ttbarjj"]); 

  TTree *t = getTree(signal, dir);
  //  string allCuts = (cuts == ""?string("tm>-1") : string(Form("tm>-1&&%s", cuts.c_str())));
  string allCuts = cuts; 
  t->Draw("m>>hS", allCuts.c_str(), "goff");
  normHist(hS, signal, LUMI);
  setTitles(hS, xtitle.c_str(), hS->GetYaxis()->GetTitle(), 0.05, 1.0, 1.5);
  hS->Draw("hist");

  string xfsignal = xfeed(signal); 
  cout << "===> signal = " << signal << " -> xf = " << xfsignal << endl;
  
  double ymax = hS->GetMaximum(); 

  shrinkPad(0.1, 0.15); 
  t = getTree(xfsignal, dir);
  t->Draw("m>>hX", allCuts.c_str(), "goff");
  normHist(hX, signal, LUMI);
  hX->Draw("samehist");
  if (hX->GetMaximum() > ymax) ymax = hX->GetMaximum();

  t = getTree("ttbarjj", dir);
  t->Draw("m>>hT", allCuts.c_str(), "goff");
  normHist(hT, "ttbarjj", LUMI);
  hT->Draw("samehist");
  if (hT->GetMaximum() > ymax) ymax = hT->GetMaximum();

  t = getTree("dy", dir);
  t->Draw("m>>hD", allCuts.c_str(), "goff");
  normHist(hT, "dy", LUMI);
  hD->Draw("samehist");
  if (hD->GetMaximum() > ymax) ymax = hD->GetMaximum();

  hS->SetMaximum(1.2*ymax); 
  
  newLegend(0.6, 0.7, 0.85, 0.85);
  legg->SetTextSize(0.03);
  legg->SetHeader(dir.c_str()); 
  legg->AddEntry(hS, Form("%s", signal.c_str()), "f");
  legg->AddEntry(hX, Form("%s", xfsignal.c_str()), "f");
  legg->AddEntry(hT, Form("%s", "ttbarjj_a"), "f");
  legg->AddEntry(hD, Form("%s", "Drell-Yan"), "f");
  
  legg->Draw();

}


// ----------------------------------------------------------------------
string plotLq::xfeed(string signal) {
  if (string::npos != signal.find("single")) {
    replaceAll(signal, "single", "pair"); 
  } else {
    replaceAll(signal, "pair", "single"); 
  }
  return signal;
}
