#include "plotLq.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"

#include "dataset.hh"
#include "util.hh"

ClassImp(plotLq)

using namespace std; 

// ----------------------------------------------------------------------
plotLq::plotLq(string dir,  string files, string setup) {

  fDirectory = dir; 

  loadFiles(files);

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  fEpsilon = 0.00001; 

  NBINS = 50; 

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);

}


// ----------------------------------------------------------------------
plotLq::~plotLq() {

}


// ----------------------------------------------------------------------
void plotLq::bookHist(string name) {
  fHists.insert(make_pair(Form("m_%s", name.c_str()), 
			  new TH1D(Form("m_%s", name.c_str()), Form("m_%s", name.c_str()), NBINS, 0, 2000.))); 

  fHists.insert(make_pair(Form("pt_%s", name.c_str()), 
			  new TH1D(Form("pt_%s", name.c_str()), Form("pt_%s", name.c_str()), 100, 0, 1000.))); 

}


// ----------------------------------------------------------------------
void plotLq::makeAll(int bitmask) {
  if (bitmask & 0x1) treeAnalysis();
}


// ----------------------------------------------------------------------
void plotLq::treeAnalysis() {

  string ds("lq_pair_01");
  fPair = true;
  fCds = ds; 
  bookHist(ds); 
  TTree *t = getTree(ds); 
  setupTree(t); 
  loopOverTree(t); 

  ds = "lq_pair_02"; 
  fPair = true;
  fCds = ds; 
  bookHist(ds); 
  t = getTree(ds); 
  setupTree(t);
  loopOverTree(t); 


  fHists["m_lq_pair_01"]->Draw(); 
  fHists["m_lq_pair_02"]->Draw("same"); 

}


// ----------------------------------------------------------------------
void plotLq::normHist(TH1D *h, double integral, string type) {
  double scale(1.); 
  if (TMath::Abs(integral - 1.) < fEpsilon) {
    // -- normalize to 1
    scale = (h->GetSumOfWeights() > 0 ? integral/h->GetSumOfWeights() : 1.); 
  } else if (TMath::Abs(integral + 1.) < fEpsilon) {
    // -- normalize to xsec*bf
    //    n = xsec * L
    //    "integral" over histogram should be xsec

    scale = (h->GetSumOfWeights() > 0 ? fDS[type]->fXsec*fDS[type]->fBf/h->Integral() : 1.); 
    setTitles(h, h->GetXaxis()->GetTitle(), "pb");
  } else {
    scale = 1.;
  }
  double c(0.), e(0.); 
  for (int i = 0; i <= h->GetNbinsX(); ++i) {
    c = h->GetBinContent(i); 
    e = h->GetBinError(i); 
    h->SetBinContent(i, c*scale);
    h->SetBinError(i, e*scale);
  }
  
}


// ----------------------------------------------------------------------
void plotLq::overlayAll() {

  // -- simple overlays
  c0->cd(1); 
  overlay("lq_pair_01", "bla", "lq_pair_02", "bla"); 

}

// ----------------------------------------------------------------------
void plotLq::overlay(TH1D* h1, string f1, TH1D* h2, string f2, bool legend) {

  bool log(false); 

  normHist(h1, -1., f1); 
  normHist(h2, -1., f2); 

  double hmax(h1->GetMaximum()); 
  if (h2->GetMaximum() > hmax) hmax = 1.2*h2->GetMaximum(); 
  if (log) {
    gPad->SetLogy(1); 
    hmax *= 2.;
  }
  h1->SetMaximum(hmax); 
  h1->DrawCopy("e"); 
  h2->DrawCopy("histsame");
  cout << "overlay(" << f1 << ", " << h1->GetName() << ", " << f2 << ", " << h2->GetName() 
       << ") legend = " << legend << " log: " << log 
       << endl;
  
  if (legend) {
    newLegend(0.40, 0.75, 0.7, 0.85); 
    legg->AddEntry(h1, fDS[f1]->fName.c_str(), "p"); 
    legg->AddEntry(h2, fDS[f2]->fName.c_str(), "l"); 
    legg->Draw();
    cout << "  drawing legend" << endl;
  }


}


// ----------------------------------------------------------------------
void plotLq::overlay(string f1, string h1name, string f2, string h2name, bool legend) {

  TH1D *h1 = fDS[f1]->getHist(Form("%s", h1name.c_str())); 
  TH1D *h2 = fDS[f2]->getHist(Form("%s", h2name.c_str())); 

  overlay(h1, f1, h2, f2, legend); 

}


// ----------------------------------------------------------------------
void plotLq::candAnalysis() {

  fGoodEvent   = false; 
  
  fGoodCandLQp = false; 
  fGoodCandLQn = false; 

  if (fRtd.ljnm > 0.) fGoodCandLQn = true; 
  if (fRtd.ljpm > 0.) fGoodCandLQp = true; 

  if (fPair) {
    if (fGoodCandLQn && fGoodCandLQp) fGoodEvent = true;
  } else {
    if (fGoodCandLQn || fGoodCandLQp) fGoodEvent = true;
  }

}

// ----------------------------------------------------------------------
void plotLq::loopFunction() {
  char cds[200];
  sprintf(cds, "%s", fCds.c_str());

  //  cout << "goodEvent: " << fGoodEvent << " m = " << fRtd.ljnm << " and " << fRtd.ljpm << endl;

  if (fGoodEvent) { 
    if (fRtd.ljnm > 0.) {
      fHists[Form("m_%s", cds)]->Fill(fRtd.ljnm); 
      fHists[Form("pt_%s", cds)]->Fill(fRtd.ljnpt); 
    }
    if (fRtd.ljpm > 0.) {
      fHists[Form("m_%s", cds)]->Fill(fRtd.ljpm); 
      fHists[Form("pt_%s", cds)]->Fill(fRtd.ljppt); 
    }
  }
}

// ----------------------------------------------------------------------
void plotLq::loopOverTree(TTree *t, int nevts, int nstart) {
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
  cout << "==> plotLq::loopOverTree> loop over dataset " << fCds << " in file " 
       << t->GetDirectory()->GetName() 
       << " with " << nentries << " entries"  << " looping from  " << nbegin << " .. " << nend
       << endl;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    
    candAnalysis();
    loopFunction();
  }

}


// ----------------------------------------------------------------------
void plotLq::setupTree(TTree *t) {

  t->SetBranchAddress("type",    &fRtd.type);
  t->SetBranchAddress("w8",      &fRtd.w8);

  t->SetBranchAddress("gpm",     &fRtd.gpm);
  t->SetBranchAddress("gpm2",    &fRtd.gpm2);
  t->SetBranchAddress("gppt",    &fRtd.gppt);
  t->SetBranchAddress("gpeta",   &fRtd.gpeta);

  t->SetBranchAddress("gnm",     &fRtd.gnm);
  t->SetBranchAddress("gnm2",    &fRtd.gnm2);
  t->SetBranchAddress("gnpt",    &fRtd.gnpt);
  t->SetBranchAddress("gneta",   &fRtd.gneta);

  t->SetBranchAddress("glqpm",   &fRtd.glqpm);
  t->SetBranchAddress("gljpm",   &fRtd.gljpm);

  t->SetBranchAddress("glqnm",   &fRtd.glqnm);
  t->SetBranchAddress("gljnm",   &fRtd.gljnm);

  
  t->SetBranchAddress("ljnm",    &fRtd.ljnm);
  t->SetBranchAddress("ljnpt",   &fRtd.ljnpt);
  t->SetBranchAddress("ljneta",  &fRtd.ljneta);

  t->SetBranchAddress("ljpm",    &fRtd.ljpm);
  t->SetBranchAddress("ljppt",   &fRtd.ljppt);
  t->SetBranchAddress("ljpeta",  &fRtd.ljpeta);

  t->SetBranchAddress("st",      &fRtd.st);
  t->SetBranchAddress("mll",     &fRtd.mll);
  t->SetBranchAddress("mljetmin",&fRtd.mljetmin);

}


// ----------------------------------------------------------------------
TTree* plotLq::getTree(string ds) {
  TTree *t(0);
  t = (TTree*)fDS[ds]->fF->Get("events"); 
  return t; 
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
    
    string sbuffer = string(buffer); 

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
    TTree *t = (TTree*)pF->Get("events"); 
    int nevt = t->GetEntries();
    if (string::npos != sname.find("dy")) {
      dataset *ds = new dataset(); 
      sdecay = "Drell-Yan";
      ds->fColor = kRed; 
      ds->fLcolor = kRed; 
      ds->fFcolor = kRed; 
      ds->fSymbol = 25; 

      ds->fF      = pF; 
      ds->fXsec   = atof(sxsec.c_str());
      ds->fBf     = 1.;
      ds->fMass   = -1.;
      ds->fLambda = -1.;
      ds->fLumi   = nevt/ds->fXsec/ds->fBf;
      ds->fName   = "MadGraph " + sdecay; 
      ds->fFillstyle = 3365; 
      ds->fSize = 1; 
      ds->fWidth = 2; 
      fDS.insert(make_pair(sname, ds)); 
      cout << "  inserted into fDS" << endl;
    }


    if (string::npos != sname.find("lq")) {
      dataset *ds = new dataset(); 
      sdecay = "LQ";
      if (string::npos != sname.find("pair")) sdecay = "LQ LQ";
      ds->fColor = kBlue; 
      ds->fLcolor = kBlue; 
      ds->fFcolor = kBlue; 
      ds->fSymbol = 24; 

      ds->fF      = pF; 
      ds->fXsec   = atof(sxsec.c_str());
      ds->fBf     = 1.;
      ds->fMass   = mass;
      ds->fLambda = lambda;
      ds->fLumi   = nevt/ds->fXsec/ds->fBf;
      ds->fName   = "MadGraph " + sdecay; 
      ds->fFillstyle = 3356; 
      ds->fSize = 1; 
      ds->fWidth = 2; 
      fDS.insert(make_pair(sname, ds)); 
    }



    // mb ub nb pb fb 
    cout << "opened MC file "  << sfile  << " as " << sname << " (" << stype << ") with xsec = " << sxsec
	 << " = " << fDS[sname]->fXsec
	 << ", equivalent lumi = " << fDS[sname]->fLumi/1000. << "/fb"
	 << endl;
    
  }

}


// ----------------------------------------------------------------------
TFile* plotLq::loadFile(string file) {
  TFile *f = TFile::Open(file.c_str());
  return f; 
}



// ----------------------------------------------------------------------
void plotLq::replaceAll(string &sInput, const string &oldString, const string &newString) {
  string::size_type foundpos = sInput.find(oldString);
  while (foundpos != string::npos)  {
    sInput.replace(sInput.begin() + foundpos, sInput.begin() + foundpos + oldString.length(), newString);
    foundpos = sInput.find(oldString);
  }
}

// ----------------------------------------------------------------------
void plotLq::newLegend(double x1, double y1, double x2, double y2, string title) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
}

// ----------------------------------------------------------------------
void plotLq::makeCanvas(int i) {
  if (i & 16) { 
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 1000);
    c5->ToggleEventStatus();
  }
  if (i & 8) { 
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    //    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1 = new TCanvas("c1", "c1", 20,  60, 1000, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) { 
    c2 = new TCanvas("c2", "c2", 300, 200, 400, 800);
    c2->ToggleEventStatus();
  }
}
