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

#include "util.hh"

ClassImp(plotLq)

using namespace std; 

// ----------------------------------------------------------------------
plotLq::plotLq(string dir,  string files, string setup): plotClass(dir, files, setup) {

  loadFiles(files);

  fHistFile = TFile::Open(Form("%s/plotLq.root", dir.c_str()), "RECREATE"); 
}


// ----------------------------------------------------------------------
plotLq::~plotLq() {
}

// ----------------------------------------------------------------------
void plotLq::bookHist(string name) {
  // -- m
  fHists.insert(make_pair(Form("m_%s", name.c_str()), 
			  new TH1D(Form("m_%s", name.c_str()), Form("m_%s", name.c_str()), 100, 0, 2000.))); 
  setTitles(fHists[Form("m_%s", name.c_str())], "m [GeV]", "Entries/bin");
  setHist(fHists[Form("m_%s", name.c_str())], fDS[name]);

  // -- st
  fHists.insert(make_pair(Form("st_%s", name.c_str()), 
			  new TH1D(Form("st_%s", name.c_str()), Form("st_%s", name.c_str()), 35, 0, 3500.))); 
  setTitles(fHists[Form("st_%s", name.c_str())], "S_{T} [GeV]", "Entries/bin");
  setHist(fHists[Form("st_%s", name.c_str())], fDS[name]);

  // -- mll
  fHists.insert(make_pair(Form("mll_%s", name.c_str()), 
			  new TH1D(Form("mll_%s", name.c_str()), Form("mll_%s", name.c_str()), 50, 0, 1500.))); 
  setTitles(fHists[Form("mll_%s", name.c_str())], "m_{l l} [GeV]", "Entries/bin");
  setHist(fHists[Form("mll_%s", name.c_str())], fDS[name]);

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
  if (bitmask & 0x1) treeAnalysis();
}

// ----------------------------------------------------------------------
void plotLq::treeAnalysis() {

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
void plotLq::optimizePairCuts(string sg, string bg, double lumi) {
  fPair = true;

  // -- signal
  fCds = sg; 
  bookHist(sg); 
  fOptMode = 1; 
  double sgScale = lumi/fDS[sg]->fLumi;
  TTree *ts = getTree(sg); 
  setupTree(ts); 
  loopOverTree(ts, 2, 1000); 

  // -- background
  fCds = bg; 
  bookHist(bg); 
  fOptMode = 2; 
  double bgScale = lumi/fDS[bg]->fLumi;
  TTree *tb = getTree(bg); 
  setupTree(tb); 
  loopOverTree(tb, 2, 1000); 


  fHistFile->cd();
  TTree *t = new TTree("opt", "opt");
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
  
  for (unsigned int i = 0; i < fSelPoints.size(); ++i) {
    s0  = fSelPoints[i].fSgCnt;
    b0  = fSelPoints[i].fBgCnt;
    s   = sgScale * fSelPoints[i].fSgCnt;
    b   = bgScale * fSelPoints[i].fBgCnt;
    ssb = (s+b>0? s/TMath::Sqrt(s+b):0.);
    sb  = (b>0?   s/TMath::Sqrt(b)  :0.);
    st  = fSelPoints[i].fLargerThan[0].second;
    mll = fSelPoints[i].fLargerThan[1].second;
    mlj = fSelPoints[i].fLargerThan[2].second;
    t->Fill();
  }
  
}


// ----------------------------------------------------------------------
void plotLq::loopFunction1() {

  // -- cuts
  fGoodEvent   = false; 
  
  if (fPair) {
    if (fRtd.m[0] > 0. && fRtd.m[1] > 0.) fGoodEvent = true;
    if (fRtd.st < 685.)       fGoodEvent = false;
    if (fRtd.mll < 150.)      fGoodEvent = false;
    if (fRtd.mljmin < 155.) fGoodEvent = false;
  } else {
    if (fRtd.m[0] > 0.) fGoodEvent = true;
  }

  char cds[200];
  sprintf(cds, "%s", fCds.c_str());

  if (fGoodEvent) { 
    fHists[Form("st_%s", cds)]->Fill(fRtd.st); 
    fHists[Form("mll_%s", cds)]->Fill(fRtd.mll); 
    fHists[Form("mljmin_%s", cds)]->Fill(fRtd.mljmin); 
      
    for (int i = 0; i < fRtd.nrec; ++i) {
      fHists[Form("m_%s", cds)]->Fill(fRtd.m[i]); 
      fHists[Form("pt_%s", cds)]->Fill(fRtd.pt[i]); 
    }
  }

}


// ----------------------------------------------------------------------
void plotLq::loopFunction2() {
  static bool first(true); 
  if (first) {
    first = false; 

    static const double stArr[] = { 300.,  320.,  340.,  360.,  380.,  400.,  420.,  440.,  460.,  480.,
				    500.,  520.,  540.,  560.,  580.,  600.,  620.,  640.,  660.,  680.,
				    700.,  750.,  800.,  850.,  900.,  950., 1000., 1050., 1100., 1150.,
                                   1200., 1250., 1300., 1350., 1400., 1450., 1500.
    };
    vector<double> stCuts(stArr, stArr + sizeof(stArr)/sizeof(stArr[0]));

    static const double mllArr[] = {100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 
				    200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 
				    300.
    };
    vector<double> mllCuts(mllArr, mllArr + sizeof(mllArr)/sizeof(mllArr[0]));


    static const double mljArr[] = {100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 
				    200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 
				    400., 450., 500., 550., 600., 650., 700., 750.
    };
    vector<double> mljCuts(mljArr, mljArr + sizeof(mljArr)/sizeof(mljArr[0]));

    for (unsigned int i = 0; i < stCuts.size(); ++i) {
      for (unsigned int j = 0; j < mljCuts.size(); ++j) {
	for (unsigned int k = 0; k < mllCuts.size(); ++k) {
	  selpoint s; 
	  s.fLargerThan.push_back(make_pair(&fRtd.st, stCuts[i])); 
	  s.fLargerThan.push_back(make_pair(&fRtd.mll, mllCuts[k])); 
	  s.fLargerThan.push_back(make_pair(&fRtd.mljmin, mljCuts[j])); 
	  fSelPoints.push_back(s); 
	}
      }
    }
  }

  for (unsigned int i = 0; i < fSelPoints.size(); ++i) {
    if (1 == fOptMode) {
      fSelPoints[i].eval(true);
    } else {
      fSelPoints[i].eval(false);
    }
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

  t->SetBranchAddress("glpt",     fRtd.glpt);
  t->SetBranchAddress("gleta",    fRtd.gleta);
  t->SetBranchAddress("glphi",    fRtd.glphi);

  t->SetBranchAddress("gqpt",     fRtd.gqpt);
  t->SetBranchAddress("gqeta",    fRtd.gqeta);
  t->SetBranchAddress("gqphi",    fRtd.gqphi);

  t->SetBranchAddress("gjpt",     fRtd.gjpt);
  t->SetBranchAddress("gjeta",    fRtd.gjeta);
  t->SetBranchAddress("gjphi",    fRtd.gjphi);

  t->SetBranchAddress("gkpt",     fRtd.gkpt);
  t->SetBranchAddress("gketa",    fRtd.gketa);
  t->SetBranchAddress("gkphi",    fRtd.gkphi);

  t->SetBranchAddress("nrec",    &fRtd.nrec);
  t->SetBranchAddress("m",        fRtd.m);
  t->SetBranchAddress("pt",       fRtd.pt);
  t->SetBranchAddress("eta" ,     fRtd.eta);
  t->SetBranchAddress("phi",      fRtd.phi);
  t->SetBranchAddress("q",        fRtd.lq);

  t->SetBranchAddress("lpt",      fRtd.lpt);
  t->SetBranchAddress("leta" ,    fRtd.leta);
  t->SetBranchAddress("lphi",     fRtd.lphi);

  t->SetBranchAddress("jpt",      fRtd.jpt);
  t->SetBranchAddress("jeta" ,    fRtd.jeta);
  t->SetBranchAddress("jphi",     fRtd.jphi);

  t->SetBranchAddress("kpt",      fRtd.kpt);
  t->SetBranchAddress("keta" ,    fRtd.keta);
  t->SetBranchAddress("kphi",     fRtd.kphi);

  t->SetBranchAddress("st",       &fRtd.st);
  t->SetBranchAddress("mll" ,     &fRtd.mll);
  t->SetBranchAddress("mljmin",   &fRtd.mljmin);



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
      ds->fXsec   = atof(sxsec.c_str());          // [xsec] = pb
      ds->fBf     = 1.;
      ds->fMass   = -1.;
      ds->fLambda = -1.;
      ds->fLumi   = nevt/ds->fXsec/ds->fBf/1000.; // [lumi] = 1/fb
      //      ds->fName   = "MadGraph " + sdecay; 
      ds->fName   = sdecay; 
      ds->fFillStyle = 3365; 
      ds->fSize = 1; 
      ds->fWidth = 2; 
      fDS.insert(make_pair(sname, ds)); 
      cout << "  inserted into fDS" << endl;
    }


    if (string::npos != sname.find("lq")) {
      dataset *ds = new dataset(); 
      sdecay = "LQ";
      if (string::npos != sname.find("pair")) sdecay = "LQ #bar{LQ}";
      sdecay = Form("%s (%.0fGeV, #Lambda=%2.1f)", sdecay.c_str(), mass, lambda);
      ds->fColor = kBlue; 
      ds->fLcolor = kBlue; 
      ds->fFcolor = kBlue; 
      ds->fSymbol = 24; 

      ds->fF      = pF; 
      ds->fXsec   = atof(sxsec.c_str());          // [xsec] = pb
      ds->fBf     = 1.;
      ds->fMass   = mass;
      ds->fLambda = lambda;
      ds->fLumi   = nevt/ds->fXsec/ds->fBf/1000.; // [lumi] = 1/fb
      ds->fName   = sdecay; 
      ds->fFillStyle = 3356; 
      ds->fSize = 1; 
      ds->fWidth = 2; 
      fDS.insert(make_pair(sname, ds)); 
    }

    // mb ub nb pb fb 
    cout << "opened MC file "  << sfile  << " as " << sname << " (" << stype << ") with xsec = " << sxsec
	 << Form(" = %8.5f", fDS[sname]->fXsec)
	 << Form(", equivalent lumi = %5.1f/fb", fDS[sname]->fLumi)
	 << endl;

  }

}


