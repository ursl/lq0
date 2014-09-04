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

  fHistFile = TFile::Open(Form("%s/plotLq.root", dir.c_str()), "RECREATE"); 
}


// ----------------------------------------------------------------------
plotLq::~plotLq() {
}

// ----------------------------------------------------------------------
void plotLq::bookHist(string name) {

  char cds[200];
  sprintf(cds, "%s", fCds.c_str());
  
  // -- start gen level
  char hist[200];
  sprintf(hist, "%s_%s", "ptl0", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 800.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(l0)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptl1", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 600.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(l1)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptq0", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 600.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(q0)", "Entries/bin");

  sprintf(hist, "%s_%s", "ptlq", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 600.))); 
  setHistTitles(fHists[hist], fDS[name], "p_{T}(lq)", "Entries/bin");


  sprintf(hist, "%s_%s", "etal0", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3., 3.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta_{T}(l0)", "Entries/bin");

  sprintf(hist, "%s_%s", "etal1", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3., 3.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta_{T}(l1)", "Entries/bin");

  sprintf(hist, "%s_%s", "etaq0", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3., 3.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta_{T}(q0)", "Entries/bin");

  sprintf(hist, "%s_%s", "etalq", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3., 3.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta_{T}(lq)", "Entries/bin");


  sprintf(hist, "%s_%s", "dRll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, 0, 6.))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta R(l,l)", "Entries/bin");

  sprintf(hist, "%s_%s", "dFll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3.15, 3.15))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta #phi(l,l)", "Entries/bin");


  sprintf(hist, "%s_%s", "dRql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, 0, 6.))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta R(q,l)", "Entries/bin");

  sprintf(hist, "%s_%s", "dFql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3.15, 3.15))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta #phi(q,l)", "Entries/bin");

  sprintf(hist, "%s_%s", "dRqk", fCds.c_str());
  fHists.insert(make_pair(hist,  new TH1D(hist, hist, 60, 0, 6.))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta R(q,k)", "Entries/bin");

  sprintf(hist, "%s_%s", "dFqk", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -3.15, 3.15))); 
  setHistTitles(fHists[hist], fDS[name], "#Delta #phi(q,k)", "Entries/bin");


  sprintf(hist, "%s_%s", "dEtall", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 60, -6., 6.))); 
  setHistTitles(fHists[hist], fDS[name], "#eta(l0) - #eta(l1)", "Entries/bin");

  sprintf(hist, "%s_%s", "dPtll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0., 500.))); 
  setHistTitles(fHists[hist], fDS[name], "pT(l0) - pT(l1) [GeV]", "Entries/bin");



  sprintf(hist, "%s_%s", "mql", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 1000.))); 
  setHistTitles(fHists[hist], fDS[name], "m(q,l) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "mqk", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 1000.))); 
  setHistTitles(fHists[hist], fDS[name], "m(q,k) [GeV]", "Entries/bin");

  sprintf(hist, "%s_%s", "mll", fCds.c_str());
  fHists.insert(make_pair(hist, new TH1D(hist, hist, 100, 0, 1000.))); 
  setHistTitles(fHists[hist], fDS[name], "m(l,k) [GeV]", "Entries/bin");
  // -- end gen level


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
  // -- pair
  fPair = true;
  fCds = "lq01_pair"; 
  bookHist(fCds); 
  TTree *tp = getTree(fCds); 
  setupTree(tp); 
  loopOverTree(tp, 3, 20000); 

  // -- single
  fPair = false;
  fCds = "lq01_single"; 
  bookHist(fCds); 
  TTree *ts = getTree(fCds); 
  setupTree(ts); 
  loopOverTree(ts, 3, 20000); 
  
  char hist[200];
  int ipad(1); 
  zone(4,4);
  sprintf(hist, "ptl0");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptl1");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dPtll");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptq0");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "ptlq");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "etal0");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "etal1");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "etaq0");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "etalq");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dFll");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dFql");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dFqk");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dRll");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dEtall");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);


  c0->cd(++ipad); 
  sprintf(hist, "dRql");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);

  c0->cd(++ipad); 
  sprintf(hist, "dRqk");
  overlay(fHists[Form("%s_%s", hist, "lq01_pair")], "lq01_pair", fHists[Form("%s_%s", hist, "lq01_single")], "lq01_single", UNITY, false, false);


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
// gen-level analysis for single vs pair production
void plotLq::loopFunction3() {

  char cds[200];
  sprintf(cds, "%s", fCds.c_str());
  
  if (fPair) {

    TLorentzVector pL[2], pQ[2]; 
    if (fRtd.glpt[0] > fRtd.glpt[1]) {
      pL[0].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pL[1].SetPtEtaPhiM(fRtd.glpt[1], fRtd.gleta[1], fRtd.glphi[1], 0);

      pQ[0].SetPtEtaPhiM(fRtd.gqpt[0], fRtd.gqeta[0], fRtd.gqphi[0], 0);
      pQ[1].SetPtEtaPhiM(fRtd.gqpt[1], fRtd.gqeta[1], fRtd.gqphi[1], 0);
    } else {
      pL[1].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pL[0].SetPtEtaPhiM(fRtd.glpt[1], fRtd.gleta[1], fRtd.glphi[1], 0);

      pQ[1].SetPtEtaPhiM(fRtd.gqpt[0], fRtd.gqeta[0], fRtd.gqphi[0], 0);
      pQ[0].SetPtEtaPhiM(fRtd.gqpt[1], fRtd.gqeta[1], fRtd.gqphi[1], 0);
    }
    
    fHists[Form("ptl0_%s", cds)]->Fill(pL[0].Pt()); 
    fHists[Form("ptl1_%s", cds)]->Fill(pL[1].Pt()); 
    fHists[Form("ptq0_%s", cds)]->Fill(pQ[0].Pt()); 
    fHists[Form("ptlq_%s", cds)]->Fill((pQ[0] + pL[0]).Pt()); 
    fHists[Form("ptlq_%s", cds)]->Fill((pQ[1] + pL[1]).Pt()); 
    
    fHists[Form("etal0_%s", cds)]->Fill(pL[0].Eta()); 
    fHists[Form("etal1_%s", cds)]->Fill(pL[1].Eta()); 
    fHists[Form("etaq0_%s", cds)]->Fill(pQ[0].Eta()); 
    fHists[Form("etalq_%s", cds)]->Fill((pQ[0] + pL[0]).Eta()); 
    fHists[Form("etalq_%s", cds)]->Fill((pQ[1] + pL[1]).Eta()); 

    fHists[Form("dFll_%s", cds)]->Fill(pL[0].DeltaPhi(pL[1])); 
    fHists[Form("dRll_%s", cds)]->Fill(pL[0].DeltaR(pL[1])); 

    fHists[Form("dRql_%s", cds)]->Fill(pQ[0].DeltaR(pL[0])); 
    fHists[Form("dRql_%s", cds)]->Fill(pQ[1].DeltaR(pL[1])); 
    fHists[Form("dFql_%s", cds)]->Fill(pQ[0].DeltaPhi(pL[0])); 
    fHists[Form("dFql_%s", cds)]->Fill(pQ[1].DeltaPhi(pL[1])); 

    fHists[Form("dRqk_%s", cds)]->Fill(pQ[0].DeltaR(pL[1])); 
    fHists[Form("dRqk_%s", cds)]->Fill(pQ[1].DeltaR(pL[0])); 
    fHists[Form("dFqk_%s", cds)]->Fill(pQ[0].DeltaPhi(pL[1])); 
    fHists[Form("dFqk_%s", cds)]->Fill(pQ[1].DeltaPhi(pL[0])); 
    
    fHists[Form("mql_%s", cds)]->Fill((pQ[0]+pL[0]).M()); 
    fHists[Form("mql_%s", cds)]->Fill((pQ[1]+pL[1]).M()); 

    fHists[Form("mqk_%s", cds)]->Fill((pQ[0]+pL[1]).M()); 
    fHists[Form("mqk_%s", cds)]->Fill((pQ[1]+pL[0]).M()); 

    fHists[Form("mll_%s", cds)]->Fill((pL[0]+pL[1]).M()); 
    
    fHists[Form("dEtall_%s", cds)]->Fill(pL[0].Eta() - pL[1].Eta()); 
    fHists[Form("dPtll_%s", cds)]->Fill(pL[0].Pt() - pL[1].Pt()); 
  } else {
    TLorentzVector pQs, pQ[2]; pQs.SetPtEtaPhiM(fRtd.gqpt[0], fRtd.gqeta[0], fRtd.gqphi[0], 0);
    TLorentzVector pLs, pL[2]; pLs.SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
    TLorentzVector pKs;        pKs.SetPtEtaPhiM(fRtd.gkpt[0], fRtd.gketa[0], fRtd.gkphi[0], 0);

    if (fRtd.glpt[0] > fRtd.gkpt[0]) {
      pL[0].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pL[1].SetPtEtaPhiM(fRtd.gkpt[0], fRtd.gketa[0], fRtd.gkphi[0], 0);

      pQ[0].SetPtEtaPhiM(fRtd.gqpt[0], fRtd.gqeta[0], fRtd.gqphi[0], 0);
    } else {
      pL[1].SetPtEtaPhiM(fRtd.glpt[0], fRtd.gleta[0], fRtd.glphi[0], 0);
      pL[0].SetPtEtaPhiM(fRtd.gkpt[0], fRtd.gketa[0], fRtd.gkphi[0], 0);

      pQ[0].SetPtEtaPhiM(fRtd.gqpt[0], fRtd.gqeta[0], fRtd.gqphi[0], 0);
    }



    fHists[Form("ptl0_%s", cds)]->Fill(pL[0].Pt()); 
    fHists[Form("ptl1_%s", cds)]->Fill(pL[1].Pt()); 
    fHists[Form("ptq0_%s", cds)]->Fill(pQ[0].Pt()); 
    fHists[Form("ptlq_%s", cds)]->Fill((pQs + pLs).Pt()); 

    fHists[Form("etal0_%s", cds)]->Fill(pL[0].Eta()); 
    fHists[Form("etal1_%s", cds)]->Fill(pL[1].Eta()); 
    fHists[Form("etalq_%s", cds)]->Fill((pQs + pLs).Eta()); 

    fHists[Form("etaq0_%s", cds)]->Fill(pQ[0].Eta()); 


    fHists[Form("dFll_%s", cds)]->Fill(pL[0].DeltaPhi(pL[1])); 
    fHists[Form("dRll_%s", cds)]->Fill(pL[0].DeltaR(pL[1])); 
    fHists[Form("dRql_%s", cds)]->Fill(pQs.DeltaR(pLs)); 
    fHists[Form("dFql_%s", cds)]->Fill(pQs.DeltaPhi(pLs)); 
    fHists[Form("dRqk_%s", cds)]->Fill(pQs.DeltaR(pKs)); 
    fHists[Form("dFqk_%s", cds)]->Fill(pQs.DeltaPhi(pKs)); 
    
    fHists[Form("mql_%s", cds)]->Fill((pQs + pLs).M()); 
    fHists[Form("mqk_%s", cds)]->Fill((pQs + pKs).M()); 
    fHists[Form("mll_%s", cds)]->Fill((pLs + pKs).M()); 
    
    fHists[Form("dEtall_%s", cds)]->Fill(pL[0].Eta() - pL[1].Eta()); 
    fHists[Form("dPtll_%s", cds)]->Fill(pL[0].Pt() - pL[1].Pt()); 
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
  if (ifunc == 3) pF = &plotLq::loopFunction3;

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
      if (string::npos != sdecay.find("bar")) {
	ds->fColor = kBlue; 
	ds->fLcolor = kBlue; 
	ds->fFcolor = kBlue; 
	ds->fSymbol = 24; 
	ds->fFillStyle = 3356; 
      } else {
	ds->fColor = kRed; 
	ds->fLcolor = kRed; 
	ds->fFcolor = kRed; 
	ds->fSymbol = 25; 
	ds->fFillStyle = 3365; 
      }
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


