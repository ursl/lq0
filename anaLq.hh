#ifndef ANAH_H
#define ANAH_H

#include <iostream>
#include <string>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#include "masses.hh"
#include "redTreeData.hh"
#include "lepton.hh"
#include "jet.hh"
#include "lq.hh"


// ----------------------------------------------------------------------
// -- class for analysis of Delphes trees
//    output is displayed by plotLq
// ----------------------------------------------------------------------


class anaLq {
public:
  anaLq(TChain *tree);
  virtual      ~anaLq();
  virtual void init();
  
  virtual void openHistFile(std::string filename);
  virtual void closeHistFile();
  virtual void bookHist();
  virtual void setupReducedTree();
  
  virtual void startAnalysis();
  virtual void endAnalysis();
  virtual int  loop(int nevents = 1, int start = -1);
  virtual void eventProcessing();
  virtual void initVariables(); 

  virtual void fillHist(int type = 0);
  virtual void fillRedTreeData(int type = 0); 

  // -- gen-level analysis
  virtual void genLevelAnalysis();
  virtual void genLQProducts(GenParticle *lq, GenParticle *l, GenParticle *q, Jet *j);
  virtual int  genIndex(GenParticle *); 
  virtual int isLeptonJet(Jet *j, double deltaR = 0.3); 
  virtual double nearestLepton(Jet *j); 

  // -- reco-level analysis
  virtual void analysis();
  virtual void leptonSelection();
  virtual double muonIso(Muon *m);
  virtual void jetSelection();
  virtual double jetMuonSeparation(Jet *j);
  virtual void preselection();
  virtual void lqlqSelection();
  virtual void lqSelection();

  // -- print utilities
  virtual void printSummary(int mode = 0); 
  virtual void dumpGenBlock(bool withGluons=false); 
  virtual void dumpGenJets(); 
  virtual bool isAncestor(GenParticle *mo, GenParticle *dau);  
  virtual void printParticle(GenParticle *); 
  virtual void dumpDaughters(GenParticle *); 

  // -- Delphes tree reader setup
  GenParticle* getParticle(int i) {return (GenParticle*)fbParticles->At(i);}
  Jet*         getGenJet(int i)   {return (Jet*)fbGenJets->At(i);}
  Jet*         getJet(int i)      {return (Jet*)fbJets->At(i);}
  Photon*      getPhoton(int i)   {return (Photon*)fbPhotons->At(i);}
  Electron*    getElectron(int i) {return (Electron*)fbElectrons->At(i);}
  Muon*        getMuon(int i)     {return (Muon*)fbMuons->At(i);}
  Track*       getTrack(int i)    {return (Track*)fbTracks->At(i);}
  HepMCEvent*  getEvent(int i)    {return (HepMCEvent*)fbEvent->At(i);}

  ExRootTreeReader *fTR;
  TClonesArray     
  *fbEvent,
    *fbParticles, 
    *fbGenJets, 
    *fbJets,
    *fbPhotons,  
    *fbTracks,
    *fbPFtracks, 
    *fbPFphotons, 
    *fbPFneutralh, 
    *fbElectrons,
    *fbMuons;


  TFile       *fpHistFile; 
  TTree       *fTree;
  redTreeData fRtd;

  TChain      *fpChain; 
  int         fNentries;

  // -- LQ(s) and daughter particles
  GenParticle *fGenLQp, *fGenLQpL, *fGenLQpQ;
  Jet         *fGenLQpJ;
  TLorentzVector fP4GenLQp, fP4GenLQpL, fP4GenLQpQ, fP4GenLQpJ
    , fP4GenLQpLQ, fP4GenLQpLJ;

  GenParticle *fGenLQn, *fGenLQnL, *fGenLQnQ; 
  Jet         *fGenLQnJ;
  TLorentzVector fP4GenLQn, fP4GenLQnL, fP4GenLQnQ, fP4GenLQnJ
    , fP4GenLQnLQ, fP4GenLQnLJ;

  // -- reco vectors
  std::vector<lepton *> fLeptons;
  std::vector<jet *> fJets;
  std::vector<lq *> fLQ;
  int         fPos, fNeg; 

  // -- cuts
  int         TYPE; // 1 = single LQ production, 2 = LQ pair production
  int         CHANNEL; // 11 = electron; 13 = muon
  double      MUISODELTAR; 
  double      L0PT, L1PT;
  double      J0PT, J1PT;

  int         fClass; 
  double      fW8;

  bool        fPreselected;
  double      fST, fMll, fMljetMin;

  static const int PTN   = 40; 
  static const int PTMAX = 800; 
};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
