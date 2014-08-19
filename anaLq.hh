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
  ~anaLq();
  void init();
  
  void openHistFile(std::string filename);
  void closeHistFile();
  void bookHist();
  void setupReducedTree();
  
  void startAnalysis();
  void endAnalysis();
  int  loop(int nevents = 1, int start = -1);
  void eventProcessing();
  void initVariables(); 
  void setCuts(std::string cuts); 
  
  void fillHist(int type = 0);
  void fillRedTreeData(int type = 0); 
  
  // -- gen-level analysis
  void genLevelAnalysis();
  void genLQProducts(GenParticle *lq, GenParticle *l, GenParticle *q, Jet *j);
  int  genIndex(GenParticle *); 
  int isLeptonJet(Jet *j, double deltaR = 0.3); 
  double nearestLepton(Jet *j); 

  // -- reco-level analysis
  void analysis();
  void leptonSelection();
  double muonIso(Muon *m);
  void jetSelection();
  double jetMuonSeparation(Jet *j);
  void preselection();
  void lqlqSelection();
  void lqSelection();

  // -- print utilities
  void printSummary(int mode = 0); 
  void dumpGenBlock(bool withGluons=false); 
  void dumpGenJets(); 
  bool isAncestor(GenParticle *mo, GenParticle *dau);  
  void printParticle(GenParticle *); 
  void dumpDaughters(GenParticle *); 

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
  int         fL0, fL1, fJ0, fJ1; 

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
