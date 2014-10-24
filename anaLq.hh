#ifndef ANAH_H
#define ANAH_H

#include <iostream>
#include <string>
#include <map>
#include <vector>

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
#include "genLq.hh"
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
  void setupReducedTree(TTree *t);
  
  void startAnalysis();
  void endAnalysis();
  int  loop(int nevents = 1, int start = -1);
  void eventProcessing();
  void initVariables(); 
  void setCuts(std::string cuts); 
  
  void fillHist();
  void fillRedTreeData(); 
  
  // -- gen-level analysis
  void truthAnalysis();
  void genLQProducts(GenParticle *lq);
  void genSelection(); 
  genLq* createGenLQ(GenParticle *l, Jet *j, GenParticle *k = 0, Jet *i = 0); 
  void genLqAnalysis();
  void genPairAnalysis();
  void genSingleAnalysis();
  void genBgPair();
  void genBgSingle();

  int  genIndex(GenParticle *); 
  int isLeptonJet(Jet *j, double deltaR = 0.3); 
  double nearestLepton(Jet *j); 

  // -- reco-level analysis
  void leptonSelection();
  double muonIso(Muon *m);
  double electronIso(Electron *e);
  void jetSelection();
  double jetMuonSeparation(Jet *j);
  void preselection();
  void lqlqSelection();
  void lqSelection();
  int  truthMatching(lq *); 

  void candAnalysis();

  // -- print utilities
  void printSummary(int mode = 0); 
  void dumpGenBlock(bool withGluons=false, int nlines = -1); 
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

  int fEvt; 

  TFile       *fpHistFile; 
  redTreeData fRtd;

  TChain      *fpChain; 
  int         fNentries;

  std::map<std::string, TH1*> fHists;

  // -- gen LQ PAIRS and daughter particles
  std::vector<genLq *> fTrueLQ, fGenLQ; 
  std::vector<GenParticle *> fGenLeptons; 
  std::vector<Jet *> fGenJets; 
  int fNGenJets; 
  // -- reco vectors
  std::vector<lepton *> fLeptons;
  std::vector<jet *> fJets;
  std::vector<lq *> fLQ;

  // -- cuts
  int         TYPE; // 1 = single LQ production, 2 = LQ pair production
  int         CHANNEL; // 11 = electron; 13 = muon

  double      MUISODELTAR, ELISODELTAR; 
  double      MUISO, ELISO; 
  double      L0PT, L1PT;
  double      J0PT, J1PT;

  std::string fName, fTypeName;
  double      fW8;

  bool        fPreselected, fGoodEvent;
  double      fST, fMll, fMljMin;

};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
