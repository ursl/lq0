#include "anaLq.hh"

#include <sstream>

#include <TProfile.h>

#include "util.hh"

using namespace std;

// ----------------------------------------------------------------------
// Run with: ./runH -c chains/bg-test -D root 
//           ./runH -f test.root 
// ----------------------------------------------------------------------

#include "anaLq.icc"


bool sortPtP(GenParticle *a, GenParticle *b) {
  return ( a->PT > b->PT );
}

bool sortPtLepton(lepton *a, lepton *b) {
  return ( a->p4.Pt() > b->p4.Pt() );
}

bool sortPtJet(jet *a, jet *b) {
  return ( a->p4.Pt() > b->p4.Pt() );
}

bool sortPtM(Muon *a, Muon *b) {
    return ( a->PT > b->PT );
}

// ----------------------------------------------------------------------
void anaLq::setCuts(string cuts) {
  cout << "==> anaLq::setCuts: " << cuts << endl;

  istringstream ss(cuts);
  string token, name, sval;

  while (getline(ss, token, ',')) {
    
    string::size_type m1 = token.find("="); 
    name = token.substr(0, m1);
    sval = token.substr(m1+1);

    if (string::npos != name.find("CHANNEL")) {
      int val; 
      val = atoi(sval.c_str()); 
      CHANNEL = val;
    }

    if (string::npos != name.find("TYPE")) {
      int val; 
      val = atoi(sval.c_str()); 
      TYPE = val;
    }

    if (string::npos != name.find("NAME")) {
      fName = sval;
    }
    
  }

}



// ----------------------------------------------------------------------
void anaLq::startAnalysis() {
  cout << "==> anaLq: startAnalysis: " << (TYPE==2?"LQ ***pair*** production":"LQ ***single*** production") 
       << " in the ***" << (CHANNEL==13?"muon":"electron") << "*** channel with name " << fName
       << endl;

  MUISODELTAR = 0.3;
  ELISODELTAR = 0.3;

  MUISO = 0.1; 
  ELISO = 0.1;

  L0PT = 45.; 
  L1PT = 45.; 

  J0PT = 125.; 
  J1PT = 45.; 

  fTypeName = "single";
  bookHist();

  fTypeName = "pair";
  bookHist();

  fTypeName = "lq";
  bookHist();

}

// ----------------------------------------------------------------------
void anaLq::endAnalysis() {
  cout << "==> anaLq: endAnalysis: ..." << endl;
}


// ----------------------------------------------------------------------
void anaLq::eventProcessing() {
  
  if (0) cout << "----------------------------------------------------------------------" << endl;

  fPreselected = false;

  // -- fish out the CORRECT TRUE LQs and their decay/associated products
  truthAnalysis(); 

  // -- fill gen-level vectors for leptons and gen-jets (leading and subleading particles)
  genSelection(); 

  //  fTypeName = "lq";
  //  genLqAnalysis(); 
  //  fillHist();

  // -- rec-level inputs: leptons and jets
  leptonSelection();
  jetSelection();


  // -- single production analysis
  fTypeName = "single";
  genSingleAnalysis();
  preselection();
  if (fPreselected) {
    lqSelection();
    candAnalysis();
  }
  fillHist(); 

  // -- pair production analysis
  fTypeName = "pair";
  genPairAnalysis();
  preselection(); 
  if (fPreselected) {
    lqlqSelection();
    candAnalysis();
  }
  fillHist(); 

}


// ----------------------------------------------------------------------
void anaLq::truthAnalysis() {

  for (unsigned int i = 0; i < fTrueLQ.size(); ++i) delete fTrueLQ[i];
  fTrueLQ.clear();


  if (0) {
    cout << "======================================================================" << endl;
    dumpGenBlock(true); 
    cout << "----------------------------------------------------------------------" << endl;
  }

  const int LQID(9000006); 
  // -- find LQ
  GenParticle *pGen(0), *pGen0(0), *pGen1(0);

  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    pGen = getParticle(i); 
    if (LQID == TMath::Abs(pGen->PID)) {
      if (0 == pGen0) {
	pGen0 = pGen;
      } else {
	pGen1 = pGen;
	break;
      }
    }
  }

  static int firstLQ(1), firstLQbar(1); 
  if (firstLQ && ((pGen0 && pGen0->PID == LQID) || (pGen1 && pGen1->PID == LQID))) {
    firstLQ = 0; 
    cout << "---------- LQ ---------- event " << fEvt << endl;
    for (int i = 0; i < 25; ++i) {
      printParticle(getParticle(i)); 
    }
  } 
  
  if (firstLQbar && ((pGen0 && pGen0->PID == -LQID) || (pGen1 && pGen1->PID == -LQID))) {
    cout << "---------- LQbar ---------- event " << fEvt << endl;
    firstLQbar = 0; 
    for (int i = 0; i < 25; ++i) {
      printParticle(getParticle(i)); 
    }
  }

  // -- off-shell t-channel single LQ production
  fOffShell = 0; 
  if (0 == pGen0) {
    GenParticle *p1(0), *p2(0); 
    for (int i = 0; i < fbParticles->GetEntries() - 1; ++i) {
      p1 = getParticle(i); 
      p2 = getParticle(i+1); 
      bool sameMother = (p1->M1 == p2->M1);
      bool quarkMother = (p1->M1 > -1?isQuark(getParticle(p1->M1)->PID): 0);
      bool gluonMother = (p1->M1 > -1?getParticle(p1->M1)->PID == 21: 0);
      bool topMother = (p1->M1 > -1?TMath::Abs(getParticle(p1->M1)->PID) == 6: 0) ;
      if (sameMother && (quarkMother || gluonMother) && !topMother) {
	if ((isLepton(p1->PID) && isQuark(p2->PID)) || (isLepton(p2->PID) && isQuark(p1->PID))) {
	  //cout << "off-shell LQ production" << endl;
	  fOffShell = 1; 
	  if (gluonMother) fOffShell = 2; 
	  pGen0 = getParticle(p1->M1);
	  break;
	}
      }
    }
  }

  // -- Fill true LQ into fTrueLQ
  if (pGen0) genLQProducts(pGen0);
  if (pGen1) genLQProducts(pGen1);

  // -- basic histogramming
  if (pGen1) {
    fpHistFile->cd("pair");
  } else {
    fpHistFile->cd("single");
  }

  if (0 && 0 == fTrueLQ.size()) {
    cout << "NO LQ found?!" << endl;
    dumpGenBlock(true); 
  }    

  for (unsigned int i = 0; i < fTrueLQ.size(); ++i) {
    ((TH1D*)gDirectory->Get("pt"))->Fill(fTrueLQ[i]->p4LQ.Pt()); 
    ((TH1D*)gDirectory->Get("phi"))->Fill(fTrueLQ[i]->p4LQ.Phi()); 
    ((TH1D*)gDirectory->Get("eta"))->Fill(fTrueLQ[i]->p4LQ.Eta()); 
    ((TH1D*)gDirectory->Get("m"))->Fill(fTrueLQ[i]->p4LQ.M()); 

    if (fTrueLQ[i]->p4LJ.Pt() < 9990.) {
      ((TH1D*)gDirectory->Get("ljpt"))->Fill(fTrueLQ[i]->p4LJ.Pt()); 
      ((TH1D*)gDirectory->Get("ljm"))->Fill(fTrueLQ[i]->p4LJ.M()); 
    }
  }
  
  fpHistFile->cd();
}



// ----------------------------------------------------------------------
// -- not used currently!
void anaLq::genLqAnalysis() {
  
  for (unsigned int i = 0; i < fGenLQ.size(); ++i) delete fGenLQ[i];
  fGenLQ.clear();
  
  if (fGenLeptons.size() < 2) return;
  if (fGenJets.size() < 2) return;

  // -- full combinatorics
  for (unsigned int il = 0; il < 2; ++il) {
    for (unsigned int ij = 0; ij < 2; ++ij) {
      fGenLQ.push_back(createGenLQ(fGenLeptons[il], fGenJets[ij], 0, 0)); 
    }
  }

}


// ----------------------------------------------------------------------
void anaLq::genPairAnalysis() {

  for (unsigned int i = 0; i < fGenLQ.size(); ++i) delete fGenLQ[i];
  fGenLQ.clear();

  if (fGenLeptons.size() < 2) return;
  if (fGenJets.size() < 2) return;

  double mdiff[2];
  TLorentzVector p4L[2], p4J[2], p4LQ[2]; 
  for (unsigned int il = 0; il < fGenLeptons.size(); ++il) {
    p4L[il].SetPtEtaPhiM(fGenLeptons[il]->PT, fGenLeptons[il]->Eta, fGenLeptons[il]->Phi, fGenLeptons[il]->Mass);
  }

  for (unsigned int ij = 0; ij < fGenJets.size(); ++ij) {
    p4J[ij].SetPtEtaPhiM(fGenJets[ij]->PT, fGenJets[ij]->Eta, fGenJets[ij]->Phi, fGenJets[ij]->Mass);
  }

  p4LQ[0] = p4L[0] + p4J[0];
  p4LQ[1] = p4L[1] + p4J[1];
  mdiff[0] = TMath::Abs(p4LQ[0].M() - p4LQ[1].M()); 

  p4LQ[0] = p4L[0] + p4J[1];
  p4LQ[1] = p4L[1] + p4J[0];
  mdiff[1] = TMath::Abs(p4LQ[0].M() - p4LQ[1].M()); 

  if (0) {
    cout << "mdiff[0] = " << mdiff[0] << " from L " << genIndex(fGenLeptons[0]) << " plus J " << fGenJets[0] 
	 << " and L " << genIndex(fGenLeptons[1]) << " plus J " << fGenJets[1] 
	 << endl;
    cout << "mdiff[1] = " << mdiff[1] << " from L " << genIndex(fGenLeptons[0]) << " plus J " << fGenJets[1] 
	 << " and L " << genIndex(fGenLeptons[1]) << " plus J " << fGenJets[0] 
	 << endl;
  }

  if (mdiff[0] < mdiff[1]) {
    fGenLQ.push_back(createGenLQ(fGenLeptons[0], fGenJets[0], 0, 0)); 
    fGenLQ.push_back(createGenLQ(fGenLeptons[1], fGenJets[1], 0, 0)); 
  } else {
    fGenLQ.push_back(createGenLQ(fGenLeptons[0], fGenJets[1], 0, 0)); 
    fGenLQ.push_back(createGenLQ(fGenLeptons[1], fGenJets[0], 0, 0)); 
  }

}


// ----------------------------------------------------------------------
void anaLq::genSingleAnalysis() {

  for (unsigned int i = 0; i < fGenLQ.size(); ++i) delete fGenLQ[i];
  fGenLQ.clear();

  if (fGenLeptons.size() < 2) return;
  if (fGenJets.size() < 1) return;

  TLorentzVector p4L[2], p4J[2], p4LQ[2]; 
  for (unsigned int il = 0; il < fGenLeptons.size(); ++il) {
    p4L[il].SetPtEtaPhiM(fGenLeptons[il]->PT, fGenLeptons[il]->Eta, fGenLeptons[il]->Phi, fGenLeptons[il]->Mass);
  }

  for (unsigned int ij = 0; ij < fGenJets.size(); ++ij) {
    p4J[ij].SetPtEtaPhiM(fGenJets[ij]->PT, fGenJets[ij]->Eta, fGenJets[ij]->Phi, fGenJets[ij]->Mass);
  }

  p4LQ[0] = p4L[0] + p4J[0];
  fGenLQ.push_back(createGenLQ(fGenLeptons[0], fGenJets[0], fGenLeptons[1], (fGenJets.size() > 1 ? fGenJets[1] : 0))); 
}


// ----------------------------------------------------------------------
genLq* anaLq::createGenLQ(GenParticle *pL, Jet *pJ, GenParticle *pK, Jet *pI) {
  genLq *lq = new genLq(); 
  lq->q = pL->Charge;
  bool tm(false); 
  for (unsigned int i = 0; i < fTrueLQ.size(); ++i) {
    if (0 == pK) {
      if (fTrueLQ[i]->pL == pL && fTrueLQ[i]->pJ == pJ) {
	tm = true; 
	break;
      }
    } else {
      if (fTrueLQ[i]->pL == pL && fTrueLQ[i]->pJ == pJ && fTrueLQ[i]->pK == pK) {
	tm = true; 
	break;
      }
    }
  }

  lq->tm = (tm?1:0); 
  lq->q = pL->Charge;
  lq->pL = pL; 
  lq->p4L.SetPtEtaPhiM(pL->PT, pL->Eta, pL->Phi, pL->Mass);

  lq->pJ = pJ;
  lq->p4J.SetPtEtaPhiM(pJ->PT, pJ->Eta,  pJ->Phi,  pJ->Mass);
  
  TLorentzVector lj = lq->p4L + lq->p4J; 
  lq->p4LJ = lj; 

  lq->pQ = 0;
  lq->p4Q = lq->p4J; 
  lq->p4LQ = lj; 

  lq->pK = pK; 
  if (pK) {
    lq->p4K.SetPtEtaPhiM(pK->PT, pK->Eta, pK->Phi, pK->Mass);
    lq->qK = pK->Charge;
  } else {
    lq->p4K.SetPtEtaPhiM(9999., 99999., 9999., 9999.);
    lq->qK = 9999;
  }

  lq->pI = pI;
  if (pI) {
    lq->p4I.SetPtEtaPhiM(pI->PT, pI->Eta, pI->Phi, pI->Mass);
  } else {
    lq->p4I.SetPtEtaPhiM(9999., 9999., 9999., 9999.);
  }

  if (0) cout << "mass = " << lq->p4LJ.M() << " tm = " << tm 
	      << " lepton: " << genIndex(pL) << " jet: " << pJ
	      << endl;
  return lq;
}


// ----------------------------------------------------------------------
void anaLq::genSelection() {
  // -- do not delete the objects pointed to. This will screw the TClonesArray!
  fGenLeptons.clear();
  fGenJets.clear();

  Jet *pJet(0); 
  double pTj0(-1.), pTj1(-1.);
  Jet *pJ0(0), *pJ1(0); 
  double dr(0.3);
  fNGenJets = 0; 
  for (int i = 0; i < fbGenJets->GetEntries(); ++i) {
    pJet = getGenJet(i); 
    // -- remove leptons from gen jets
    int lidx = isLeptonJet(pJet, dr); 
    if (lidx > 0)  continue;
    
    if (pJet->PT > pTj0) {
      pTj1 = pTj0; 
      pJ1  = pJ0; 
      pTj0 = pJet->PT;
      pJ0  = pJet;
    } else if (pJet->PT > pTj1) {
      pTj1 = pJet->PT;
      pJ1  = pJet;
    }
    ++fNGenJets;
  }

  GenParticle *pGen(0), *pL0(0), *pL1(0); 
  double pTl0(-1.), pTl1(-1.);
  TLorentzVector m4Gen, m4L0, m4L1; 
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    pGen = getParticle(i); 
    if (13 == TMath::Abs(pGen->PID)) {
      m4Gen.SetPtEtaPhiM(pGen->PT, pGen->Eta, pGen->Phi, MMUON); 
      if (pTl0 > -1. && m4Gen.DeltaR(m4L0) < 0.1) continue;
      if (pTl1 > -1. && m4Gen.DeltaR(m4L1) < 0.1) continue;
      if (pGen->PT > pTl0) {
	pL1  = pL0; 
	pTl1 = pTl0; 
	m4L1 = m4L0; 
	pL0  = pGen;
	pTl0 = pGen->PT; 
	m4L0 = m4Gen; 
      } else if (pGen->PT > pTl1) {
	pL1  = pGen;
	pTl1 = pGen->PT; 
	m4L1 = m4Gen; 
      } 
    }
  }


  if (pL0) fGenLeptons.push_back(pL0); 
  if (pL1) fGenLeptons.push_back(pL1);
  if (pJ0) fGenJets.push_back(pJ0);
  if (pJ1) fGenJets.push_back(pJ1);

}

// ----------------------------------------------------------------------
void anaLq::leptonSelection() {

  for (unsigned int i = 0; i < fLeptons.size(); ++i) delete fLeptons[i];
  fLeptons.clear();
  
  if (13 == CHANNEL) {
    Muon *pM(0); 
    TLorentzVector m4; 
    for (int i = 0; i < fbMuons->GetEntries(); ++i) {
      pM = getMuon(i); 
      if (pM->PT < L1PT) continue;
      if (TMath::Abs(pM->Eta) > 2.4) continue;
      if (muonIso(pM) > MUISO) continue;
      lepton *l = new lepton;
      m4.SetPtEtaPhiM(pM->PT, pM->Eta, pM->Phi, MMUON); 
      l->p4        = m4;
      l->q         = pM->Charge;
      l->pMuon     = pM;
      l->pElectron = 0; 
      fLeptons.push_back(l);
    }      
  }

  if (11 == CHANNEL) {
    Electron *pE(0); 
    TLorentzVector e4; 
    for (int i = 0; i < fbElectrons->GetEntries(); ++i) {
      pE = getElectron(i); 
      if (pE->PT < L1PT) continue;
      if (TMath::Abs(pE->Eta) > 2.5) continue;
      if (electronIso(pE) > ELISO) continue;
      lepton *l = new lepton;
      e4.SetPtEtaPhiM(pE->PT, pE->Eta, pE->Phi, MELECTRON); 
      l->p4        = e4;
      l->q         = pE->Charge;
      l->pMuon     = 0;
      l->pElectron = pE; 
      fLeptons.push_back(l);
    }      
    

  }

  if (fLeptons.size() > 1) {
    sort(fLeptons.begin(), fLeptons.end(), sortPtLepton);
  }


}


// ----------------------------------------------------------------------
void anaLq::jetSelection() {

  for (unsigned int i = 0; i < fJets.size(); ++i) delete fJets[i];
  fJets.clear();

  Jet *pJ(0); 
  TLorentzVector j4; 
  for (int i = 0; i < fbJets->GetEntries(); ++i) {
    pJ = getJet(i); 
    if (pJ->PT < J1PT) continue;
    if (TMath::Abs(pJ->Eta) > 2.4) continue;
    if (jetMuonSeparation(pJ) < 0.3) continue;
    jet *j = new jet;
    j4.SetPtEtaPhiM(pJ->PT, pJ->Eta, pJ->Phi, pJ->Mass); 
    j->p4   = j4;
    j->pJet = pJ; 
    fJets.push_back(j);
  }
  
  if (fJets.size() > 1) {
    sort(fJets.begin(), fJets.end(), sortPtJet);
  }

}

// ----------------------------------------------------------------------
void anaLq::preselection() {

  fGoodEvent   = false; 

  fPreselected = false;
  fST = 9999.; 

  if ("single" == fTypeName) {
    if (fLeptons.size() < 2) return;
    if (fJets.size() < 1)    return;

    if (fLeptons[0]->p4.Pt() < L0PT) return;
    if (fJets[0]->p4.Pt() < J0PT)    return;

    if (fLeptons[0]->p4.DeltaR(fLeptons[1]->p4) < 0.3) return;

    if (fLeptons[0]->q*fLeptons[1]->q > 0) return;

    TLorentzVector ll4; 
    ll4 = fLeptons[0]->p4 + fLeptons[1]->p4;
    if (ll4.M() < 100.)  return;

    fST = fLeptons[0]->p4.Pt() + fLeptons[1]->p4.Pt() + fJets[0]->p4.Pt();
    if (fST < 0) return;

    fPreselected = true;
    return;
  }

  if ("pair" == fTypeName) {
    if (fLeptons.size() < 2) return;
    if (fJets.size() < 2)    return;

    if (fLeptons[0]->p4.Pt() < L0PT) return;
    if (fJets[0]->p4.Pt() < J0PT)    return;

    if (fLeptons[0]->p4.DeltaR(fLeptons[1]->p4) < 0.3) return;

    TLorentzVector ll4; 
    ll4 = fLeptons[0]->p4 + fLeptons[1]->p4;
    if (ll4.M() < 100.)  return;

    fST = fLeptons[0]->p4.Pt() + fLeptons[1]->p4.Pt() + fJets[0]->p4.Pt() + fJets[1]->p4.Pt();
    if (fST < 300) return;

    fPreselected = true;
    return;
  }



}

// ----------------------------------------------------------------------
// -- this is CMS' pair selection
void anaLq::lqlqSelection() {
  // -- clear LQ vector
  for (unsigned int i = 0; i < fLQ.size(); ++i) delete fLQ[i];
  fLQ.clear();
  
  int iBest(-1);
  double mdiff(99999.), mdiffBest(99999.); 
  TLorentzVector lq0, lq1; 
  for (int i = 0; i < 2; ++i) {
    lq0 = fLeptons[0]->p4 + fJets[i]->p4; 
    lq1 = fLeptons[1]->p4 + fJets[1-i]->p4; 
    //    cout << "i = " << i << " LQ0: " << lq0.M() << " LQ1: " << lq1.M() << endl;
    mdiff = TMath::Abs(lq0.M() - lq1.M()); 
    if (mdiff < mdiffBest) {
      mdiffBest = mdiff; 
      iBest = i; 
    }
  }

  lq0 = fLeptons[0]->p4 + fJets[iBest]->p4; 
  lq1 = fLeptons[1]->p4 + fJets[1-iBest]->p4; 
  //  cout << "iBest = " << iBest << " LQ0 = " << lq0.M() << "  LQ1 = " << lq1.M() << endl;

  fMll = (fLeptons[0]->p4 + fLeptons[1]->p4).M();

  // -- I hope the following is correct
  double mlq0 = lq0.M();
  double mlq1 = lq1.M();
  if (mlq0 < mlq1) {
    fMljMin = mlq0;
  } else {
    fMljMin = mlq1;
  }  


  lq *Lq = new lq;
  Lq->p4   = lq0; 
  Lq->q    = fLeptons[0]->q; 
  Lq->idxL = 0; 
  Lq->idxK = -1; 
  Lq->idxJ = iBest; 
  Lq->tm   = truthMatching(Lq);
  fLQ.push_back(Lq); 

  Lq = new lq;
  Lq->p4   = lq1; 
  Lq->q    = fLeptons[1]->q; 
  Lq->idxL = 1; 
  Lq->idxK = -1; 
  Lq->idxJ = 1-iBest; 
  Lq->tm   = truthMatching(Lq);
  fLQ.push_back(Lq); 
      
}

// ----------------------------------------------------------------------
void anaLq::lqSelection() {

  // -- clear LQ vector
  for (unsigned int i = 0; i < fLQ.size(); ++i) delete fLQ[i];
  fLQ.clear();


  int iBest(-1);
  TLorentzVector lq0 = fLeptons[0]->p4 + fJets[0]->p4; 
  TLorentzVector lq1 = fLeptons[1]->p4 + fJets[0]->p4; 
  if (lq0.M() > lq1.M()) {
    iBest = 0; 
  } else {
    iBest = 1; 
  }

  lq0 = fLeptons[iBest]->p4 + fJets[0]->p4; 
  lq *Lq = new lq;
  Lq->p4   = lq0; 
  Lq->q    = fLeptons[iBest]->q; 
  Lq->idxL = iBest; 
  Lq->idxK = 1-iBest;
  Lq->idxJ = 0;
  Lq->tm   = truthMatching(Lq);
  fLQ.push_back(Lq); 
}


// ----------------------------------------------------------------------
// -- this is the identical selection as in plotLq.cc:candAnalysis, coded w/o the reduced tree for consistency checks!
void anaLq::candAnalysis() {

  fGoodEvent   = false; 
  int idxI(0), idxJ(1); 

  if (1 == TYPE) {
    if (fLQ[idxI]->p4.M() > 0.) fGoodEvent = true;
  }

  if (2 == TYPE) {
    if (fLQ[idxI]->p4.M() > 0. && fLQ[idxJ]->p4.M() > 0.) fGoodEvent = true;
    if (fST < 685.)       fGoodEvent = false;
    if (fMll < 150.)      fGoodEvent = false;
    if (fMljMin < 155.)   fGoodEvent = false;
  } 

}


// ----------------------------------------------------------------------
void anaLq::fillHist() {

  fpHistFile->cd(fTypeName.c_str()); 
  fillRedTreeData(); 
  TTree *t = (TTree*)gDirectory->Get("events");
  t->Fill();

  if (fPreselected) { 
    fHists["pre_l0pt"]->Fill(fLeptons[0]->p4.Pt()); 
    fHists["pre_l1pt"]->Fill(fLeptons[1]->p4.Pt()); 

    fHists["pre_j0pt"]->Fill(fJets[0]->p4.Pt()); 
    if (fJets.size() > 1) fHists["pre_j1pt"]->Fill(fJets[1]->p4.Pt()); 

    fHists["pre_st"]->Fill(fRtd.st); 
    fHists["pre_mll"]->Fill(fMll); 
    fHists["pre_mljetmin"]->Fill(fMljMin); 
      
    if (fLQ[0]->p4.M() > 0.) {
      fHists["pre_m"]->Fill(fLQ[0]->p4.M()); 
      fHists["pre_pt"]->Fill(fLQ[0]->p4.Pt()); 
      if (2 == TYPE) {
	fHists["pre_m"]->Fill(fLQ[1]->p4.M()); 
	fHists["pre_pt"]->Fill(fLQ[1]->p4.Pt()); 
      }
    }

    if (fGoodEvent) { 
      fHists["sel_st"]->Fill(fRtd.st); 
      fHists["sel_mll"]->Fill(fMll); 
      fHists["sel_mljetmin"]->Fill(fMljMin); 
      
      fHists["sel_m"]->Fill(fLQ[0]->p4.M()); 
      fHists["sel_pt"]->Fill(fLQ[0]->p4.Pt()); 
      if (2 == TYPE) {
	fHists["sel_m"]->Fill(fLQ[1]->p4.M()); 
	fHists["sel_pt"]->Fill(fLQ[1]->p4.Pt()); 
      }
    }

  }

}



// ----------------------------------------------------------------------
void anaLq::bookHist() {
  fpHistFile->cd();			   
  cout << "==> anaLq: bookHist in directory " << fTypeName << endl;

  TH1D *h1(0); 
  (void)h1;

  fpHistFile->mkdir(Form("%s", fTypeName.c_str())); 
  fpHistFile->cd(Form("%s", fTypeName.c_str()));

  // -- Reduced Tree
  TTree *t = new TTree("events", "events");
  setupReducedTree(t);

  // -- Histograms
  h1 = new TH1D("pt",  "lq gen pt", 150, 0., 1500.); 
  h1 = new TH1D("phi", "lq gen phi", 40, -3.15, 3.15); 
  h1 = new TH1D("eta", "lq gen eta", 40, -5., 5.); 
  h1 = new TH1D("m",   "lq gen m", 100, 0., 2000.); 

  h1 = new TH1D("ljpt",  "l+j gen pt", 150, 0., 1500); 
  h1 = new TH1D("ljm",   "l+j gen m", 100, 0., 2000.); 


  vector<string> levels; 
  levels.push_back("pre"); 
  levels.push_back("sel"); 

  for (unsigned int i = 0; i < levels.size(); ++i) {
    // -- histograms as in plotLq:
    // -- m
    fHists.insert(make_pair(Form("%s_m", levels[i].c_str()), 
			    new TH1D(Form("%s_m", levels[i].c_str()), 
				     Form("%s_m", levels[i].c_str()), 
				     40, 0, 2000.))); 
    setTitles(fHists[Form("%s_m", levels[i].c_str())], "m [GeV]", "Entries/bin");
    
    // -- st
    fHists.insert(make_pair(Form("%s_st", levels[i].c_str()), 
			    new TH1D(Form("%s_st", levels[i].c_str()), 
				     Form("%s_st", levels[i].c_str()), 
				     14, 0, 3500.))); 
    setTitles(fHists[Form("%s_st", levels[i].c_str())], "S_{T} [GeV]", "Entries/bin");
    
    // -- mll
    fHists.insert(make_pair(Form("%s_mll", levels[i].c_str()), 
			    new TH1D(Form("%s_mll", levels[i].c_str()), 
				     Form("%s_mll", levels[i].c_str()), 
				     30, 0, 1500.))); 
    setTitles(fHists[Form("%s_mll", levels[i].c_str())], "m_{l l} [GeV]", "Entries/bin");
    
    // -- mljetmin
    fHists.insert(make_pair(Form("%s_mljetmin", levels[i].c_str()), 
			    new TH1D(Form("%s_mljetmin", levels[i].c_str()), 
				     Form("%s_mljetmin", levels[i].c_str()),
				     15, 0, 1500.))); 
    setTitles(fHists[Form("%s_mljetmin", levels[i].c_str())], "m_{l jet}^{min} [GeV]", "Entries/bin");
    
    
    // -- pt
    fHists.insert(make_pair(Form("%s_pt", levels[i].c_str()), 
			    new TH1D(Form("%s_pt", levels[i].c_str()), 
				     Form("%s_pt", levels[i].c_str()), 
				     100, 0, 1000.))); 
    setTitles(fHists[Form("%s_pt", levels[i].c_str())], "p_{T} [GeV]", "Entries/bin");

    // -- lepton pt
    fHists.insert(make_pair(Form("%s_l0pt", levels[i].c_str()), 
			    new TH1D(Form("%s_l0pt", levels[i].c_str()), 
				     Form("%s_l0pt", levels[i].c_str()), 
				     64, 0, 1600.))); 
    setTitles(fHists[Form("%s_l0pt", levels[i].c_str())], "p_{T}(l_{0}) [GeV]", "Entries/bin");

    fHists.insert(make_pair(Form("%s_l1pt", levels[i].c_str()), 
			    new TH1D(Form("%s_l1pt", levels[i].c_str()), 
				     Form("%s_l1pt", levels[i].c_str()), 
				     40, 0, 800.))); 
    setTitles(fHists[Form("%s_l1pt", levels[i].c_str())], "p_{T}(l_{1}) [GeV]", "Entries/bin");

    // -- jet pt
    fHists.insert(make_pair(Form("%s_j0pt", levels[i].c_str()), 
			    new TH1D(Form("%s_j0pt", levels[i].c_str()), 
				     Form("%s_j0pt", levels[i].c_str()), 
				     64, 0, 1600.))); 
    setTitles(fHists[Form("%s_j0pt", levels[i].c_str())], "p_{T}(j_{0}) [GeV]", "Entries/bin");

    fHists.insert(make_pair(Form("%s_j1pt", levels[i].c_str()), 
			    new TH1D(Form("%s_j1pt", levels[i].c_str()), 
				     Form("%s_j1pt", levels[i].c_str()), 
				     40, 0, 800.))); 
    setTitles(fHists[Form("%s_j1pt", levels[i].c_str())], "p_{T}(j_{1}) [GeV]", "Entries/bin");
  }
  
}


// ----------------------------------------------------------------------
void anaLq::initVariables() {
  //  fW8 = getEvent(0)->Weight;

}

// ----------------------------------------------------------------------
void anaLq::printParticle(GenParticle *p) {
  if (p) {
    for (int i = 0; i < fbParticles->GetEntries(); ++i) {
      if (p == getParticle(i)) {
	cout << Form("%4d %+8d  M: %4d %4d D: %4d %4d P: %+9.3f %+9.3f %+9.3f PT: %+9.3f %+9.3f %+9.3f" , 
		     i, p->PID, p->M1, p->M2, p->D1, p->D2, 
		     p->Px, p->Py, p->Pz,
		     p->PT, p->Eta, p->Phi
		     ) 
	     << endl;
	break;
      }  
    }
  }
}


// ----------------------------------------------------------------------
int anaLq::genIndex(GenParticle *p) {
  int idx(-1); 
  if (p) {
    for (int i = 0; i < fbParticles->GetEntries(); ++i) {
      if (p == getParticle(i)) {
	idx = i; 
	break;
      } 
    }
  }
  return idx; 
}


// ----------------------------------------------------------------------
// FIXME: This does not work for jets!
bool anaLq::isAncestor(GenParticle *mo, GenParticle *dau) {

  if (mo == dau)  return true;

  GenParticle *pGen(dau); 
  while ((pGen = getParticle(pGen->M1))) {
    if (mo == pGen) {
      return true;
    }
    if (pGen->M1 < 0) break;
    if (pGen->M1 > fbParticles->GetEntries()) break;
  }
  return false;
}

// ----------------------------------------------------------------------
void anaLq::dumpGenBlock(bool withGluons, int nlines) {

  if (nlines < 0) nlines = fbParticles->GetEntries();
  cout << "=== GenBlock for event " << fEvt << " ===" << endl;
  for (int i = 0; i < nlines; ++i) {
    if (21 == getParticle(i)->PID) {
      if (withGluons) printParticle(getParticle(i)); 
    } else {
      printParticle(getParticle(i)); 
    }
  }

}

// ----------------------------------------------------------------------
void anaLq::dumpGenJets() {

  Jet *pJet(0); 
  int il(0); 
  string sl("");
  cout << "=== GenJets ===" << endl;
  for (int i = 0; i < fbGenJets->GetEntries(); ++i) {
    pJet = getGenJet(i); 
    il = isLeptonJet(pJet);
    if (il > 0) {
      sl = Form("lepton %d", il); 
    } else {
      sl = Form("0x%lx", (void*)pJet); 
    }

    cout << Form("%4d PT: %+9.3f %+9.3f %+9.3f Mass: %7.3f dE: %5.4f dF: %5.4f n+: %d n0: %d %s" , 
		 i, pJet->PT, pJet->Eta, pJet->Phi, pJet->Mass,
		 pJet->DeltaEta, pJet->DeltaPhi, pJet->NCharged, pJet->NNeutrals, 
		 sl.c_str()
		 ) 
	 << endl;
  }
}



// ----------------------------------------------------------------------
void anaLq::dumpDaughters(GenParticle *pMom) {
  // dysfunctional with the MadGrap/DELPHES gen block
  
}

// ----------------------------------------------------------------------
void anaLq::setupReducedTree(TTree *t) {
  t->Branch("type",    &fRtd.type,            "type/I");
  t->Branch("w8",      &fRtd.w8,              "w8/D");

  t->Branch("ngen",    &fRtd.ngen,            "ngen/I");
  t->Branch("gm",       fRtd.gm,              "gm[ngen]/D");
  t->Branch("gpt",      fRtd.gpt,             "gpt[ngen]/D");
  t->Branch("geta",     fRtd.geta,            "geta[ngen]/D");
  t->Branch("gphi",     fRtd.gphi,            "gphi[ngen]/D");
  t->Branch("gmlj",     fRtd.gmlj,            "gmlj[ngen]/D");
  t->Branch("glq",      fRtd.glq,             "glq[ngen]/I");
  t->Branch("gtm",      fRtd.gtm,             "gtm[ngen]/I");

  t->Branch("glpt",     fRtd.glpt,           "glpt[ngen]/D");
  t->Branch("gleta",    fRtd.gleta,          "gleta[ngen]/D");
  t->Branch("glphi",    fRtd.glphi,          "glphi[ngen]/D");

  t->Branch("gqpt",     fRtd.gqpt,           "gqpt[ngen]/D");
  t->Branch("gqeta",    fRtd.gqeta,          "gqeta[ngen]/D");
  t->Branch("gqphi",    fRtd.gqphi,          "gqphi[ngen]/D");

  t->Branch("gjpt",     fRtd.gjpt,           "gjpt[ngen]/D");
  t->Branch("gjeta",    fRtd.gjeta,          "gjeta[ngen]/D");
  t->Branch("gjphi",    fRtd.gjphi,          "gjphi[ngen]/D");

  t->Branch("gkpt",     fRtd.gkpt,           "gkpt[ngen]/D");
  t->Branch("gketa",    fRtd.gketa,          "gketa[ngen]/D");
  t->Branch("gkphi",    fRtd.gkphi,          "gkphi[ngen]/D");
  t->Branch("gkq",      fRtd.gkq,            "gkq[ngen]/I");

  t->Branch("gipt",     fRtd.gipt,           "gipt[ngen]/D");
  t->Branch("gieta",    fRtd.gieta,          "gieta[ngen]/D");
  t->Branch("giphi",    fRtd.giphi,          "giphi[ngen]/D");

  t->Branch("nrec",    &fRtd.nrec,            "nrec/I");
  t->Branch("m",        fRtd.m,               "m[nrec]/D");
  t->Branch("pt",       fRtd.pt,              "pt[nrec]/D");
  t->Branch("eta",      fRtd.eta,             "eta[nrec]/D");
  t->Branch("phi",      fRtd.phi,             "phi[nrec]/D");
  t->Branch("lq",       fRtd.lq,              "lq[nrec]/I");
  t->Branch("tm",       fRtd.tm,              "tm[nrec]/I");

  t->Branch("lpt",      fRtd.lpt,              "lpt[nrec]/D");
  t->Branch("leta",     fRtd.leta,             "leta[nrec]/D");
  t->Branch("lphi",     fRtd.lphi,             "lphi[nrec]/D");

  t->Branch("jpt",      fRtd.jpt,              "jpt[nrec]/D");
  t->Branch("jeta",     fRtd.jeta,             "jeta[nrec]/D");
  t->Branch("jphi",     fRtd.jphi,             "jphi[nrec]/D");

  t->Branch("kpt",      fRtd.kpt,              "kpt[nrec]/D");
  t->Branch("keta",     fRtd.keta,             "keta[nrec]/D");
  t->Branch("kphi",     fRtd.kphi,             "kphi[nrec]/D");

  t->Branch("st",       &fRtd.st,              "st/D");
  t->Branch("mll",      &fRtd.mll,             "mll/D");
  t->Branch("mljmin",   &fRtd.mljmin,          "mljmin/D");



}


// ----------------------------------------------------------------------
void anaLq::fillRedTreeData() {

  fRtd.type    = TYPE; 
  fRtd.channel = CHANNEL; 
  fRtd.w8      = fW8;

  for (unsigned int i = 0; i < fTrueLQ.size(); ++i) {
    fRtd.gm[i]   = fTrueLQ[i]->p4LQ.M(); 
    fRtd.gpt[i]  = fTrueLQ[i]->p4LQ.Pt(); 
    fRtd.geta[i] = fTrueLQ[i]->p4LQ.Eta(); 
    fRtd.gphi[i] = fTrueLQ[i]->p4LQ.Phi(); 

    fRtd.glq[i]  = fTrueLQ[i]->q; 
    fRtd.gtm[i]  = fTrueLQ[i]->tm; 

    fRtd.glpt[i]  = fTrueLQ[i]->p4L.Pt(); 
    fRtd.gleta[i] = fTrueLQ[i]->p4L.Eta(); 
    fRtd.glphi[i] = fTrueLQ[i]->p4L.Phi(); 

    fRtd.gqpt[i]  = fTrueLQ[i]->p4Q.Pt(); 
    fRtd.gqeta[i] = fTrueLQ[i]->p4Q.Eta(); 
    fRtd.gqphi[i] = fTrueLQ[i]->p4Q.Phi(); 

    if (fTrueLQ[i]->pJ) {
      fRtd.gmlj[i]  = fTrueLQ[i]->p4LJ.M(); 
      fRtd.gjpt[i]  = fTrueLQ[i]->p4J.Pt(); 
      fRtd.gjeta[i] = fTrueLQ[i]->p4J.Eta(); 
      fRtd.gjphi[i] = fTrueLQ[i]->p4J.Phi(); 
    } else {
      fRtd.gmlj[i]  = 9999.;
    }

    if (fTrueLQ[i]->pK) {
      fRtd.gkpt[i]  = fTrueLQ[i]->p4K.Pt(); 
      fRtd.gketa[i] = fTrueLQ[i]->p4K.Eta(); 
      fRtd.gkphi[i] = fTrueLQ[i]->p4K.Phi(); 
      fRtd.gkq[i]   = fTrueLQ[i]->qK; 
    } else {
      fRtd.gkpt[i]  = 9999.;
      fRtd.gketa[i] = 9999.;
      fRtd.gkphi[i] = 9999.;
      fRtd.gkq[i]   = 9999;
    }
  }
  fRtd.ngen = fTrueLQ.size(); 

  if (!fPreselected) {
    fRtd.nrec = 0; 
    return;
  }

  // -- and now fill all LQs (pair production characterized by kpt < 0; single production with kpt > 0)
  for (unsigned int i = 0; i < fLQ.size(); ++i) {
    fRtd.m[i]    = fLQ[i]->p4.M();
    fRtd.pt[i]   = fLQ[i]->p4.Pt();
    fRtd.eta[i]  = fLQ[i]->p4.Eta();
    fRtd.phi[i]  = fLQ[i]->p4.Phi();
    
    fRtd.tm[i]    = fLQ[i]->tm;
    fRtd.lq[i]    = fLeptons[fLQ[i]->idxL]->q;
    fRtd.lpt[i]   = fLeptons[fLQ[i]->idxL]->p4.Pt();
    fRtd.leta[i]  = fLeptons[fLQ[i]->idxL]->p4.Phi();
    fRtd.lphi[i]  = fLeptons[fLQ[i]->idxL]->p4.Eta();
    
    fRtd.jpt[i]   = fJets[fLQ[i]->idxJ]->p4.Pt();
    fRtd.jeta[i]  = fJets[fLQ[i]->idxJ]->p4.Phi();
    fRtd.jphi[i]  = fJets[fLQ[i]->idxJ]->p4.Eta();
    
    if (fLQ[i]->idxK > -1) {
      fRtd.kpt[i]   = fLeptons[fLQ[i]->idxK]->p4.Pt();
      fRtd.keta[i]  = fLeptons[fLQ[i]->idxK]->p4.Phi();
      fRtd.kphi[i]  = fLeptons[fLQ[i]->idxK]->p4.Eta();
    } else {
      fRtd.kpt[i]   = 9999.;
      fRtd.keta[i]  = 9999.;
      fRtd.kphi[i]  = 9999.;
    }
  }

  fRtd.nrec = fLQ.size(); 
  fRtd.st     = fST; 
  fRtd.mll    = fMll; 
  fRtd.mljmin = fMljMin; 
}


// ----------------------------------------------------------------------
void anaLq::genLQProducts(GenParticle *lq) {

  GenParticle *l(0); 
  GenParticle *q(0); 
  Jet *j(0); 
  int lqIdx = genIndex(lq); 

  // -- find gen-level direct daughters of LQ
  GenParticle *pGen(0);
  for (int i = lqIdx+1; i < fbParticles->GetEntries(); ++i) {
    pGen = getParticle(i); 
    int pid = pGen->PID; 
    if (pGen->M1 == lqIdx) {
      if (isLepton(pid)) {
	l = pGen; 
      } else {
	q = pGen;
      }
    }
    if (l != 0 && q != 0) break; // both lepton and quark have been found
  }

  if (0 == q || 0 == l) {
    cout << "XXXXX LQ found, but decay quark or lepton missing?!" << endl;
    return;
  }

  TLorentzVector p4Q; 
  p4Q.SetPtEtaPhiM(q->PT, q->Eta, q->Phi, q->Mass); 

  // -- find gen-jet for quark
  // FIXME use "TRefArray Particles" instead of dR matching!
  Jet *pJet(0); 
  double dRMin(9999.);
  int dRBest(-1); 
  double dR(0.);
  TLorentzVector p4J; 
  double dr(0.3);
  for (int i = 0; i < fbGenJets->GetEntries(); ++i) {
    pJet = getGenJet(i); 
    // -- remove leptons from gen jets
    if (isLeptonJet(pJet, dr)) {
      continue;
    }

    p4J.SetPtEtaPhiM(pJet->PT, pJet->Eta, pJet->Phi, pJet->Mass);
    dR = p4Q.DeltaR(p4J); 
    if (dR < dRMin) {
      dRMin = dR; 
      dRBest = i; 
    }
    
    // FIXME: Need to fix isAncestor to work with jets!
    //     GenParticle *pGen = (GenParticle*)pJet->Particles.At(0);
    //     cout << "genJet " << i << " has npart = " << pJet->Particles.GetEntriesFast()  
    // 	 << ", 0 is at index " << genIndex(pGen)
    // 	 << endl;
    //     if (isAncestor(lq, pGen)) cout << "  THIS IS FROM LQ " << lqIdx << endl;

  }    

  // -- so far choose on best delta(R) match!
  if (dRBest > -1 && dRMin < 0.5) {
    j = getGenJet(dRBest);
    //    cout << " xx choosing jet = " << dRBest << " with ET = " << j->PT << " for jet of LQ " << lqIdx << endl;
  }

  // -- create and fill struct
  genLq *gen = new genLq(); 
  gen->q  = l->Charge;
  gen->tm = fOffShell; 

  gen->pL = l; 
  gen->p4L.SetPtEtaPhiM(l->PT, l->Eta, l->Phi, l->Mass);
  
  gen->pQ = q; 
  gen->p4Q.SetPtEtaPhiM(q->PT, q->Eta, q->Phi, q->Mass);

  // -- this also covers off-shell LQ production, where the "mother" of the L and Q is an off-shell quark, which has not much to do with the (L+Q) system)
  TLorentzVector LQ = gen->p4L + gen->p4Q; 
  gen->pLQ = lq;
  gen->p4LQ = LQ;

  if (j) {
    gen->pJ = j;
    gen->p4J.SetPtEtaPhiM(j->PT, j->Eta,  j->Phi,  j->Mass);
    
    TLorentzVector lj = gen->p4L + gen->p4J; 
    gen->p4LJ = lj; 
  } else {
    gen->pJ = 0;
    gen->p4J.SetPtEtaPhiM(9999., 9999., 9999., 9999.); 
    gen->p4LJ.SetPtEtaPhiM(9999., 9999., 9999., 9999.); 
  }

  // -- find POSSIBLE gen-level bachelor lepton K to LQ
  int lqMomIdx = lq->M1;
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    if (i == lqIdx) continue;
    pGen = getParticle(i); 
    if (isLepton(pGen->PID) && ((pGen->M1 == lqMomIdx) || (pGen->M1 == lqIdx + 1)  || (pGen->M1 == lqIdx - 1))) {
      break;
    } else {
      pGen = 0; 
    }
  }

  if (pGen) {
    gen->pK = pGen;
    gen->p4K.SetPtEtaPhiM(pGen->PT, pGen->Eta, pGen->Phi, pGen->Mass);
    gen->qK = pGen->Charge;
  } else {
    gen->pK = 0;
    gen->p4K.SetPtEtaPhiM(9999., 9999., 9999., 9999.); 
    gen->qK = 9999;
  }

  fTrueLQ.push_back(gen);

  if (0) cout << "Truth LQ lepton = " << genIndex(l) << " quark = " << genIndex(q) << " jet = " << j 
	      << " and mass = " << gen->p4LQ.M() << (fOffShell == 1? " off-shell": "") 
	      << endl;

  // -- printout for validation of background processing
  if (0 && 1 == fOffShell) {
    cout << "---------- ?? ---------- event " << fEvt << " LQ at " << lqIdx << " lepton = " << genIndex(l) << " quark = " << genIndex(q) << endl;
    for (int i = 0; i < 100; ++i) {
      printParticle(getParticle(i)); 
    }
  }    

}



// ----------------------------------------------------------------------
void anaLq::genBgPair() {
  
  Jet *pJet(0); 
  double pTj0(-1.), pTj1(-1.);
  double pTl0(-1.), pTl1(-1.);
  int j0(-1), j1(-1);
  int l0(-1), l1(-1);
  double dr(0.3);
  for (int i = 0; i < fbGenJets->GetEntries(); ++i) {
    pJet = getGenJet(i); 
    // -- remove leptons from gen jets
    int lidx = isLeptonJet(pJet, dr); 
    if (lidx > 0) {
      GenParticle *pGen = getParticle(lidx); 
      if (13 == TMath::Abs(pGen->PID)) {
	if (pGen->PT > pTl0) {
	  pTl0 = pGen->PT;
	  l0 = lidx;
	} else if (pGen->PT > pTl1) {
	  pTl1 = pGen->PT;
	  l1 = lidx; 
	}
      }
      continue;
    }
    
    if (pJet->PT > pTj0) {
      pTj0 = pJet->PT;
      j0 = i; 
    } else if (pJet->PT > pTj1) {
      pTj1 = pJet->PT;
      j1 = i; 
    }
  }

  if (l1 > -1 && j1 > -1) {
    if (0) {
      cout << "Leptons: " << endl;
      if (l0 > 0) printParticle(getParticle(l0));
      if (l1 > 0) printParticle(getParticle(l1));
      
      cout << "Jets: " << endl;
      if (j0 > -1) cout << getGenJet(j0)->PT << " " << getGenJet(j0)->Eta << " " << getGenJet(j0)->Phi << endl;
      if (j1 > -1) cout << getGenJet(j1)->PT << " " << getGenJet(j1)->Eta << " " << getGenJet(j1)->Phi << endl;
    }
  } else {
    return;
  }

  TLorentzVector p4L0, p4L1, p4J0, p4J1; 
  GenParticle *pL0 = getParticle(l0);
  p4L0.SetPtEtaPhiM(pL0->PT, pL0->Eta, pL0->Phi, pL0->Mass);
  GenParticle *pL1 = getParticle(l1);
  p4L1.SetPtEtaPhiM(pL1->PT, pL1->Eta, pL1->Phi, pL1->Mass);

  Jet *pJ0 = getGenJet(j0);
  p4J0.SetPtEtaPhiM(pJ0->PT, pJ0->Eta, pJ0->Phi, pJ0->Mass);

  Jet *pJ1 = getGenJet(j1);
  p4J1.SetPtEtaPhiM(pJ1->PT, pJ1->Eta, pJ1->Phi, pJ1->Mass);

  // -- combine the leptons with the best-matching jet
  TLorentzVector p4Lq0, p4Lq1; 
  double mdiff0(99999.), mdiff1(99999.); 
  TLorentzVector lq0, lq1; 
  lq0 = p4L0 + p4J0;
  lq1 = p4L1 + p4J1;
  mdiff0 = TMath::Abs(lq0.M() - lq1.M()); 

  lq0 = p4L0 + p4J1;
  lq1 = p4L1 + p4J0;
  mdiff1 = TMath::Abs(lq0.M() - lq1.M()); 

  genLq *fake0 = new genLq(); 
  fake0->q = pL0->Charge;

  fake0->pL = pL0; 
  fake0->p4L.SetPtEtaPhiM(pL0->PT, pL0->Eta, pL0->Phi, pL0->Mass);

  fake0->pK = 0; 
  fake0->p4K.SetPtEtaPhiM(9999., 9999., 9999., 9999.);

  genLq *fake1 = new genLq(); 
  fake1->q = pL1->Charge;

  fake1->pL = pL1; 
  fake1->p4L.SetPtEtaPhiM(pL1->PT, pL1->Eta, pL1->Phi, pL1->Mass);

  fake1->pK = 0; 
  fake1->p4K.SetPtEtaPhiM(9999., 9999., 9999., 9999.);

  if (mdiff0 < mdiff1) {
    p4Lq0 = p4L0 + p4J0;
    p4Lq1 = p4L1 + p4J1;

    fake0->pLQ = 0;
    fake0->p4LQ = p4Lq0;
    fake0->p4LJ = p4Lq0;

    fake0->pQ = 0; 
    fake0->p4Q = p4J0;

    fake0->pJ = pJ0;
    fake0->p4J = p4J0;

    fake1->pLQ = 0;
    fake1->p4LQ = p4Lq1;
    fake1->p4LJ = p4Lq1;

    fake1->pQ = 0; 
    fake1->p4Q = p4J1;

    fake1->pJ = pJ1;
    fake1->p4J = p4J1;
  } else {
    p4Lq0 = p4L0 + p4J1;
    p4Lq1 = p4L1 + p4J0;

    fake0->pLQ = 0;
    fake0->p4LQ = p4Lq0;
    fake0->p4LJ = p4Lq0;

    fake0->pQ = 0; 
    fake0->p4Q = p4J1;

    fake0->pJ = pJ0;
    fake0->p4J = p4J1;

    fake1->pLQ = 0;
    fake1->p4LQ = p4Lq1;
    fake1->p4LJ = p4Lq1;

    fake1->pQ = 0; 
    fake1->p4Q = p4J0;

    fake1->pJ = pJ1;
    fake1->p4J = p4J0;
  }

  fGenLQ.push_back(fake0);
  fGenLQ.push_back(fake1);



}


// ----------------------------------------------------------------------
void anaLq::genBgSingle() {

}





// ----------------------------------------------------------------------
int anaLq::isLeptonJet(Jet *j, double deltaR) {

  TLorentzVector p4J; 
  p4J.SetPtEtaPhiM(j->PT, j->Eta, j->Phi, j->Mass);

  GenParticle *pGen(0);
  TLorentzVector p4; 
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    pGen = getParticle(i); 
    int pid = pGen->PID; 
    if (isLepton(pid)) {
      p4.SetPtEtaPhiM(pGen->PT, pGen->Eta, pGen->Phi, pGen->Mass); 
      if (TMath::Abs(p4J.Pt() - p4.Pt())/p4J.Pt() > 0.5) continue; // skip leptons with less than 50% of the jet pT
      if (p4J.DeltaR(p4) < deltaR) {
	return i;
      }
    }
  }
  return 0;
}


// ----------------------------------------------------------------------
double anaLq::nearestLepton(Jet *j) {

  TLorentzVector p4J; 
  p4J.SetPtEtaPhiM(j->PT, j->Eta, j->Phi, j->Mass);

  GenParticle *pGen(0);
  TLorentzVector p4; 
  double dr(99.), drMin(99.); 
  int ibest(-1);
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    pGen = getParticle(i); 
    int pid = pGen->PID; 
    if (isLepton(pid)) {
      p4.SetPtEtaPhiM(pGen->PT, pGen->Eta, pGen->Phi, pGen->Mass); 
      dr = p4J.DeltaR(p4);
      if (dr < drMin) {
	drMin = dr;
	ibest = i;
      }
    }
  }
  if (drMin < 0.5) cout << "best lepton matching at bar code " << ibest << " with dr = " << drMin << endl;
  return drMin;
}


// ----------------------------------------------------------------------
void anaLq::printSummary(int mode) {

  if (mode & 0x1) dumpGenBlock(); 
  if (mode & 0x2) dumpGenJets();

}

// ----------------------------------------------------------------------
double anaLq::muonIso(Muon *m) {

  double iso(0.);
  Track *pT(0);
  TLorentzVector m4, t4; 
  m4.SetPtEtaPhiM(m->PT, m->Eta, m->Phi, MMUON); 
  for (int i = 0; i < fbTracks->GetEntries(); ++i) {
    pT = getTrack(i); 
    t4.SetPtEtaPhiM(pT->PT, pT->Eta, pT->Phi, MPION); 
    if (pT->Particle == m->Particle) {
      continue;
    }
    if (m4.DeltaR(t4) < MUISODELTAR) {
      iso += t4.Pt();
    }
  }
  iso /= m->PT; 
  return iso;
}

// ----------------------------------------------------------------------
double anaLq::electronIso(Electron *e) {

  double iso(0.);
  Track *pT(0);
  TLorentzVector e4, t4; 
  e4.SetPtEtaPhiM(e->PT, e->Eta, e->Phi, MELECTRON); 
  for (int i = 0; i < fbTracks->GetEntries(); ++i) {
    pT = getTrack(i); 
    t4.SetPtEtaPhiM(pT->PT, pT->Eta, pT->Phi, MPION); 
    if (pT->Particle == e->Particle) {
      continue;
    }
    if (e4.DeltaR(t4) < ELISODELTAR) {
      iso += t4.Pt();
    }
  }
  iso /= e->PT; 
  return iso;
}


// ----------------------------------------------------------------------
double anaLq::jetMuonSeparation(Jet *j) {
  Muon *pM(0); 
  TLorentzVector j4, m4; 
  j4.SetPtEtaPhiM(j->PT, j->Eta, j->Phi, j->Mass); 

  double sep(99.), sepMin(99.); 
  for (int i = 0; i < fbMuons->GetEntries(); ++i) {
    pM = getMuon(i); 
    if (TMath::Abs(pM->Eta) > 2.1) continue;

    m4.SetPtEtaPhiM(pM->PT, pM->Eta, pM->Phi, MMUON); 
    sep = j4.DeltaR(m4); 
    if (sep < sepMin) {
      sepMin = sep; 
    }
  }

  return sepMin;
}


// ----------------------------------------------------------------------
int anaLq::truthMatching(lq *Lq) {
  
  TLorentzVector j = fJets[Lq->idxJ]->p4;
  TLorentzVector l = fLeptons[Lq->idxL]->p4;
  TLorentzVector k = (Lq->idxK > -1? fLeptons[Lq->idxK]->p4: TLorentzVector(0., 0., 0., 0.));

  for (unsigned int i = 0; i < fTrueLQ.size(); ++i) {
    TLorentzVector gl = fTrueLQ[i]->p4L;
    TLorentzVector gj = fTrueLQ[i]->p4Q;
    TLorentzVector gk = fTrueLQ[i]->p4K;

    if (k.Pt() > 0.1) {
      if (l.DeltaR(gl) < 0.3 && j.DeltaR(gj) < 0.3 &&k.DeltaR(gk) < 0.3) {
	//	cout << "truthmatched to true single LQ " << i << endl;
	return i;
      } else {
	// 	cout << "L: pT = " << l.Pt() << "/" << gl.Pt() << ", eta = " << l.Eta() << "/" << gl.Eta() << ", phi = " << l.Phi() << "/" << gl.Phi() << endl;
	// 	cout << "J: pT = " << j.Pt() << "/" << gj.Pt() << ", eta = " << j.Eta() << "/" << gj.Eta() << ", phi = " << j.Phi() << "/" << gj.Phi() << endl;
	// 	cout << "K: pT = " << k.Pt() << "/" << gk.Pt() << ", eta = " << k.Eta() << "/" << gk.Eta() << ", phi = " << k.Phi() << "/" << gk.Phi() << endl;
      }
    } else {
      //      cout << "K: pT = " << k.Pt() << endl;
      if (l.DeltaR(gl) < 0.3 && j.DeltaR(gj) < 0.3) {
	//      cout << "truthmatched to true pair LQ " << i << endl;
	return i;
      }
    }
  }
  return -1;
}
