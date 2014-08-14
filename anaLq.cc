#include "anaLq.hh"

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

bool sortPtE(Electron *a, Electron *b) {
    return ( a->PT > b->PT );
}

bool sortPtM(Muon *a, Muon *b) {
    return ( a->PT > b->PT );
}

// ----------------------------------------------------------------------
void anaLq::startAnalysis() {
  cout << "anaLq: startAnalysis: ..." << endl;
}

// ----------------------------------------------------------------------
void anaLq::endAnalysis() {
  cout << "anaLq: endAnalysis: ..." << endl;
}


// ----------------------------------------------------------------------
void anaLq::eventProcessing() {

  initVariables();

  //   cout << "w8: " <<  getEvent(0)->Weight 
  //        << " processid: " << getEvent(0)->ProcessID
  //        << " MPI: " << getEvent(0)->MPI
  //        << " X1: " << getEvent(0)->X1 
  //        << " X2: " << getEvent(0)->X2
  //        << " qcd: " << getEvent(0)->AlphaQCD
  //        << " nevents = " << fbEvent->GetEntries() 
  //        << endl;

  genLevelAnalysis(); 
  fillHist(); 

}


// ----------------------------------------------------------------------
void anaLq::genLevelAnalysis() {

  static int first(1); 
  if (0) {
    cout << "======================================================================" << endl;
    dumpGenBlock(true); 
    cout << "----------------------------------------------------------------------" << endl;
  }

  if (first) {
    first = 0; 
    for (int i = 0; i < 15; ++i) {
      printParticle(getParticle(i)); 
    }
  }
    

  fGenLQp = fGenLQn = 0; 
  const int LQID(9000006); 
  // -- find LQ
  GenParticle *pGen(0);
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    pGen = getParticle(i); 
    if (LQID == TMath::Abs(pGen->PID)) {
      if (pGen->PID < 0) {
	fGenLQn = pGen; 
	genLQProducts(fGenLQn, fGenLQnL, fGenLQnQ, fGenLQnJ); 
      } else {
	fGenLQp = pGen; 
	genLQProducts(fGenLQp, fGenLQpL, fGenLQpQ, fGenLQpJ); 
      }
    }
  }
}



// ----------------------------------------------------------------------
void anaLq::fillHist(int cat) {

  fillRedTreeData(); 
  fTree->Fill();

  string histdir("");
  //  histdir = Form("class%d", cat); 
  //  cout << histdir << endl;
 

}



// ----------------------------------------------------------------------
void anaLq::bookHist() {
  fpHistFile->cd();			   
  cout << "==> anaLq: bookHist"  << endl;

  TH1D *h1(0); 
  TH2D *h2(0); 
  TProfile *hp(0);

  //  fpHistFile->mkdir(Form("class%d", i)); 
  //  fpHistFile->cd(Form("class%d", i));

  h1 = new TH1D("pt",  "lq gen pt", PTN, 0., PTMAX); 
  h1 = new TH1D("eta", "lq gen eta", PTN, -3., 3.); 
  h1 = new TH1D("phi", "lq gen phi", PTN, -3.15, 3.15); 
  h1 = new TH1D("m",   "lq gen m", PTN, 400., 1500.); 

  h1 = new TH1D("lqpt",  "l+q gen pt", PTN, 0., PTMAX); 
  h1 = new TH1D("lqeta", "l+q gen eta", PTN, -3., 3.); 
  h1 = new TH1D("lqphi", "l+q gen phi", PTN, -3.15, 3.15); 
  h1 = new TH1D("lqm",   "l+q gen m", PTN, 400., 1500.); 

  h1 = new TH1D("ljpt",  "l+j gen pt", PTN, 0., PTMAX); 
  h1 = new TH1D("ljeta", "l+j gen eta", PTN, -3., 3.); 
  h1 = new TH1D("ljphi", "l+j gen phi", PTN, -3.15, 3.15); 
  h1 = new TH1D("ljm",   "l+j gen m", PTN, 400., 1500.); 

  
  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupReducedTree();
  
  
}


// ----------------------------------------------------------------------
void anaLq::initVariables() {
  fW8 = 0.;
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
void anaLq::dumpGenBlock(bool withGluons) {

  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
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
      sl = ""; 
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
void anaLq::setupReducedTree() {
  fTree->Branch("type",    &fRtd.type,            "type/I");
  fTree->Branch("w8",      &fRtd.w8,              "w8/D");

  fTree->Branch("gpm",     &fRtd.gpm,             "gpm/D");
  fTree->Branch("gpm2",    &fRtd.gpm2,            "gpm2/D");
  fTree->Branch("gppt",    &fRtd.gppt,            "gppt/D");
  fTree->Branch("gpeta",   &fRtd.gpeta,           "gpeta/D");

  fTree->Branch("gnm",     &fRtd.gnm,             "gnm/D");
  fTree->Branch("gnm2",    &fRtd.gnm2,            "gnm2/D");
  fTree->Branch("gnpt",    &fRtd.gnpt,            "gnpt/D");
  fTree->Branch("gneta",   &fRtd.gneta,           "gneta/D");

  fTree->Branch("glqpm",   &fRtd.glqpm,           "glqpm/D");
  fTree->Branch("gljpm",   &fRtd.gljpm,           "gljpm/D");
  fTree->Branch("ljpm",    &fRtd.ljpm,            "ljpm/D");

  fTree->Branch("glqnm",   &fRtd.glqnm,           "glqnm/D");
  fTree->Branch("gljnm",   &fRtd.gljnm,           "gljnm/D");
  fTree->Branch("ljnm",    &fRtd.ljnm,            "ljnm/D");

}


// ----------------------------------------------------------------------
void anaLq::fillRedTreeData(int type) {

  fRtd.type = fClass; 
  fRtd.w8   = fW8;

  if (fGenLQp) {
    fRtd.gpm   = fGenLQp->Mass;
    fRtd.gpm2  = fP4GenLQp.M();
    fRtd.gppt  = fGenLQp->PT;
    fRtd.gpeta = fGenLQp->Eta;

    fRtd.glqpm = fP4GenLQpLQ.M();  
    fRtd.gljpm = fP4GenLQpLJ.M();  
    fRtd.ljpm;  

  } else {
    fRtd.gpm   = -9999.;
    fRtd.gpm2  = -9999.;
    fRtd.gppt  = -9999.;
    fRtd.gpeta = -9999.;
  }




  if (fGenLQn) {
    fRtd.gnm   = fGenLQn->Mass;
    fRtd.gnm2  = fP4GenLQn.M();
    fRtd.gnpt  = fGenLQn->PT;
    fRtd.gneta = fGenLQn->Eta;

    fRtd.glqnm = fP4GenLQnLQ.M();  
    fRtd.gljnm = fP4GenLQnLJ.M();  
    fRtd.ljnm;  

  } else {
    fRtd.gnm   = -9999.;
    fRtd.gnm2  = -9999.;
    fRtd.gnpt  = -9999.;
    fRtd.gneta = -9999.;
  }


  //  cout << "  fillRedTreeData: " << fW8 << endl;

}


// ----------------------------------------------------------------------
void anaLq::genLQProducts(GenParticle *lq, GenParticle *l, GenParticle *q, Jet *j) {

  l = q = 0; 
  j = 0; 
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
    cout << "LQ found, but decay quark or lepton missing?!" << endl;
    return;
  }

  TLorentzVector p4Q; 
  p4Q.SetPtEtaPhiM(q->PT, q->Eta, q->Phi, q->Mass); 

  //   printParticle(lq); 
  //   printParticle(l); 
  //   printParticle(q); 

  // -- find gen-jet for quark
  Jet *pJet(0); 
  double dRMin(9999.), dPtMin(9999.);
  int dRBest(-1), dPtBest(-1); 
  double dR(0.), dPt(0.);
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
    dPt = TMath::Abs(p4J.Pt() - p4Q.Pt());
    if (dR < dRMin) {
      dRMin = dR; 
      dRBest = i; 
    }
    if (dPt < dPtMin) {
      dPtMin = dPt; 
      dPtBest = i; 
    }
  }    

  // -- so far choose on best delta(R) match!
  if (dRBest > -1 && dRMin < 0.3) {
    j = getGenJet(dRBest);
  }

  // -- now fill data members, chosing according to charge of lepton
  if (l->Charge < 0) {
    fGenLQn  = lq; 
    fP4GenLQn.SetPtEtaPhiM(lq->PT, lq->Eta, lq->Phi, lq->Mass);

    fGenLQnL = l;
    fP4GenLQnL.SetPtEtaPhiM(l->PT, l->Eta,  l->Phi,  l->Mass);

    fGenLQnQ = q; 
    fP4GenLQnQ.SetPtEtaPhiM(q->PT, q->Eta,  q->Phi,  q->Mass);

    fP4GenLQnLQ = fP4GenLQnL + fP4GenLQnQ;    

    if (0 != j) {
      fGenLQnJ = j;
      fP4GenLQnJ.SetPtEtaPhiM(j->PT, j->Eta,  j->Phi,  j->Mass);
      fP4GenLQnLJ = fP4GenLQnL + fP4GenLQnJ;

      if (0 && fP4GenLQnLJ.M() < 300) {
	printSummary(3); 
	cout << "dRMin = " << dRMin << " gen jet " << dRBest 
	     << " and q " << genIndex(q) 
	     << " and l " << genIndex(l) 
	     << endl;
      }
    } else {
      fGenLQnJ = 0; 
      fP4GenLQnJ.SetPtEtaPhiM(-9999., -9999., -9999., -9999.);
      fP4GenLQnLJ.SetPtEtaPhiM(-9999., -9999., -9999., -9999.);
    }
  } else {
    fGenLQp  = lq; 
    fP4GenLQp.SetPtEtaPhiM(lq->PT, lq->Eta, lq->Phi, lq->Mass);

    fGenLQpL = l;
    fP4GenLQpL.SetPtEtaPhiM(l->PT, l->Eta,  l->Phi,  l->Mass);

    fGenLQpQ = q; 
    fP4GenLQpQ.SetPtEtaPhiM(q->PT, q->Eta,  q->Phi,  q->Mass);

    fP4GenLQpLQ = fP4GenLQpL + fP4GenLQpQ;    

    if (0 != j) {
      fGenLQpJ = j;
      fP4GenLQpJ.SetPtEtaPhiM(j->PT, j->Eta,  j->Phi,  j->Mass);
      fP4GenLQpLJ = fP4GenLQpL + fP4GenLQpJ;
      if (0 && fP4GenLQpLJ.M() < 300) { 
	printSummary(3); 
	cout << "dRMin = " << dRMin << " gen jet " << dRBest 
	     << " and q " << genIndex(q) 
	     << " and l " << genIndex(l) 
	     << endl;
      }
    } else {
      fGenLQpJ = 0; 
      fP4GenLQpJ.SetPtEtaPhiM(-9999., -9999., -9999., -9999.);
      fP4GenLQpLJ.SetPtEtaPhiM(-9999., -9999., -9999., -9999.);
    }
  }
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

