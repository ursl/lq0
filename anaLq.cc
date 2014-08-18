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
  return ( a->fP4.Pt() > b->fP4.Pt() );
}

bool sortPtJet(jet *a, jet *b) {
  return ( a->fP4.Pt() > b->fP4.Pt() );
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
    
  }

}



// ----------------------------------------------------------------------
void anaLq::startAnalysis() {
  cout << "==> anaLq: startAnalysis: " << (TYPE==2?"LQ ***pair*** production":"LQ ***single*** production") 
       << " in the ***" << (CHANNEL==13?"muon":"electron") << "*** channel"
       << endl;

  MUISODELTAR = 0.3;

  L0PT = 45.; 
  L1PT = 45.; 

  J0PT = 125.; 
  J1PT = 45.; 
}

// ----------------------------------------------------------------------
void anaLq::endAnalysis() {
  cout << "==> anaLq: endAnalysis: ..." << endl;
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
  analysis();
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
void anaLq::analysis() {
  leptonSelection();
  jetSelection();
  preselection();

  if (!fPreselected) return;

  if (2 == TYPE) lqlqSelection();
  if (1 == TYPE) lqSelection();
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
      if (TMath::Abs(pM->Eta) > 2.1) continue;
      if (muonIso(pM) > 0.1) continue;
      lepton *l = new lepton;
      m4.SetPtEtaPhiM(pM->PT, pM->Eta, pM->Phi, MMUON); 
      l->fP4 = m4;
      l->q   = pM->Charge;
      l->fpMuon = pM;
      l->fpElectron = 0; 
      fLeptons.push_back(l);
    }      

    if (fLeptons.size() > 1) {
      sort(fLeptons.begin(), fLeptons.end(), sortPtLepton);
    }

  }

  if (11 == CHANNEL) {
    cout << "implement electron selection!!" << endl;
    exit(1); 
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
    j->fP4 = j4;
    j->fpJet = pJ; 
    fJets.push_back(j);
  }
  
  if (fJets.size() > 1) {
    sort(fJets.begin(), fJets.end(), sortPtJet);
  }

}

// ----------------------------------------------------------------------
void anaLq::preselection() {

  fPreselected = false;
  fST = -9999.; 

  if (2 == TYPE) {
    if (fLeptons.size() < 2) return;
    if (fJets.size() < 2)    return;

    if (fLeptons[0]->fP4.Pt() < L0PT) return;
    if (fJets[0]->fP4.Pt() < J0PT)    return;

    TLorentzVector ll4; 
    ll4 = fLeptons[0]->fP4 + fLeptons[1]->fP4;
    if (ll4.M() < 50.)  return;

    fST = fLeptons[0]->fP4.Pt() + fLeptons[1]->fP4.Pt() + fJets[0]->fP4.Pt() + fJets[1]->fP4.Pt();
    if (fST < 300) return;

    fPreselected = true;
    return;
  }


  if (1 == TYPE) {
    if (fLeptons.size() < 1) return;
    if (fJets.size() < 1)    return;

    if (fLeptons[0]->fP4.Pt() < L0PT) return;
    if (fJets[0]->fP4.Pt() < J0PT)    return;

    fST = fLeptons[0]->fP4.Pt() + fJets[0]->fP4.Pt();
    fPreselected = true;
    return;
  }

}

// ----------------------------------------------------------------------
void anaLq::lqlqSelection() {

  for (unsigned int i = 0; i < fLQ.size(); ++i) delete fLQ[i];
  fLQ.clear();

  
  if (!fPreselected) return;

  int iBest(-1);
  double mdiff(99999.), mdiffBest(99999.); 
  TLorentzVector lq0, lq1; 
  for (int i = 0; i < 2; ++i) {
    lq0 = fLeptons[0]->fP4 + fJets[i]->fP4; 
    lq1 = fLeptons[1]->fP4 + fJets[1-i]->fP4; 
    //    cout << "i = " << i << " LQ0: " << lq0.M() << " LQ1: " << lq1.M() << endl;
    mdiff = TMath::Abs(lq0.M() - lq1.M()); 
    if (mdiff < mdiffBest) {
      mdiffBest = mdiff; 
      iBest = i; 
    }
  }

  lq0 = fLeptons[0]->fP4 + fJets[iBest]->fP4; 
  lq1 = fLeptons[1]->fP4 + fJets[1-iBest]->fP4; 
  //  cout << "iBest = " << iBest << " LQ0 = " << lq0.M() << "  LQ1 = " << lq1.M() << endl;

  fMll = (fLeptons[0]->fP4 + fLeptons[1]->fP4).M();

  // -- I hope the following is correct
  double mlq0 = lq0.M();
  double mlq1 = lq1.M();
  if (mlq0 < mlq1) {
    fMljetMin = mlq0;
  } else {
    fMljetMin = mlq1;
  }  

  // -- charge determination
  fPos = (fLeptons[0]->q > 0?0:1);  
  fNeg = (fLeptons[1]->q > 0?0:1);  
  if (fPos == fNeg) {
    fNeg = 1 - fPos; 
    cout << "XXXXXXXX problem with charge assignment, using positive muon only, i.e. pos = " << fPos << " neg = " << fNeg << endl;
  }

  // -- leading lepton/jet determination
  fL0 = (fLeptons[0]->fP4.Pt() > fLeptons[1]->fP4.Pt()? 0: 1);  
  fL1 = (fLeptons[0]->fP4.Pt() > fLeptons[1]->fP4.Pt()? 1: 0);  

  fJ0 = (fJets[0]->fP4.Pt() > fJets[1]->fP4.Pt()? 0: 1);  
  fJ1 = (fJets[0]->fP4.Pt() > fJets[1]->fP4.Pt()? 1: 0);  


  lq *Lq = new lq;
  Lq->fP4 = lq0; 
  fLQ.push_back(Lq); 

  Lq = new lq;
  Lq->fP4 = lq1; 
  fLQ.push_back(Lq); 
      
}

// ----------------------------------------------------------------------
void anaLq::lqSelection() {
  
  for (unsigned int i = 0; i < fLQ.size(); ++i) delete fLQ[i];
  fLQ.clear();
  
  if (!fPreselected) return;

  TLorentzVector lq0 = fLeptons[0]->fP4 + fJets[0]->fP4; 

  // -- charge determination
  fPos = (fLeptons[0]->q > 0?0:-1);  
  fNeg = (fLeptons[0]->q < 0?0:-1);  

  lq *Lq = new lq;
  Lq->fP4 = lq0; 
  fLQ.push_back(Lq); 
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
  (void)h1;

  //  fpHistFile->mkdir(Form("class%d", i)); 
  //  fpHistFile->cd(Form("class%d", i));

  h1 = new TH1D("pt",  "lq gen pt", PTN, 0., PTMAX); 
  h1 = new TH1D("phi", "lq gen phi", PTN, -3.15, 3.15); 
  h1 = new TH1D("m",   "lq gen m", PTN, 400., 1500.); 

  h1 = new TH1D("lqpt",  "l+q gen pt", PTN, 0., PTMAX); 
  h1 = new TH1D("lqphi", "l+q gen phi", PTN, -3.15, 3.15); 
  h1 = new TH1D("lqm",   "l+q gen m", PTN, 400., 1500.); 

  h1 = new TH1D("ljpt",  "l+j gen pt", PTN, 0., PTMAX); 
  h1 = new TH1D("ljphi", "l+j gen phi", PTN, -3.15, 3.15); 
  h1 = new TH1D("ljm",   "l+j gen m", PTN, 400., 1500.); 

  
  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupReducedTree();
  
  
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

  fTree->Branch("gnm",     &fRtd.gnm,             "gnm/D");
  fTree->Branch("gnm2",    &fRtd.gnm2,            "gnm2/D");
  fTree->Branch("gnpt",    &fRtd.gnpt,            "gnpt/D");

  fTree->Branch("glqpm",   &fRtd.glqpm,           "glqpm/D");
  fTree->Branch("gljpm",   &fRtd.gljpm,           "gljpm/D");

  fTree->Branch("glqnm",   &fRtd.glqnm,           "glqnm/D");
  fTree->Branch("gljnm",   &fRtd.gljnm,           "gljnm/D");

  
  fTree->Branch("ljnm",    &fRtd.ljnm,            "ljnm/D");
  fTree->Branch("ljnpt",   &fRtd.ljnpt,           "ljnpt/D");

  fTree->Branch("ljpm",    &fRtd.ljpm,            "ljpm/D");
  fTree->Branch("ljppt",   &fRtd.ljppt,           "ljppt/D");

  fTree->Branch("l0pt",    &fRtd.l0pt,            "l0pt/D");
  fTree->Branch("l1pt",    &fRtd.l1pt,            "l1pt/D");

  fTree->Branch("j0pt",    &fRtd.j0pt,            "j0pt/D");
  fTree->Branch("j1pt",    &fRtd.j1pt,            "j1pt/D");


  fTree->Branch("st",      &fRtd.st,              "st/D");
  fTree->Branch("mll",     &fRtd.mll,             "mll/D");
  fTree->Branch("mljetmin",&fRtd.mljetmin,        "mljetmin/D");


}


// ----------------------------------------------------------------------
void anaLq::fillRedTreeData(int type) {

  fRtd.type    = TYPE; 
  fRtd.channel = CHANNEL; 
  fRtd.w8      = fW8;

  if (fGenLQp) {
    fRtd.gpm   = fGenLQp->Mass;
    fRtd.gpm2  = fP4GenLQp.M();
    fRtd.gppt  = fGenLQp->PT;

    fRtd.glqpm = fP4GenLQpLQ.M();  
    fRtd.gljpm = fP4GenLQpLJ.M();  

  } else {
    fRtd.gpm   = -9999.;
    fRtd.gpm2  = -9999.;
    fRtd.gppt  = -9999.;

    fRtd.glqpm = -9999.;
    fRtd.gljpm = -9999.;
  }

  if (fGenLQn) {
    fRtd.gnm   = fGenLQn->Mass;
    fRtd.gnm2  = fP4GenLQn.M();
    fRtd.gnpt  = fGenLQn->PT;

    fRtd.glqnm = fP4GenLQnLQ.M();  
    fRtd.gljnm = fP4GenLQnLJ.M();  

  } else {
    fRtd.gnm   = -9999.;
    fRtd.gnm2  = -9999.;
    fRtd.gnpt  = -9999.;

    fRtd.glqnm = -9999.;
    fRtd.gljnm = -9999.;
  }


  if (fPreselected) {
    if (fPos > -1) {
      fRtd.ljpm = fLQ[fPos]->fP4.M();  
      fRtd.ljppt = fLQ[fPos]->fP4.Pt();  
    } else {
      fRtd.ljpm = -9999.;
      fRtd.ljppt = -9999.;  
    }
    
    if (fNeg > -1) {
      fRtd.ljnm = fLQ[fNeg]->fP4.M();  
      fRtd.ljnpt = fLQ[fNeg]->fP4.Pt();  
    } else {
      fRtd.ljnm = -9999.;  
      fRtd.ljnpt = -9999.;  
    }
    
    fRtd.l0pt = fLeptons[fL0]->fP4.Pt();
    fRtd.j0pt = fJets[fL0]->fP4.Pt();

    if (2 == TYPE) {
      fRtd.l1pt = fLeptons[fL1]->fP4.Pt();
      fRtd.j1pt = fJets[fL1]->fP4.Pt();

      fRtd.mljetmin = fMljetMin; 
      fRtd.mll      = fMll; 
    } 
    fRtd.st       = fST; 
  } else {
    fRtd.ljpm = -9999.;
    fRtd.ljppt = -9999.;  

    fRtd.ljnm = -9999.;  
    fRtd.ljnpt = -9999.;  

    fRtd.l0pt = -9999.;  
    fRtd.l1pt = -9999.;  

    fRtd.j0pt = -9999.;  
    fRtd.j1pt = -9999.;  

    fRtd.mljetmin = -9999.; 
    fRtd.mll      = -9999.; 
    fRtd.st       = -9999.;
  }

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
    cout << "XXXXX LQ found, but decay quark or lepton missing?!" << endl;
    return;
  }

  TLorentzVector p4Q; 
  p4Q.SetPtEtaPhiM(q->PT, q->Eta, q->Phi, q->Mass); 

  //   printParticle(lq); 
  //   printParticle(l); 
  //   printParticle(q); 

  // -- find gen-jet for quark
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
