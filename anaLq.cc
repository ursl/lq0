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

  for (unsigned int i = 0; i < fGenLQ.size(); ++i) delete fGenLQ[i];
  fGenLQ.clear();


  if (0) {
    cout << "======================================================================" << endl;
    dumpGenBlock(true); 
    cout << "----------------------------------------------------------------------" << endl;
  }


  static int first(1); 
  if (first) {
    first = 0; 
    for (int i = 0; i < 15; ++i) {
      printParticle(getParticle(i)); 
    }
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

  if (pGen1) {
    genLQProducts(pGen0);
    genLQProducts(pGen1);
  } else {
    genLQSingle(pGen0);
  }
}


// ----------------------------------------------------------------------
void anaLq::analysis() {
  leptonSelection();
  jetSelection();
  preselection();

  if (fPreselected) {
    if (2 == TYPE) lqlqSelection();
    if (1 == TYPE) lqSelection();
    candAnalysis();
  }
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
      l->p4        = m4;
      l->q         = pM->Charge;
      l->pMuon     = pM;
      l->pElectron = 0; 
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
  fGoodCandLQp = false; 
  fGoodCandLQn = false; 

  fPreselected = false;
  fST = -9999.; 

  if (2 == TYPE) {
    if (fLeptons.size() < 2) return;
    if (fJets.size() < 2)    return;

    if (fLeptons[0]->p4.Pt() < L0PT) return;
    if (fJets[0]->p4.Pt() < J0PT)    return;

    if (fLeptons[0]->p4.DeltaR(fLeptons[1]->p4) < 0.3) return;

    TLorentzVector ll4; 
    ll4 = fLeptons[0]->p4 + fLeptons[1]->p4;
    if (ll4.M() < 50.)  return;

    fST = fLeptons[0]->p4.Pt() + fLeptons[1]->p4.Pt() + fJets[0]->p4.Pt() + fJets[1]->p4.Pt();
    if (fST < 300) return;

    fPreselected = true;
    return;
  }


  if (1 == TYPE) {
    if (fLeptons.size() < 1) return;
    if (fJets.size() < 1)    return;

    if (fLeptons[0]->p4.Pt() < L0PT) return;
    if (fJets[0]->p4.Pt() < J0PT)    return;

    fST = fLeptons[0]->p4.Pt() + fJets[0]->p4.Pt();
    fPreselected = true;
    return;
  }

}

// ----------------------------------------------------------------------
void anaLq::lqlqSelection() {

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
  fLQ.push_back(Lq); 

  Lq = new lq;
  Lq->p4   = lq1; 
  Lq->q    = fLeptons[1]->q; 
  Lq->idxL = 1; 
  Lq->idxK = -1; 
  Lq->idxJ = 1-iBest; 
  fLQ.push_back(Lq); 
      
}

// ----------------------------------------------------------------------
void anaLq::lqSelection() {
  // FIXME this has to be improved!
  for (unsigned int i = 0; i < fLQ.size(); ++i) delete fLQ[i];
  fLQ.clear();
  
  TLorentzVector lq0 = fLeptons[0]->p4 + fJets[0]->p4; 

  lq *Lq = new lq;
  Lq->p4 = lq0; 
  Lq->q    = fLeptons[0]->q; 
  Lq->idxL = 0; 
  Lq->idxK = (fLeptons.size() > 1? 1:-1); 
  Lq->idxJ = 0;
  fLQ.push_back(Lq); 

}


// ----------------------------------------------------------------------
// -- this is the identical selection as in plotLq.cc:candAnalysis, coded w/o the reduced tree for consistency checks!
void anaLq::candAnalysis() {

  fGoodEvent   = false; 
  
  if (2 == TYPE) {
    if (2 == fLQ.size() && fLQ[0]->p4.M() > 0. && fLQ[1]->p4.M() > 0.) fGoodEvent = true;
    if (fST < 685.)       fGoodEvent = false;
    if (fMll < 150.)      fGoodEvent = false;
    if (fMljMin < 155.)   fGoodEvent = false;
  } else {
    if (1 == fLQ.size() && fLQ[0]->p4.M() > 0.) fGoodEvent = true;
  }

}


// ----------------------------------------------------------------------
void anaLq::fillHist() {

  fillRedTreeData(); 
  fTree->Fill();

  char cds[200]; 
  sprintf(cds, "%s", fName.c_str());

  if (fPreselected) { 
    fHists[Form("pre_l0pt_%s", cds)]->Fill(fLeptons[0]->p4.Pt()); 
    if (2 == TYPE) fHists[Form("pre_l1pt_%s", cds)]->Fill(fLeptons[1]->p4.Pt()); 

    fHists[Form("pre_j0pt_%s", cds)]->Fill(fJets[0]->p4.Pt()); 
    if (2 == TYPE) fHists[Form("pre_j1pt_%s", cds)]->Fill(fJets[1]->p4.Pt()); 

    fHists[Form("pre_st_%s", cds)]->Fill(fRtd.st); 
    fHists[Form("pre_mll_%s", cds)]->Fill(fMll); 
    fHists[Form("pre_mljetmin_%s", cds)]->Fill(fMljMin); 
      
    if (fLQ[0]->p4.M() > 0.) {
      fHists[Form("pre_m_%s", cds)]->Fill(fLQ[0]->p4.M()); 
      fHists[Form("pre_pt_%s", cds)]->Fill(fLQ[0]->p4.Pt()); 
    }
    if (2 == fLQ.size() && fLQ[1]->p4.M() > 0.) {
      fHists[Form("pre_m_%s", cds)]->Fill(fLQ[1]->p4.M()); 
      fHists[Form("pre_pt_%s", cds)]->Fill(fLQ[1]->p4.Pt()); 
    }


    if (fGoodEvent) { 
      fHists[Form("sel_st_%s", cds)]->Fill(fRtd.st); 
      fHists[Form("sel_mll_%s", cds)]->Fill(fMll); 
      fHists[Form("sel_mljetmin_%s", cds)]->Fill(fMljMin); 
      
      if (fLQ[0]->p4.M() > 0.) {
	fHists[Form("sel_m_%s", cds)]->Fill(fLQ[0]->p4.M()); 
	fHists[Form("sel_pt_%s", cds)]->Fill(fLQ[0]->p4.Pt()); 
      }
      if (2 == fLQ.size() && fLQ[1]->p4.M() > 0.) {
	fHists[Form("sel_m_%s", cds)]->Fill(fLQ[1]->p4.M()); 
	fHists[Form("sel_pt_%s", cds)]->Fill(fLQ[1]->p4.Pt()); 
      }
    }

  }

}



// ----------------------------------------------------------------------
void anaLq::bookHist() {
  fpHistFile->cd();			   
  cout << "==> anaLq: bookHist"  << endl;

  TH1D *h1(0); 
  (void)h1;

  //  fpHistFile->mkdir(Form("class%d", i)); 
  //  fpHistFile->cd(Form("class%d", i));

  h1 = new TH1D("pt",  "lq gen pt", 40, 0., 400.); 
  h1 = new TH1D("phi", "lq gen phi", 40, -3.15, 3.15); 
  h1 = new TH1D("m",   "lq gen m", 22, 400., 1500.); 

  h1 = new TH1D("lqpt",  "l+q gen pt", 40, 0., 400); 
  h1 = new TH1D("lqphi", "l+q gen phi", 40, -3.15, 3.15); 
  h1 = new TH1D("lqm",   "l+q gen m", 22, 400., 1500.); 

  h1 = new TH1D("ljpt",  "l+j gen pt", 40, 0., 400); 
  h1 = new TH1D("ljphi", "l+j gen phi", 40, -3.15, 3.15); 
  h1 = new TH1D("ljm",   "l+j gen m", 40, 400., 1500.); 


  vector<string> levels; 
  levels.push_back("pre"); 
  levels.push_back("sel"); 

  for (unsigned int i = 0; i < levels.size(); ++i) {
  
    // -- histograms as in plotLq:

    // -- m
    fHists.insert(make_pair(Form("%s_m_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_m_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_m_%s", levels[i].c_str(), fName.c_str()), 
				     40, 0, 2000.))); 
    setTitles(fHists[Form("%s_m_%s", levels[i].c_str(), fName.c_str())], "m [GeV]", "Entries/bin");
    
    // -- st
    fHists.insert(make_pair(Form("%s_st_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_st_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_st_%s", levels[i].c_str(), fName.c_str()), 
				     14, 0, 3500.))); 
    setTitles(fHists[Form("%s_st_%s", levels[i].c_str(), fName.c_str())], "S_{T} [GeV]", "Entries/bin");
    
    // -- mll
    fHists.insert(make_pair(Form("%s_mll_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_mll_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_mll_%s", levels[i].c_str(), fName.c_str()), 
				     30, 0, 1500.))); 
    setTitles(fHists[Form("%s_mll_%s", levels[i].c_str(), fName.c_str())], "m_{l l} [GeV]", "Entries/bin");
    
    // -- mljetmin
    fHists.insert(make_pair(Form("%s_mljetmin_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_mljetmin_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_mljetmin_%s", levels[i].c_str(), fName.c_str()),
				     15, 0, 1500.))); 
    setTitles(fHists[Form("%s_mljetmin_%s", levels[i].c_str(), fName.c_str())], "m_{l jet}^{min} [GeV]", "Entries/bin");
    
    
    // -- pt
    fHists.insert(make_pair(Form("%s_pt_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_pt_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_pt_%s", levels[i].c_str(), fName.c_str()), 
				     100, 0, 1000.))); 
    setTitles(fHists[Form("%s_pt_%s", levels[i].c_str(), fName.c_str())], "p_{T} [GeV]", "Entries/bin");

    // -- lepton pt
    fHists.insert(make_pair(Form("%s_l0pt_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_l0pt_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_l0pt_%s", levels[i].c_str(), fName.c_str()), 
				     64, 0, 1600.))); 
    setTitles(fHists[Form("%s_l0pt_%s", levels[i].c_str(), fName.c_str())], "p_{T}(l_{0}) [GeV]", "Entries/bin");

    fHists.insert(make_pair(Form("%s_l1pt_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_l1pt_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_l1pt_%s", levels[i].c_str(), fName.c_str()), 
				     40, 0, 800.))); 
    setTitles(fHists[Form("%s_l1pt_%s", levels[i].c_str(), fName.c_str())], "p_{T}(l_{1}) [GeV]", "Entries/bin");

    // -- jet pt
    fHists.insert(make_pair(Form("%s_j0pt_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_j0pt_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_j0pt_%s", levels[i].c_str(), fName.c_str()), 
				     64, 0, 1600.))); 
    setTitles(fHists[Form("%s_j0pt_%s", levels[i].c_str(), fName.c_str())], "p_{T}(j_{0}) [GeV]", "Entries/bin");

    fHists.insert(make_pair(Form("%s_j1pt_%s", levels[i].c_str(), fName.c_str()), 
			    new TH1D(Form("%s_j1pt_%s", levels[i].c_str(), fName.c_str()), 
				     Form("%s_j1pt_%s", levels[i].c_str(), fName.c_str()), 
				     40, 0, 800.))); 
    setTitles(fHists[Form("%s_j1pt_%s", levels[i].c_str(), fName.c_str())], "p_{T}(j_{1}) [GeV]", "Entries/bin");
  }
  
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

  fTree->Branch("ngen",    &fRtd.ngen,            "ngen/I");
  fTree->Branch("gm",       fRtd.gm,              "gm[ngen]/D");
  fTree->Branch("gpt",      fRtd.gpt,             "gpt[ngen]/D");
  fTree->Branch("geta",     fRtd.geta,            "geta[ngen]/D");
  fTree->Branch("gphi",     fRtd.gphi,            "gphi[ngen]/D");
  fTree->Branch("gmlj",     fRtd.gmlj,            "gmlj[ngen]/D");
  fTree->Branch("glq",      fRtd.glq,             "glq[ngen]/I");

  fTree->Branch("glpt",     &fRtd.glpt,           "glpt[ngen]/D");
  fTree->Branch("gleta",    &fRtd.gleta,          "gleta[ngen]/D");
  fTree->Branch("glphi",    &fRtd.glphi,          "glphi[ngen]/D");

  fTree->Branch("gqpt",     &fRtd.gqpt,           "gqpt[ngen]/D");
  fTree->Branch("gqeta",    &fRtd.gqeta,          "gqeta[ngen]/D");
  fTree->Branch("gqphi",    &fRtd.gqphi,          "gqphi[ngen]/D");

  fTree->Branch("gjpt",     &fRtd.gjpt,           "gjpt[ngen]/D");
  fTree->Branch("gjeta",    &fRtd.gjeta,          "gjeta[ngen]/D");
  fTree->Branch("gjphi",    &fRtd.gjphi,          "gjphi[ngen]/D");

  fTree->Branch("gkpt",     &fRtd.gkpt,           "gkpt[ngen]/D");
  fTree->Branch("gketa",    &fRtd.gketa,          "gketa[ngen]/D");
  fTree->Branch("gkphi",    &fRtd.gkphi,          "gkphi[ngen]/D");

  fTree->Branch("nrec",    &fRtd.nrec,            "nrec/I");
  fTree->Branch("m",       &fRtd.m,               "ljm[nrec]/D");
  fTree->Branch("pt",      &fRtd.pt,              "pt[nrec]/D");
  fTree->Branch("eta" ,    &fRtd.eta,             "eta[nrec]/D");
  fTree->Branch("phi",     &fRtd.phi,             "phi[nrec]/D");
  fTree->Branch("q",       &fRtd.lq,              "lq[nrec]/I");

  fTree->Branch("lpt",      &fRtd.lpt,              "lpt[nrec]/D");
  fTree->Branch("leta" ,    &fRtd.leta,             "leta[nrec]/D");
  fTree->Branch("lphi",     &fRtd.lphi,             "lphi[nrec]/D");

  fTree->Branch("jpt",      &fRtd.jpt,              "jpt[nrec]/D");
  fTree->Branch("jeta" ,    &fRtd.jeta,             "jeta[nrec]/D");
  fTree->Branch("jphi",     &fRtd.jphi,             "jphi[nrec]/D");

  fTree->Branch("kpt",      &fRtd.kpt,              "kpt[nrec]/D");
  fTree->Branch("keta" ,    &fRtd.keta,             "keta[nrec]/D");
  fTree->Branch("kphi",     &fRtd.kphi,             "kphi[nrec]/D");

  fTree->Branch("st",       &fRtd.st,              "st/D");
  fTree->Branch("mll" ,     &fRtd.mll,             "mll/D");
  fTree->Branch("mljmin",   &fRtd.mljmin,          "mljmin/D");



}


// ----------------------------------------------------------------------
void anaLq::fillRedTreeData() {

  fRtd.type    = TYPE; 
  fRtd.channel = CHANNEL; 
  fRtd.w8      = fW8;

  for (unsigned int i = 0; i < fGenLQ.size(); ++i) {
    fRtd.gm[i]   = fGenLQ[i]->p4LQ.M(); 
    fRtd.gpt[i]  = fGenLQ[i]->p4LQ.Pt(); 
    fRtd.geta[i] = fGenLQ[i]->p4LQ.Eta(); 
    fRtd.gphi[i] = fGenLQ[i]->p4LQ.Phi(); 

    fRtd.glq[i]  = fGenLQ[i]->q; 

    fRtd.glpt[i]  = fGenLQ[i]->p4L.Pt(); 
    fRtd.gleta[i] = fGenLQ[i]->p4L.Eta(); 
    fRtd.glphi[i] = fGenLQ[i]->p4L.Phi(); 

    fRtd.gqpt[i]  = fGenLQ[i]->p4Q.Pt(); 
    fRtd.gqeta[i] = fGenLQ[i]->p4Q.Eta(); 
    fRtd.gqphi[i] = fGenLQ[i]->p4Q.Phi(); 

    if (fGenLQ[i]->pJ) {
      fRtd.gmlj[i]  = fGenLQ[i]->p4LJ.M(); 
      fRtd.gjpt[i]  = fGenLQ[i]->p4J.Pt(); 
      fRtd.gjeta[i] = fGenLQ[i]->p4J.Eta(); 
      fRtd.gjphi[i] = fGenLQ[i]->p4J.Phi(); 
    } else {
      fRtd.gmlj[i]  = -9999.;
    }

    if (fGenLQ[i]->pK) {
      fRtd.gkpt[i]  = fGenLQ[i]->p4K.Pt(); 
      fRtd.gketa[i] = fGenLQ[i]->p4K.Eta(); 
      fRtd.gkphi[i] = fGenLQ[i]->p4K.Phi(); 
    } else {
      fRtd.gkpt[i]  = -9999.;
      fRtd.gketa[i] = -9999.;
      fRtd.gkphi[i] = -9999.;
    }
  }
  fRtd.ngen = fGenLQ.size(); 

  if (fPreselected) {
    for (unsigned int i = 0; i < fLQ.size(); ++i) {
      fRtd.m[i]    = fLQ[i]->p4.M();
      fRtd.pt[i]   = fLQ[i]->p4.Pt();
      fRtd.eta[i]  = fLQ[i]->p4.Eta();
      fRtd.phi[i]  = fLQ[i]->p4.Phi();

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
	fRtd.kpt[i]   = -9999.;
	fRtd.keta[i]  = -9999.;
	fRtd.kphi[i]  = -9999.;
      }
    }
  }
  fRtd.nrec = fLQ.size(); 


}


// ----------------------------------------------------------------------
void anaLq::genLQProducts(GenParticle *lq) {

  genLq *gen = new genLq(); 

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
  if (dRBest > -1 && dRMin < 0.3) {
    j = getGenJet(dRBest);
    //    cout << " xx choosing jet = " << dRBest << " with ET = " << j->PT << " for jet of LQ " << lqIdx << endl;
  }

  // -- fill struct
  gen->q = l->Charge;

  gen->pLQ = lq;
  gen->p4LQ.SetPtEtaPhiM(lq->PT, lq->Eta, lq->Phi, lq->Mass);
  
  gen->pL = l; 
  gen->p4L.SetPtEtaPhiM(l->PT, l->Eta, l->Phi, l->Mass);
  
  gen->pQ = q; 
  gen->p4Q.SetPtEtaPhiM(q->PT, q->Eta, q->Phi, q->Mass);

  if (j) {
    gen->pJ = j;
    gen->p4J.SetPtEtaPhiM(j->PT, j->Eta,  j->Phi,  j->Mass);
    
    TLorentzVector lj = gen->p4L + gen->p4J; 
    gen->p4LJ = lj; 
  } else {
    gen->pJ = 0;
    gen->p4J.SetPtEtaPhiM(-9999., -9999., -9999., -9999.); 
    gen->p4LJ.SetPtEtaPhiM(-9999., -9999., -9999., -9999.); 
  }

  // -- this will be reset in genLqSingle for single LQ processing
  gen->pK = 0;
  gen->p4K.SetPtEtaPhiM(-9999., -9999., -9999., -9999.); 

  fGenLQ.push_back(gen);
}


// ----------------------------------------------------------------------
void anaLq::genLQSingle(GenParticle *lq) {

  // -- Fill LQ into single variables
  genLQProducts(lq); 

  // -- find gen-level bachelor lepton K to LQ
  int lqIdx = genIndex(lq);
  int lqMomIdx = lq->M1;
  GenParticle *pGen(0);
  for (int i = 0; i < fbParticles->GetEntries(); ++i) {
    if (i == lqIdx) continue;
    pGen = getParticle(i); 
    if (pGen->M1 == lqMomIdx) {
      break;
    }
  }

  genLq *gen = fGenLQ.back(); 
  gen->pK = pGen;
  gen->p4K.SetPtEtaPhiM(pGen->PT, pGen->Eta, pGen->Phi, pGen->Mass);
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
