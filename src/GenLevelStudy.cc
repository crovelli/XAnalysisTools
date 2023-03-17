#include "../include/GenLevelStudy.hh"
#include "TLorentzVector.h"
#include <TCanvas.h>
#include <TStyle.h>

using namespace std;

GenLevelStudy::GenLevelStudy(TTree *tree)
  : XBase(tree) {

}

GenLevelStudy::~GenLevelStudy() {} 

bool GenLevelStudy::contains(vector<int> vec, const int & elem) {

  bool result = false;
  if( find(vec.begin(), vec.end(), elem) != vec.end() ) {
    result = true;
  }
  return result;
}

void GenLevelStudy::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // Mass plots
  TH1F *myHmassXrecoAll = new TH1F("myHmassXrecoAll", "myHmassXrecoAll", 100, 3.4, 4.1);

  // Mass plots separate per decay
  TH1F *myHmassXgen1  = new TH1F("myHmassXgen1",  "myHmassXgen1",  100, 2.5, 5.5);
  TH1F *myHmassXreco1 = new TH1F("myHmassXreco1", "myHmassXreco1", 100, 3.4, 4.1);
  TH1F *myHmassXgen2  = new TH1F("myHmassXgen2",  "myHmassXgen2",  100, 2.5, 5.5);
  TH1F *myHmassXreco2 = new TH1F("myHmassXreco2", "myHmassXreco2", 100, 3.4, 4.1);
  TH1F *myHmassXgen3a = new TH1F("myHmassXgen3a", "myHmassXgen3a", 100, 2.5, 5.5);
  TH1F *myHmassXgen3b = new TH1F("myHmassXgen3b", "myHmassXgen3b", 100, 3.4, 4.1);
  TH1F *myHmassXgen3c = new TH1F("myHmassXgen3c", "myHmassXgen3c", 100, 2.5, 5.5);
  TH1F *myHmassXreco3 = new TH1F("myHmassXreco3", "myHmassXreco3", 100, 3.4, 4.1);
  TH1F *myHmassXgen4a = new TH1F("myHmassXgen4a", "myHmassXgen4a", 100, 2.5, 5.5);
  TH1F *myHmassXgen4b = new TH1F("myHmassXgen4b", "myHmassXgen4b", 100, 3.4, 4.1);
  TH1F *myHmassXgen4c = new TH1F("myHmassXgen4c", "myHmassXgen4c", 100, 2.5, 5.5);
  TH1F *myHmassXreco4 = new TH1F("myHmassXreco4", "myHmassXreco4", 100, 3.4, 4.1);
  TH1F *myHmassXgen5  = new TH1F("myHmassXgen5",  "myHmassXgen5",  100, 2.5, 5.5);
  TH1F *myHmassXreco5 = new TH1F("myHmassXreco5", "myHmassXreco5", 100, 3.4, 4.1);
  TH1F *myHmassXgen6  = new TH1F("myHmassXgen6",  "myHmassXgen6",  100, 2.5, 5.5);
  TH1F *myHmassXreco6 = new TH1F("myHmassXreco6", "myHmassXreco6", 100, 3.4, 4.1);
  TH1F *myHmassXgen7a = new TH1F("myHmassXgen7a", "myHmassXgen7a", 100, 2.5, 5.5);
  TH1F *myHmassXgen7b = new TH1F("myHmassXgen7b", "myHmassXgen7b", 100, 3.4, 4.1);
  TH1F *myHmassXreco7 = new TH1F("myHmassXreco7", "myHmassXreco7", 100, 3.4, 4.1);
  
  // To check the JPsi mother
  TH1F *myHmothId_500    = new TH1F("myHmothId_500",    "myHmothId_500",       500,       0.,     500.);
  TH1F *myHmothId_520    = new TH1F("myHmothId_520",    "myHmothId_520",        20,     500.,     520.);
  TH1F *myHmothId_1000   = new TH1F("myHmothId_1000",   "myHmothId_1000",      480,     520.,    1000.);
  TH1F *myHmothId_10000  = new TH1F("myHmothId_10000",  "myHmothId_10000",    9000,    1000.,   10000.);
  TH1F *myHmothId_20000  = new TH1F("myHmothId_20000",  "myHmothId_20000",   10000,   10000.,   20000.);
  TH1F *myHmothId_100000 = new TH1F("myHmothId_100000", "myHmothId_100000",  10000,   20000.,  100000.);
  TH1F *myHmothId_110000 = new TH1F("myHmothId_110000", "myHmothId_110000",  10000,  100000.,  110000.);

  


  // To select a specific event
  int thisIsThe = 0;

  // Vectors with mother types
  vector<int> mothNum;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    thisIsThe++;
    
    // Events with >=1 reconstructed B with Bmass and K0 in the signal region
    // mX should be in the enlarged Psi2s - X range  
    int bestBmass = -1;
    float bestBmassDelta = -1;
    bool atLeastOneB = false;
    for (int nB=0; nB<nB0; nB++) {
      //if (B0_K0s_nmcFitted_mass[nB]<499.30000 && B0_K0s_nmcFitted_mass[nB]<495.50000) {
      //if (B0_finalFit_mass[nB]<5.2802900 && B0_finalFit_mass[nB]>5.2785100) {
      //if (B0_finalFit_X_mass[nB]<4.1 && B0_finalFit_X_mass[nB]>3.4) { 
      if (B0_K0s_nmcFitted_mass[nB]<0.541 && B0_K0s_nmcFitted_mass[nB]>0.454) {
	if (B0_finalFit_mass[nB]<5.35 && B0_finalFit_mass[nB]>5.20) {
	  if (B0_finalFit_X_mass[nB]<4.1 && B0_finalFit_X_mass[nB]>3.4) { // X: 3.81-3.93; Psi2s: 3.65-3.72
	    myHmassXrecoAll->Fill(B0_finalFit_X_mass[nB]);
	    atLeastOneB = true;
	    if (fabs(B0_finalFit_mass[nB]-5.27966)<bestBmassDelta) {
	      bestBmassDelta = fabs(B0_finalFit_mass[nB]-5.27966);
	      bestBmass = nB;
	    }
	  }}}
    }
    if (!atLeastOneB) continue;
    
    // Gen level: which is the jpsi mother
    int theBMothIdx = -1;
    int theChic0MothIdx = -1;                                                         
    int theChic1MothIdx = -1;
    int theChic2MothIdx = -1;                                                                                                                                                         
    for (int ii=0; ii<nGenPart; ii++){
      int pdgId = GenPart_pdgId[ii];
      int mothIdx = GenPart_genPartIdxMother[ii];
      if (mothIdx<0) continue;
      int mothId = GenPart_pdgId[mothIdx];

      if (pdgId==443 && abs(mothId)<500) myHmothId_500 -> Fill(mothId);
      if (pdgId==443 && abs(mothId)>=500 && abs(mothId)<520) myHmothId_520 -> Fill(mothId);
      if (pdgId==443 && abs(mothId)>=520 && abs(mothId)<1000) myHmothId_1000 -> Fill(mothId);
      if (pdgId==443 && abs(mothId)>=1000 && abs(mothId)<10000) myHmothId_10000 -> Fill(mothId);
      if (pdgId==443 && abs(mothId)>=10000 && abs(mothId)<20000) myHmothId_20000 -> Fill(mothId);
      if (pdgId==443 && abs(mothId)>=20000 && abs(mothId)<100000) myHmothId_100000 -> Fill(mothId);
      if (pdgId==443 && abs(mothId)>=100000 && abs(mothId)<110000) myHmothId_110000 -> Fill(mothId);
      if (pdgId==443 && abs(mothId)>=110000) cout << "XXXXXXXXXXXXXx, mothId>110000" << endl;

      if (pdgId==443) {
	int thisMum = abs(mothId);
	if (!contains(mothNum, thisMum)) mothNum.push_back(thisMum);
      }

      // Mother = B0 (511)
      if (abs(mothId)==511 && abs(pdgId)==443) theBMothIdx = mothIdx;

      // Mother = Xc0(1p) (10441)
      if (abs(mothId)==10441 && abs(pdgId)==443) theChic0MothIdx = mothIdx;

      // Mother = Xc1(1p) (20443)
      if (abs(mothId)==20443 && abs(pdgId)==443) theChic1MothIdx = mothIdx;

      // Mother = Xc2(1p) (445)
      if (abs(mothId)==445 && abs(pdgId)==443) theChic2MothIdx = mothIdx;      // per ora trascuriamo
    }

    
    // Detailed study: when we have B0->JPsi + anything
    bool B0sonJPsi = false;
    bool B0sonK0   = false;
    bool B0sonRho  = false;
    bool B0sonPi   = false;
    bool B0sonKstar0 = false;
    bool B0sonKstarC = false;
    bool B0sonK1     = false;
    bool B0sonKstar2 = false;
    int B0idxJPsi   = -1;
    int B0idxK0     = -1;
    int B0idxRho    = -1;
    int B0idxPi1    = -1;
    int B0idxPi2    = -1;
    int B0idxKstar0 = -1;
    int B0idxKstarC = -1;
    int B0idxK1     = -1;
    int B0idxKstar2 = -1;
    
    if (theBMothIdx>=0) {                                                                                                                                                         
      for (int ii=0; ii<nGenPart; ii++){                                                                                                                                             
	int pdgId = GenPart_pdgId[ii];                                                                                                                                               
	int mothIdx = GenPart_genPartIdxMother[ii];                                                                                                                                  
	if (mothIdx==theBMothIdx) {                                                                                                                                               
	  int mothId = GenPart_pdgId[mothIdx];                                                                                                                                       
	  if (abs(mothId)!=511) continue;
	  
	  // uncomment to check all the possible decays
	  // cou << "Event = " << event << ": MotherId = " << mothId << ", ID = " << pdgId << endl;
	  

	  // which decay
	  if (abs(pdgId)==443) { B0sonJPsi=true; B0idxJPsi = ii; }
	  if (abs(pdgId)==310 || abs(pdgId)==311) { B0sonK0=true; B0idxK0 = ii; }
	  if (abs(pdgId)==113) { B0sonRho=true; B0idxRho = ii; }
	  if (abs(pdgId)==211) { B0sonPi=true;
	    if (B0idxPi1==-1) B0idxPi1 = ii;
	    else B0idxPi2 = ii;
	  }
	  if (abs(pdgId)==313) { B0sonKstar0=true; B0idxKstar0 = ii; }
	  if (abs(pdgId)==323) { B0sonKstarC=true; B0idxKstarC = ii; }
	  if (abs(pdgId)==10313) { B0sonK1=true; B0idxK1 = ii; }
	  if (abs(pdgId)==315) { B0sonKstar2=true; B0idxKstar2 = ii; }
	}
      }

      // K*(892) daughters
      int B0idxKstarCdauK = -1;
      int B0idxKstarCdauP = -1;
      if (B0sonKstarC) {
	for (int ii=0; ii<nGenPart; ii++){   
	  int pdgId = GenPart_pdgId[ii];                 
	  int mothIdx = GenPart_genPartIdxMother[ii]; 
	  if (mothIdx>-1) {                                                                 
	    int mothId = GenPart_pdgId[mothIdx];                            
	    if (abs(mothId)!=323) continue;
	    if (abs(pdgId)==311) B0idxKstarCdauK = ii;
	    if (abs(pdgId)==211) B0idxKstarCdauP = ii;
	  }
	}
      }

      // K1(1270) daughters
      int B0idxK1KstarCdauK = -1;
      int B0idxK1KstarCdauP = -1;
      int B0idxK1dauP = -1;
      if (B0sonK1) {
	for (int ii=0; ii<nGenPart; ii++){   
	  int pdgId = GenPart_pdgId[ii];                 
	  int mothIdx = GenPart_genPartIdxMother[ii]; 
	  if (mothIdx>-1) {                                                                 
	    int mothId = GenPart_pdgId[mothIdx];                            
	    int gmothIdx = GenPart_genPartIdxMother[mothIdx]; 	    
	    if (gmothIdx>-1) {                                                                 	    
	      int gmothId = GenPart_pdgId[gmothIdx];                            
	      if (abs(gmothId)==10313 && abs(mothId)==323 && abs(pdgId)==311) B0idxK1KstarCdauK = ii;
	      if (abs(mothId)==10313 && abs(pdgId)==311) B0idxK1KstarCdauK = ii;
	      if (abs(gmothId)==10313 && abs(mothId)==323 && abs(pdgId)==211) B0idxK1KstarCdauP = ii;
	      if (abs(mothId)==10313 && abs(pdgId)==211 && B0idxK1KstarCdauP==-1) B0idxK1KstarCdauP = ii;
	      if (abs(mothId)==10313 && abs(pdgId)==211 && B0idxK1KstarCdauP>-1) B0idxK1dauP = ii;
	    }
	  }
	}
      }
      
      // uncomment to check all the possible decays
      // cout << endl;

      // Plots decay per decay
      if (B0sonJPsi && B0sonK0 && B0sonPi) {
	// cout << "JPsi PiPi K0: " << B0idxJPsi << " " << B0idxK0 << " " << B0idxPi1 << " " << B0idxPi2 << endl;
	if (B0idxJPsi<0 || B0idxK0<0 || B0idxPi1<0 || B0idxPi2<0) cout << "Error JPsi K0 PiPi" << endl;

	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);

	TLorentzVector P1TLV(0,0,0,0);
	P1TLV.SetPtEtaPhiM(GenPart_pt[B0idxPi1], GenPart_eta[B0idxPi1], GenPart_phi[B0idxPi1], GenPart_mass[B0idxPi1]);

	TLorentzVector P2TLV(0,0,0,0);
	P2TLV.SetPtEtaPhiM(GenPart_pt[B0idxPi2], GenPart_eta[B0idxPi2], GenPart_phi[B0idxPi2], GenPart_mass[B0idxPi2]);

	float genXmass = (JPsiTLV + P1TLV + P2TLV).M();
	  
        myHmassXreco1->Fill(B0_finalFit_X_mass[0]);
        myHmassXgen1->Fill(genXmass);
      }

      if (B0sonJPsi && B0sonK0 && B0sonRho) {
	// cout << "JPsi Rho K0: " << B0idxJPsi << " " << B0idxK0 << " " << B0idxRho << endl;
	if (B0idxJPsi<0 || B0idxK0<0 || B0idxRho<0) cout << "Error JPsi K0 Rho" << endl;

	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);

	TLorentzVector RhoTLV(0,0,0,0);
	RhoTLV.SetPtEtaPhiM(GenPart_pt[B0idxRho], GenPart_eta[B0idxRho], GenPart_phi[B0idxRho], GenPart_mass[B0idxRho]);

	float genXmass = (JPsiTLV + RhoTLV).M();
	  
        myHmassXreco2->Fill(B0_finalFit_X_mass[0]);
        myHmassXgen2->Fill(genXmass);
      }

      if (B0sonJPsi && B0sonKstarC && B0sonPi) {
	// cout << "JPsi K*+(892) pi: " << B0idxJPsi << " " << B0idxKstarC << " " << B0idxPi1 << " " << B0idxKstarCdauK << " " << B0idxKstarCdauP << endl;
	if (B0idxJPsi<0 || B0idxKstarC<0 || B0idxPi1<0 || B0idxKstarCdauK<0 || B0idxKstarCdauP<0) cout << "Error JPsi K*+(892) pi" << endl;

	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);
	
	TLorentzVector P1TLV(0,0,0,0);
	P1TLV.SetPtEtaPhiM(GenPart_pt[B0idxPi1], GenPart_eta[B0idxPi1], GenPart_phi[B0idxPi1], GenPart_mass[B0idxPi1]);
	
	TLorentzVector P2TLV(0,0,0,0);
	P2TLV.SetPtEtaPhiM(GenPart_pt[B0idxPi2], GenPart_eta[B0idxPi2], GenPart_phi[B0idxPi2], GenPart_mass[B0idxPi2]);

	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[B0idxKstarCdauK], GenPart_eta[B0idxKstarCdauK], GenPart_phi[B0idxKstarCdauK], GenPart_mass[B0idxKstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[B0idxKstarCdauP], GenPart_eta[B0idxKstarCdauP], GenPart_phi[B0idxKstarCdauP], GenPart_mass[B0idxKstarCdauP]);

	float genXmassa = (JPsiTLV + DKTLV + DPTLV).M();
	float genXmassb = (JPsiTLV + P1TLV + DKTLV).M();
	float genXmassc = (JPsiTLV + DPTLV + P1TLV).M();
	
        myHmassXreco3->Fill(B0_finalFit_X_mass[0]);
        myHmassXgen3a->Fill(genXmassa);
        myHmassXgen3b->Fill(genXmassb);
        myHmassXgen3c->Fill(genXmassc);
      }
      
      if (B0sonJPsi && B0sonK1) {
	// cout << "ThisIsThe = " << thisIsThe << " => JPsi K1(1270): " << B0idxJPsi << " " << B0idxK1 << " " << B0idxK1KstarCdauK  << " " << B0idxK1KstarCdauP << " " << B0idxK1dauP << endl;
	if (B0idxJPsi<0 || B0idxK1<0 || B0idxK1KstarCdauK<0 || B0idxK1KstarCdauP<0 || B0idxK1dauP<0) cout << "Error JPsi K1" << endl;
	if (B0idxJPsi<0 || B0idxK1<0 || B0idxK1KstarCdauK<0 || B0idxK1KstarCdauP<0 || B0idxK1dauP<0) continue;

      	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[B0idxK1KstarCdauK], GenPart_eta[B0idxK1KstarCdauK], GenPart_phi[B0idxK1KstarCdauK], GenPart_mass[B0idxK1KstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[B0idxK1KstarCdauP], GenPart_eta[B0idxK1KstarCdauP], GenPart_phi[B0idxK1KstarCdauP], GenPart_mass[B0idxK1KstarCdauP]);

	TLorentzVector DK1PTLV(0,0,0,0);
	DK1PTLV.SetPtEtaPhiM(GenPart_pt[B0idxK1dauP], GenPart_eta[B0idxK1dauP], GenPart_phi[B0idxK1dauP], GenPart_mass[B0idxK1dauP]);

	float genXmassa = (JPsiTLV + DKTLV + DPTLV).M();
	float genXmassb = (JPsiTLV + DK1PTLV + DKTLV).M();
	float genXmassc = (JPsiTLV + DK1PTLV + DPTLV).M();
	
        myHmassXreco4->Fill(B0_finalFit_X_mass[0]);
        myHmassXgen4a->Fill(genXmassa);
        myHmassXgen4b->Fill(genXmassb);
        myHmassXgen4c->Fill(genXmassc);
      }
            
    } // Cases with B0 -> JPsi xxx

    
    
    // Detailed study: when we have Xc0->JPsi + anything
    bool Xc0sonJPsi   = false;
    bool Xc0sonKstarC = false;
    int Xc0idxJPsi   = -1;
    int Xc0idxKstarC = -1;
    
    if (theChic0MothIdx>=0) {                                              
      for (int ii=0; ii<nGenPart; ii++){                                           
	int pdgId = GenPart_pdgId[ii];                      
	int mothIdx = GenPart_genPartIdxMother[ii];                        
	if (mothIdx>-1) {                                            
	  int mothId = GenPart_pdgId[mothIdx];
	  
	  // uncomment to check all the possible decays
	  // cout << "Event = " << event << ": MotherId = " << mothId << ", ID = " << pdgId << endl;
	  
	  // which decay
	  if (abs(mothId)==10441 && abs(pdgId)==443) { Xc0sonJPsi=true; Xc0idxJPsi = ii; }
	  if (abs(mothId)==521 && abs(pdgId)==323)   { Xc0sonKstarC=true; Xc0idxKstarC = ii; }
	}
      }
      
      // K*(892) daughters
      int Xc0idxKstarCdauK = -1;
      int Xc0idxKstarCdauP = -1;
      if (Xc0sonKstarC) {
	for (int ii=0; ii<nGenPart; ii++){   
	  int pdgId = GenPart_pdgId[ii];                 
	  int mothIdx = GenPart_genPartIdxMother[ii]; 
	  if (mothIdx>-1) {                                                                 
	    int mothId = GenPart_pdgId[mothIdx];                            
	    if (abs(mothId)!=323) continue;
	    if (abs(pdgId)==311) Xc0idxKstarCdauK = ii;
	    if (abs(pdgId)==211) Xc0idxKstarCdauP = ii;
	  }
	}
      }
      
      // uncomment to check all the possible decays
      // cout << endl;
      if (Xc0sonJPsi && Xc0sonKstarC) {
        // cout << "JPsi K*+(892): " << Xc0idxJPsi << " " << Xc0idxKstarC << " " << Xc0idxKstarCdauK << " " << Xc0idxKstarCdauP << endl;
	if (Xc0idxJPsi<0 || Xc0idxKstarC<0 || Xc0idxKstarCdauK<0 || Xc0idxKstarCdauP<0) cout << "Error JPsi K*+(892)" << endl;

	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[Xc0idxJPsi], GenPart_eta[Xc0idxJPsi], GenPart_phi[Xc0idxJPsi], GenPart_mass[Xc0idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[Xc0idxKstarCdauK], GenPart_eta[Xc0idxKstarCdauK], GenPart_phi[Xc0idxKstarCdauK], GenPart_mass[Xc0idxKstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[Xc0idxKstarCdauP], GenPart_eta[Xc0idxKstarCdauP], GenPart_phi[Xc0idxKstarCdauP], GenPart_mass[Xc0idxKstarCdauP]);

	float genXmass = (JPsiTLV + DKTLV + DPTLV).M();

        myHmassXreco5->Fill(B0_finalFit_X_mass[0]);
        myHmassXgen5->Fill(genXmass);
      }
                 
    } // Cases with Xc0 -> JPsi xxx



    // ----------------------------------------
    // Detailed study: when we have Xc1->JPsi + anything
    bool Xc1sonJPsi   = false;
    bool Xc1sonK      = false;
    bool Xc1sonPi     = false;
    bool Xc1sonKstarC = false;
    int Xc1idxJPsi    = -1;
    int Xc1idxK       = -1;
    int Xc1idxPi      = -1;
    int Xc1idxKstarC  = -1;
    
    if (theChic1MothIdx>=0) {                                              
      for (int ii=0; ii<nGenPart; ii++){                                           
	int pdgId = GenPart_pdgId[ii];                      
	int mothIdx = GenPart_genPartIdxMother[ii];                        
	if (mothIdx>-1) {                                            
	  int mothId = GenPart_pdgId[mothIdx];
	  
	  // uncomment to check all the possible decays
	  // cout << "Event = " << event << ": MotherId = " << mothId << ", ID = " << pdgId << endl;
	  
	  // which decay
	  if (abs(mothId)==20443 && abs(pdgId)==443) { Xc1sonJPsi=true; Xc1idxJPsi = ii; }
	  if ((abs(mothId)==521 || abs(mothId)==511) && abs(pdgId)==321) { Xc1sonK=true; Xc1idxK = ii; }
	  if ((abs(mothId)==521 || abs(mothId)==511) && abs(pdgId)==211) { Xc1sonPi=true; Xc1idxPi = ii; }
	  if ((abs(mothId)==521 || abs(mothId)==511) && abs(pdgId)==323) { Xc1sonKstarC=true; Xc1idxKstarC = ii; }
	}
      }
      
      // K*(892) daughters
      int Xc1idxKstarCdauK = -1;
      int Xc1idxKstarCdauP = -1;
      if (Xc1sonKstarC) {
	for (int ii=0; ii<nGenPart; ii++){   
	  int pdgId = GenPart_pdgId[ii];                 
	  int mothIdx = GenPart_genPartIdxMother[ii]; 
	  if (mothIdx>-1) {                                                                 
	    int mothId = GenPart_pdgId[mothIdx];                            
	    if (abs(mothId)!=323) continue;
	    if (abs(pdgId)==311) Xc1idxKstarCdauK = ii;
	    if (abs(pdgId)==211) Xc1idxKstarCdauP = ii;
	  }
	}
      }
      
      // uncomment to check all the possible decays
      // cout << endl;
      if (Xc1sonJPsi && Xc1sonKstarC) {
	// cout << "Xc1 => JPsi, K*+(892): " << Xc1idxJPsi << " " << Xc1idxKstarC << " " << Xc1idxKstarCdauK << " " << Xc1idxKstarCdauP << endl;
	// if (Xc1idxJPsi<0 || Xc1idxKstarC<0 || Xc1idxKstarCdauK<0 || Xc1idxKstarCdauP<0) cout << "Error JPsi K*+(892)" << endl;

	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxJPsi], GenPart_eta[Xc1idxJPsi], GenPart_phi[Xc1idxJPsi], GenPart_mass[Xc1idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxKstarCdauK], GenPart_eta[Xc1idxKstarCdauK], GenPart_phi[Xc1idxKstarCdauK], GenPart_mass[Xc1idxKstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxKstarCdauP], GenPart_eta[Xc1idxKstarCdauP], GenPart_phi[Xc1idxKstarCdauP], GenPart_mass[Xc1idxKstarCdauP]);

	float genXmass = (JPsiTLV + DKTLV + DPTLV).M();

        myHmassXreco6->Fill(B0_finalFit_X_mass[0]);
        myHmassXgen6->Fill(genXmass);
      }

      if (Xc1sonJPsi && (Xc1sonK || Xc1sonPi)) {
	// cout << "Xc1 => JPsi, K, Pi: " << Xc1idxJPsi << " " << Xc1idxK << " " << Xc1idxPi << endl;

	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxJPsi], GenPart_eta[Xc1idxJPsi], GenPart_phi[Xc1idxJPsi], GenPart_mass[Xc1idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	if (Xc1sonK) DKTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxK], GenPart_eta[Xc1idxK], GenPart_phi[Xc1idxK], GenPart_mass[Xc1idxK]);

	TLorentzVector DPTLV(0,0,0,0);
	if (Xc1sonPi) DPTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxPi], GenPart_eta[Xc1idxPi], GenPart_phi[Xc1idxPi], GenPart_mass[Xc1idxPi]);

	float genXmassa = -1.;
	float genXmassb = -1.;
	if (Xc1sonK && Xc1sonPi) { genXmassa = (JPsiTLV + DKTLV + DPTLV).M(); cout << "con 3: " << genXmassa << endl; }
	else if (Xc1sonK && !Xc1sonPi) { genXmassb = (JPsiTLV + DKTLV).M(); cout << "con 2: " << genXmassb << endl; }
	else if (!Xc1sonK && Xc1sonPi) { genXmassb = (JPsiTLV + DPTLV).M(); cout << "con 2: " << genXmassb << endl; }

        myHmassXreco7->Fill(B0_finalFit_X_mass[0]);
        myHmassXgen7a->Fill(genXmassa);
        myHmassXgen7b->Fill(genXmassb);
      }
                 
    } // Cases with Xc1 -> JPsi xxx
    // ----------------------------------------

  } // entries

    // Which mothers:
  cout << "mothNum.size = " << mothNum.size() << endl;
  for (int ii=0; ii<mothNum.size(); ii++) cout << "Mothers: ii = " << ii << ", " << mothNum[ii] << endl;


  gStyle->SetOptStat(111111);

  TCanvas c0("c0","",1);  
  myHmassXrecoAll->Draw();
  c0.SaveAs("myHmassXrecoAll.png");

  // 
  TCanvas ca1("ca1","",1);
  myHmassXreco2->Draw();
  ca1.SaveAs("myHmassXreco_JPsi_K0_Rho.png");

  TCanvas ca2("ca2","",1);
  myHmassXgen2->Draw();
  ca2.SaveAs("myHmassXgen_JPsi_K0_Rho.png");

  TCanvas cb1("cb1","",1);
  myHmassXreco1->Draw();
  cb1.SaveAs("myHmassXreco_JPsi_K0_PiPi.png");

  TCanvas cb2("cb2","",1);
  myHmassXgen1->Draw();
  cb2.SaveAs("myHmassXgen_JPsi_K0_PiPi.png");

  TCanvas cc1("cc1","",1);
  myHmassXreco3->Draw();
  cc1.SaveAs("myHmassXreco_JPsi_KstarC_Pi.png");

  TCanvas cc2("cc2","",1);
  cc2.Divide(3,1);
  cc2.cd(1); myHmassXgen3a->Draw();
  cc2.cd(2); myHmassXgen3b->Draw("same");
  cc2.cd(3); myHmassXgen3c->Draw("same");
  cc2.SaveAs("myHmassXgen_JPsi_KstarC_Pi.png");

  TCanvas cc3("cc3","",1);
  myHmassXgen3c->Draw("same");
  cc3.SaveAs("myHmassXgen_JPsi_KstarC_Pi___WithPi.png");  

  TCanvas cd1("cd1","",1);
  myHmassXreco4->Draw();
  cd1.SaveAs("myHmassXreco_JPsi_K1.png");

  TCanvas cd2("cd2","",1);
  cd2.Divide(3,1);
  cd2.cd(1); myHmassXgen4a->Draw();
  cd2.cd(2); myHmassXgen4b->Draw("same");
  cd2.cd(3); myHmassXgen4c->Draw("same");
  cd2.SaveAs("myHmassXgen_JPsi_K1.png");

  TCanvas cd3("cd3","",1);
  myHmassXgen4c->Draw("same");
  cd3.SaveAs("myHmassXgen_JPsi_K1___WithPi.png");  

  TCanvas ce1("ce1","",1);
  myHmassXreco5->Draw();
  ce1.SaveAs("myHmassXreco_Xc0_JPsi_KstarC.png");

  TCanvas ce2("ce2","",1);
  myHmassXgen5->Draw();
  ce2.SaveAs("myHmassXgen_Xc0_JPsi_KstarC.png");

  TCanvas cf1("cf1","",1);
  myHmassXreco6->Draw();
  cf1.SaveAs("myHmassXreco_Xc1_JPsi_KstarC.png");

  TCanvas cf2("cf2","",1);
  myHmassXgen6->Draw();
  cf2.SaveAs("myHmassXgen_Xc1_JPsi_KstarC.png");
  
  TCanvas cg1("cg1","",1);
  myHmassXreco7->Draw();
  cg1.SaveAs("myHmassXreco_Xc1_JPsi_KoP.png");

  TCanvas cg2("cg2","",1);
  cg2.Divide(2,1);
  cg2.cd(1); myHmassXgen7a->Draw();
  cg2.cd(2); myHmassXgen7b->Draw();
  cg2.SaveAs("myHmassXgen_Xc1_JPsi_KoP.png");

  //TCanvas cd3("cd3","",1);
  //myHmassXgen4c->Draw("same");
  //cd3.SaveAs("myHmassXgen_JPsi_K1___WithPi.png");  

  
  
  TCanvas c2("c2","",1);
  myHmothId_500->Draw();
  c2.SaveAs("mothId_500.png");

  TCanvas c3("c3","",1);
  myHmothId_520->Draw();
  c3.SaveAs("mothId_520.png");

  TCanvas c4("c4","",1);
  myHmothId_1000->Draw();
  c4.SaveAs("mothId_1000.png");

  TCanvas c5("c5","",1);
  myHmothId_10000->Draw();
  c5.SaveAs("mothId_10000.png");

  TCanvas c6("c6","",1);
  myHmothId_20000->Draw();
  c6.SaveAs("mothId_20000.png");

  TCanvas c7("c7","",1);
  myHmothId_100000->Draw();
  c7.SaveAs("mothId_100000.png");

  TCanvas c8("c8","",1);
  myHmothId_110000->Draw();
  c8.SaveAs("mothId_110000.png");
}

