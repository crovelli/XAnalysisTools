#include "../include/GenLevelStudy.hh"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TCanvas.h>
#include <TH2.h>
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
  TH1F *myHmassXrecoAll = new TH1F("myHmassXrecoAll", "myHmassXrecoAll", 70, 3.4, 4.1);

  // Mass plots separate per decay
  TH1F *myHmassXgen1      = new TH1F("myHmassXgen1",      "myHmassXgen1",     100, 2.5, 5.5);
  TH1F *myHmassXgen1z     = new TH1F("myHmassXgen1z",     "myHmassXgen1z",     70, 3.4, 4.1);
  TH1F *myHmassXreco1     = new TH1F("myHmassXreco1",     "myHmassXreco1",     70, 3.4, 4.1);
  TH1F *myHmassXreco1m    = new TH1F("myHmassXreco1m",    "myHmassXreco1m",    70, 3.4, 4.1);  
  TH1F *myHmassXreco1pre  = new TH1F("myHmassXreco1pre",  "myHmassXreco1pre",  70, 3.4, 4.1);
  TH1F *myHmassXreco1mpre = new TH1F("myHmassXreco1mpre", "myHmassXreco1mpre", 70, 3.4, 4.1);
  // 
  TH1F *myHmassXgen2     = new TH1F("myHmassXgen2",     "myHmassXgen2",    100, 2.5, 5.5);
  TH1F *myHmassXgen2z    = new TH1F("myHmassXgen2z",    "myHmassXgen2",     70, 3.4, 4.1);
  TH1F *myHmassXreco2    = new TH1F("myHmassXreco2",    "myHmassXreco2",    70, 3.4, 4.1);
  TH1F *myHmassXreco2pre = new TH1F("myHmassXreco2pre", "myHmassXreco2pre", 70, 3.4, 4.1);
  //
  TH1F *myHmassXgen3a    = new TH1F("myHmassXgen3a",    "myHmassXgen3a",   100, 2.5, 5.5);
  TH1F *myHmassXgen3b    = new TH1F("myHmassXgen3b",    "myHmassXgen3b",   100, 2.5, 5.5);
  TH1F *myHmassXgen3c    = new TH1F("myHmassXgen3c",    "myHmassXgen3c",   100, 2.5, 5.5);
  TH1F *myHmassXgen3az   = new TH1F("myHmassXgen3az",   "myHmassXgen3az",   70, 3.4, 4.1);
  TH1F *myHmassXgen3bz   = new TH1F("myHmassXgen3bz",   "myHmassXgen3bz",   70, 3.4, 4.1);
  TH1F *myHmassXgen3cz   = new TH1F("myHmassXgen3cz",   "myHmassXgen3cz",   70, 3.4, 4.1);
  TH1F *myHmassXreco3    = new TH1F("myHmassXreco3",    "myHmassXreco3",    70, 3.4, 4.1);
  TH1F *myHmassXreco3pre = new TH1F("myHmassXreco3pre", "myHmassXreco3pre", 70, 3.4, 4.1);
  // 
  TH1F *myHmassXgen4a    = new TH1F("myHmassXgen4a",    "myHmassXgen4a",   100, 2.5, 5.5);
  TH1F *myHmassXgen4b    = new TH1F("myHmassXgen4b",    "myHmassXgen4b",   100, 2.5, 5.5);
  TH1F *myHmassXgen4c    = new TH1F("myHmassXgen4c",    "myHmassXgen4c",   100, 2.5, 5.5);
  TH1F *myHmassXgen4az   = new TH1F("myHmassXgen4az",   "myHmassXgen4a",    70, 3.4, 4.1);
  TH1F *myHmassXgen4bz   = new TH1F("myHmassXgen4bz",   "myHmassXgen4b",    70, 3.4, 4.1);
  TH1F *myHmassXgen4cz   = new TH1F("myHmassXgen4cz",   "myHmassXgen4c",    70, 3.4, 4.1);
  TH1F *myHmassXreco4    = new TH1F("myHmassXreco4",    "myHmassXreco4",    70, 3.4, 4.1);
  TH1F *myHmassXreco4pre = new TH1F("myHmassXreco4pre", "myHmassXreco4pre", 70, 3.4, 4.1);
  // 
  TH1F *myHmassXgen5     = new TH1F("myHmassXgen5",     "myHmassXgen5",    100, 2.5, 5.5);
  TH1F *myHmassXgen5z    = new TH1F("myHmassXgen5z",    "myHmassXgen5z",    70, 3.4, 4.1);
  TH1F *myHmassXreco5    = new TH1F("myHmassXreco5",    "myHmassXreco5",    70, 3.4, 4.1);
  TH1F *myHmassXreco5pre = new TH1F("myHmassXreco5pre", "myHmassXreco5pre", 70, 3.4, 4.1);
  //
  TH1F *myHmassXgen6     = new TH1F("myHmassXgen6",     "myHmassXgen6",    100, 2.5, 5.5);
  TH1F *myHmassXgen6z    = new TH1F("myHmassXgen6z",    "myHmassXgen6z",    70, 3.4, 4.1);
  TH1F *myHmassXreco6    = new TH1F("myHmassXreco6",    "myHmassXreco6",    70, 3.4, 4.1);
  TH1F *myHmassXreco6pre = new TH1F("myHmassXreco6pre", "myHmassXreco6pre", 70, 3.4, 4.1);
  //
  TH1F *myHmassXgen7a     = new TH1F("myHmassXgen7a",     "myHmassXgen7a",    100, 2.5, 5.5);
  TH1F *myHmassXgen7b     = new TH1F("myHmassXgen7b",     "myHmassXgen7b",    100, 2.5, 5.5);
  TH1F *myHmassXgen7az    = new TH1F("myHmassXgen7az",    "myHmassXgen7az",    70, 3.4, 4.1);
  TH1F *myHmassXgen7bz    = new TH1F("myHmassXgen7bz",    "myHmassXgen7bz",    70, 3.4, 4.1);
  TH1F *myHmassXreco7     = new TH1F("myHmassXreco7",     "myHmassXreco7",     70, 3.4, 4.1);
  TH1F *myHmassXreco7pre  = new TH1F("myHmassXreco7pre" , "myHmassXreco7pre",  70, 3.4, 4.1);

  TH2F *myH2minDr1  = new TH2F("myH2minDr1", "myH2minDr1", 50,0.,0.2,50, 0.,0.2);
  TH2F *myH2minDr1z = new TH2F("myH2minDr1z","myH2minDr1z",50,0.,0.02,50, 0.,0.02);

  // For summary plots
  TH1F *myHmassXrecoB0        = new TH1F("myHmassXrecoB0" ,       "myHmassXrecoB0",        35, 3.4, 4.1);
  TH1F *myHmassXrecoB0Pre     = new TH1F("myHmassXrecoB0pre",     "myHmassXrecoB0pre",     35, 3.4, 4.1);
  TH1F *myHmassXrecoBPlus     = new TH1F("myHmassXrecoBPlus" ,    "myHmassXrecoBPlus",     35, 3.4, 4.1);
  TH1F *myHmassXrecoBPlusPre  = new TH1F("myHmassXrecoBPluspre" , "myHmassXrecoBPluspre",  35, 3.4, 4.1);
  TH1F *myHmassXrecoChic0     = new TH1F("myHmassXrecoChic0" ,    "myHmassXrecoChic0",     35, 3.4, 4.1);
  TH1F *myHmassXrecoChic0Pre  = new TH1F("myHmassXrecoChic0pre" , "myHmassXrecoChic0pre",  35, 3.4, 4.1);
  TH1F *myHmassXrecoChic1     = new TH1F("myHmassXrecoChic1" ,    "myHmassXrecoChic1",     35, 3.4, 4.1);
  TH1F *myHmassXrecoChic1Pre  = new TH1F("myHmassXrecoChic1pre" , "myHmassXrecoChic1pre",  35, 3.4, 4.1);
  TH1F *myHmassXrecoChic2     = new TH1F("myHmassXrecoChic2" ,    "myHmassXrecoChic2",     35, 3.4, 4.1);
  TH1F *myHmassXrecoChic2Pre  = new TH1F("myHmassXrecoChic2pre" , "myHmassXrecoChic2pre",  35, 3.4, 4.1);
  
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
    bool atLeastOneB = false;
    for (int nB=0; nB<nB0; nB++) {
      if (B0_K0s_nmcFitted_mass[nB]<0.541 && B0_K0s_nmcFitted_mass[nB]>0.454) {
	if (B0_finalFit_mass[nB]<5.35 && B0_finalFit_mass[nB]>5.20) {
	  if (B0_finalFit_X_mass[nB]<4.1 && B0_finalFit_X_mass[nB]>3.4) { // X: 3.81-3.93; Psi2s: 3.65-3.72
	    myHmassXrecoAll->Fill(B0_finalFit_X_mass[nB]);
	    atLeastOneB = true;
	  }}}
    }
    if (!atLeastOneB) continue;

    // Reco level pre-fit
    TLorentzVector PR1TLV(0,0,0,0);
    TLorentzVector PR2TLV(0,0,0,0);
    TLorentzVector MR1TLV(0,0,0,0);
    TLorentzVector MR2TLV(0,0,0,0);
    PR1TLV.SetPtEtaPhiM(B0_PiPi_prefit_pi1_pt[0], B0_PiPi_prefit_pi1_eta[0], B0_PiPi_prefit_pi1_phi[0], 0.139571);
    PR2TLV.SetPtEtaPhiM(B0_PiPi_prefit_pi2_pt[0], B0_PiPi_prefit_pi2_eta[0], B0_PiPi_prefit_pi2_phi[0], 0.139571);
    MR1TLV.SetPtEtaPhiM(B0_MuMu_prefit_mu1_pt[0], B0_MuMu_prefit_mu1_eta[0], B0_MuMu_prefit_mu1_phi[0], 0.10565837);
    MR2TLV.SetPtEtaPhiM(B0_MuMu_prefit_mu2_pt[0], B0_MuMu_prefit_mu2_eta[0], B0_MuMu_prefit_mu2_phi[0], 0.10565837);
    TVector3 PR1TV3(0,0,0);
    TVector3 PR2TV3(0,0,0);
    PR1TV3.SetPtEtaPhi(B0_PiPi_prefit_pi1_pt[0], B0_PiPi_prefit_pi1_eta[0], B0_PiPi_prefit_pi1_phi[0]);
    PR2TV3.SetPtEtaPhi(B0_PiPi_prefit_pi2_pt[0], B0_PiPi_prefit_pi2_eta[0], B0_PiPi_prefit_pi2_phi[0]);
    //
    float recoXmassPreFit = (PR1TLV + PR2TLV + MR1TLV + MR2TLV).M();
    
    // Gen level: which is the jpsi mother
    int theBMothIdx = -1;
    int theBPlusMothIdx = -1;
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

      // Mother = B+ (521)
      if (abs(mothId)==521 && abs(pdgId)==443) theBPlusMothIdx = mothIdx;

      // Mother = Xc0(1p) (10441)
      if (abs(mothId)==10441 && abs(pdgId)==443) theChic0MothIdx = mothIdx;

      // Mother = Xc1(1p) (20443)
      if (abs(mothId)==20443 && abs(pdgId)==443) theChic1MothIdx = mothIdx;

      // Mother = Xc2(1p) (445)
      if (abs(mothId)==445 && abs(pdgId)==443) theChic2MothIdx = mothIdx;      // per ora trascuriamo
    }


    // ----------------------------------------------------
    // Detailed study: when we have B0->JPsi + anything
    bool B0sonJPsi   = false;
    bool B0sonK0     = false;
    bool B0sonRho    = false;
    bool B0sonPi     = false;
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

      // B-> JPsi PiPi K0
      if (B0sonJPsi && B0sonK0 && B0sonPi) {
	if (B0idxJPsi<0 || B0idxK0<0 || B0idxPi1<0 || B0idxPi2<0) cout << "Error JPsi K0 PiPi" << endl;
	
	// gen level
	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);

	TLorentzVector P1TLV(0,0,0,0);
	P1TLV.SetPtEtaPhiM(GenPart_pt[B0idxPi1], GenPart_eta[B0idxPi1], GenPart_phi[B0idxPi1], GenPart_mass[B0idxPi1]);
	TVector3 P1TV3(0,0,0);
	P1TV3.SetPtEtaPhi(GenPart_pt[B0idxPi1], GenPart_eta[B0idxPi1], GenPart_phi[B0idxPi1]);

	TLorentzVector P2TLV(0,0,0,0);
	P2TLV.SetPtEtaPhiM(GenPart_pt[B0idxPi2], GenPart_eta[B0idxPi2], GenPart_phi[B0idxPi2], GenPart_mass[B0idxPi2]);
	TVector3 P2TV3(0,0,0);
	P2TV3.SetPtEtaPhi(GenPart_pt[B0idxPi2], GenPart_eta[B0idxPi2], GenPart_phi[B0idxPi2]);

	// gen-reco match for pions
	float deltaR1 = 999;
	float deltaR2 = 999;
	float dR11 = PR1TV3.DeltaR(P1TV3);
	float dR21 = PR2TV3.DeltaR(P1TV3);
	float dR12 = PR1TV3.DeltaR(P2TV3);
	float dR22 = PR2TV3.DeltaR(P2TV3);
	if (dR11<dR21) deltaR1 = dR11;
	else deltaR1 = dR21;
	if (dR12<dR22) deltaR2 = dR12;
	else deltaR2 = dR22;

	// histos
	float genXmass = (JPsiTLV + P1TLV + P2TLV).M();
        myHmassXreco1->Fill(B0_finalFit_X_mass[0]);
        myHmassXreco1pre->Fill(recoXmassPreFit);

        if (deltaR1<0.01 && deltaR2<0.01) {
	  myHmassXreco1m->Fill(B0_finalFit_X_mass[0]);
	  myHmassXreco1mpre->Fill(recoXmassPreFit);
	}
        myHmassXgen1->Fill(genXmass);
        myHmassXgen1z->Fill(genXmass);
	myH2minDr1->Fill(deltaR1,deltaR2);
	myH2minDr1z->Fill(deltaR1,deltaR2);
      }


      // B-> JPsi Rho K0
      if (B0sonJPsi && B0sonK0 && B0sonRho) {
	if (B0idxJPsi<0 || B0idxK0<0 || B0idxRho<0) cout << "Error JPsi K0 Rho" << endl;

	// gen level
	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);

	TLorentzVector RhoTLV(0,0,0,0);
	RhoTLV.SetPtEtaPhiM(GenPart_pt[B0idxRho], GenPart_eta[B0idxRho], GenPart_phi[B0idxRho], GenPart_mass[B0idxRho]);

	// histos
	float genXmass = (JPsiTLV + RhoTLV).M();
	myHmassXreco2->Fill(B0_finalFit_X_mass[0]);
        myHmassXreco2pre->Fill(recoXmassPreFit);
        myHmassXgen2->Fill(genXmass);
        myHmassXgen2z->Fill(genXmass);
      }


      // JPsi K*+(892) -> Kpi pi (+ gamma)
      if (B0sonJPsi && B0sonKstarC && B0sonPi) {
	if (B0idxJPsi<0 || B0idxKstarC<0 || B0idxPi1<0 || B0idxKstarCdauK<0 || B0idxKstarCdauP<0) cout << "Error JPsi K*+(892) pi" << endl;

	// gen level
	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);
	
	TLorentzVector P1TLV(0,0,0,0);
	P1TLV.SetPtEtaPhiM(GenPart_pt[B0idxPi1], GenPart_eta[B0idxPi1], GenPart_phi[B0idxPi1], GenPart_mass[B0idxPi1]);
	
	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[B0idxKstarCdauK], GenPart_eta[B0idxKstarCdauK], GenPart_phi[B0idxKstarCdauK], GenPart_mass[B0idxKstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[B0idxKstarCdauP], GenPart_eta[B0idxKstarCdauP], GenPart_phi[B0idxKstarCdauP], GenPart_mass[B0idxKstarCdauP]);

	// histos
	float genXmassa = (JPsiTLV + DKTLV + DPTLV).M();
	float genXmassb = (JPsiTLV + P1TLV + DKTLV).M();
	float genXmassc = (JPsiTLV + DPTLV + P1TLV).M();

        myHmassXreco3->Fill(B0_finalFit_X_mass[0]);
        myHmassXreco3pre->Fill(recoXmassPreFit);
        myHmassXgen3a->Fill(genXmassa);
        myHmassXgen3b->Fill(genXmassb);
        myHmassXgen3c->Fill(genXmassc);
        myHmassXgen3az->Fill(genXmassa);
        myHmassXgen3bz->Fill(genXmassb);
        myHmassXgen3cz->Fill(genXmassc);
      }
      

      if (B0sonJPsi && B0sonK1) {
	if (B0idxJPsi<0 || B0idxK1<0 || B0idxK1KstarCdauK<0 || B0idxK1KstarCdauP<0 || B0idxK1dauP<0) cout << "Error JPsi K1" << endl;
	if (B0idxJPsi<0 || B0idxK1<0 || B0idxK1KstarCdauK<0 || B0idxK1KstarCdauP<0 || B0idxK1dauP<0) continue;

	// gen level
      	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[B0idxJPsi], GenPart_eta[B0idxJPsi], GenPart_phi[B0idxJPsi], GenPart_mass[B0idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[B0idxK1KstarCdauK], GenPart_eta[B0idxK1KstarCdauK], GenPart_phi[B0idxK1KstarCdauK], GenPart_mass[B0idxK1KstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[B0idxK1KstarCdauP], GenPart_eta[B0idxK1KstarCdauP], GenPart_phi[B0idxK1KstarCdauP], GenPart_mass[B0idxK1KstarCdauP]);

	TLorentzVector DK1PTLV(0,0,0,0);
	DK1PTLV.SetPtEtaPhiM(GenPart_pt[B0idxK1dauP], GenPart_eta[B0idxK1dauP], GenPart_phi[B0idxK1dauP], GenPart_mass[B0idxK1dauP]);

	// histos
	float genXmassa = (JPsiTLV + DKTLV + DPTLV).M();
	float genXmassb = (JPsiTLV + DK1PTLV + DKTLV).M();
	float genXmassc = (JPsiTLV + DK1PTLV + DPTLV).M();
	
        myHmassXreco4->Fill(B0_finalFit_X_mass[0]);
        myHmassXreco4pre->Fill(recoXmassPreFit);
        myHmassXgen4a->Fill(genXmassa);
        myHmassXgen4b->Fill(genXmassb);
        myHmassXgen4c->Fill(genXmassc);
        myHmassXgen4az->Fill(genXmassa);
        myHmassXgen4bz->Fill(genXmassb);
        myHmassXgen4cz->Fill(genXmassc);
      }
            
    } // Cases with B0 -> JPsi xxx

    
    // -------------------------------------------------
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
	if (Xc0idxJPsi<0 || Xc0idxKstarC<0 || Xc0idxKstarCdauK<0 || Xc0idxKstarCdauP<0) cout << "Error JPsi K*+(892)" << endl;

	// Gen level
	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[Xc0idxJPsi], GenPart_eta[Xc0idxJPsi], GenPart_phi[Xc0idxJPsi], GenPart_mass[Xc0idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[Xc0idxKstarCdauK], GenPart_eta[Xc0idxKstarCdauK], GenPart_phi[Xc0idxKstarCdauK], GenPart_mass[Xc0idxKstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[Xc0idxKstarCdauP], GenPart_eta[Xc0idxKstarCdauP], GenPart_phi[Xc0idxKstarCdauP], GenPart_mass[Xc0idxKstarCdauP]);

	// histos
	float genXmass = (JPsiTLV + DKTLV + DPTLV).M();
        myHmassXreco5->Fill(B0_finalFit_X_mass[0]);
        myHmassXreco5pre->Fill(recoXmassPreFit);
        myHmassXgen5->Fill(genXmass);
        myHmassXgen5z->Fill(genXmass);
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

	// Gen Level
	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxJPsi], GenPart_eta[Xc1idxJPsi], GenPart_phi[Xc1idxJPsi], GenPart_mass[Xc1idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	DKTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxKstarCdauK], GenPart_eta[Xc1idxKstarCdauK], GenPart_phi[Xc1idxKstarCdauK], GenPart_mass[Xc1idxKstarCdauK]);

	TLorentzVector DPTLV(0,0,0,0);
	DPTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxKstarCdauP], GenPart_eta[Xc1idxKstarCdauP], GenPart_phi[Xc1idxKstarCdauP], GenPart_mass[Xc1idxKstarCdauP]);

	// histos
	float genXmass = (JPsiTLV + DKTLV + DPTLV).M();
        myHmassXreco6->Fill(B0_finalFit_X_mass[0]);
        myHmassXreco6pre->Fill(recoXmassPreFit);
        myHmassXgen6->Fill(genXmass);
        myHmassXgen6z->Fill(genXmass);
      }

      if (Xc1sonJPsi && (Xc1sonK || Xc1sonPi)) {
	// cout << "Xc1 => JPsi, K, Pi: " << Xc1idxJPsi << " " << Xc1idxK << " " << Xc1idxPi << endl;

	// Gen Level
	TLorentzVector JPsiTLV(0,0,0,0);
	JPsiTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxJPsi], GenPart_eta[Xc1idxJPsi], GenPart_phi[Xc1idxJPsi], GenPart_mass[Xc1idxJPsi]);
	
	TLorentzVector DKTLV(0,0,0,0);
	if (Xc1sonK) DKTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxK], GenPart_eta[Xc1idxK], GenPart_phi[Xc1idxK], GenPart_mass[Xc1idxK]);

	TLorentzVector DPTLV(0,0,0,0);
	if (Xc1sonPi) DPTLV.SetPtEtaPhiM(GenPart_pt[Xc1idxPi], GenPart_eta[Xc1idxPi], GenPart_phi[Xc1idxPi], GenPart_mass[Xc1idxPi]);

	// histos
	float genXmassa = -1.;
	float genXmassb = -1.;
	if (Xc1sonK && Xc1sonPi) { genXmassa = (JPsiTLV + DKTLV + DPTLV).M(); }
	else if (Xc1sonK && !Xc1sonPi) { genXmassb = (JPsiTLV + DKTLV).M(); }
	else if (!Xc1sonK && Xc1sonPi) { genXmassb = (JPsiTLV + DPTLV).M(); }

        myHmassXreco7->Fill(B0_finalFit_X_mass[0]);
        myHmassXreco7pre->Fill(recoXmassPreFit);
        myHmassXgen7a->Fill(genXmassa);
        myHmassXgen7b->Fill(genXmassb);
        myHmassXgen7az->Fill(genXmassa);
        myHmassXgen7bz->Fill(genXmassb);
      }
                 
    } // Cases with Xc1 -> JPsi xxx


    // ----------------------------------------
    // Summary plots
    if (theBMothIdx>=0) {                                                                                                                                                         
      myHmassXrecoB0->Fill(B0_finalFit_X_mass[0]);
      myHmassXrecoB0Pre->Fill(recoXmassPreFit);
    } else if (theBPlusMothIdx>=0) {                                                                                                                                                         
      myHmassXrecoBPlus->Fill(B0_finalFit_X_mass[0]);
      myHmassXrecoBPlusPre->Fill(recoXmassPreFit);
    } else if (theChic0MothIdx>=0) {                                                                                                                                                         
      myHmassXrecoChic0->Fill(B0_finalFit_X_mass[0]);
      myHmassXrecoChic0Pre->Fill(recoXmassPreFit);
    } else if (theChic1MothIdx>=0) {                                                                                                                                                         
      myHmassXrecoChic1->Fill(B0_finalFit_X_mass[0]);
      myHmassXrecoChic1Pre->Fill(recoXmassPreFit);
    } else if (theChic2MothIdx>=0) {                                                                                                                                                         
      myHmassXrecoChic2->Fill(B0_finalFit_X_mass[0]);
      myHmassXrecoChic2Pre->Fill(recoXmassPreFit);
    }


  } // entries

    // Which mothers:
  cout << "mothNum.size = " << mothNum.size() << endl;
  for (int ii=0; ii<mothNum.size(); ii++) cout << "Mothers: ii = " << ii << ", " << mothNum[ii] << endl;


  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);

  TCanvas c0("c0","",1);  
  myHmassXrecoAll->SetTitle("");
  myHmassXrecoAll->GetXaxis()->SetTitle("m('X')");
  myHmassXrecoAll->Draw();
  c0.SaveAs("myHmassXrecoAll.png");

  // --------------------------------
  TCanvas ca1("ca1","",1);
  myHmassXreco2->Draw();
  ca1.SaveAs("myHmassXreco_JPsi_K0_Rho.png");

  TCanvas ca1p("ca1p","",1);
  myHmassXreco2pre->Draw();
  ca1p.SaveAs("myHmassXrecoPrefit_JPsi_K0_Rho.png");

  TCanvas ca2("ca2","",1);
  myHmassXgen2->Draw();
  ca2.SaveAs("myHmassXgen_JPsi_K0_Rho.png");

  TCanvas ca2z("ca2z","",1);
  myHmassXgen2z->Draw();
  ca2z.SaveAs("myHmassXgenZoom_JPsi_K0_Rho.png");

  // --------------------------------
  TCanvas cb1("cb1","",1);
  myHmassXreco1->Draw();
  cb1.SaveAs("myHmassXreco_JPsi_K0_PiPi.png");

  TCanvas cb1p("cb1p","",1);
  myHmassXreco1pre->Draw();
  cb1p.SaveAs("myHmassXrecoPrefit_JPsi_K0_PiPi.png");

  TCanvas cb1mp("cb1mp","",1);
  myHmassXreco1mpre->Draw();
  cb1mp.SaveAs("myHmassXrecoMatchedPrefit_JPsi_K0_PiPi.png");

  TCanvas cb1m("cb1m","",1);
  myHmassXreco1m->Draw();
  cb1m.SaveAs("myHmassXrecoMatched_JPsi_K0_PiPi.png");

  TCanvas cb2("cb2","",1);
  myHmassXgen1->Draw();
  cb2.SaveAs("myHmassXgen_JPsi_K0_PiPi.png");

  TCanvas cb2z("cb2z","",1);
  myHmassXgen1z->Draw();
  cb2z.SaveAs("myHmassXgenZoom_JPsi_K0_PiPi.png");

  TCanvas cb3("cb3","",1);
  myH2minDr1->Draw("colz");
  cb3.SaveAs("myMinDr_JPsi_K0_PiPi.png");

  TCanvas cb3z("cb3z","",1);
  myH2minDr1z->Draw("colz");
  cb3z.SaveAs("myMinDrZoom_JPsi_K0_PiPi.png");


  // --------------------------------
  TCanvas cc1("cc1","",1);
  myHmassXreco3->Draw();
  cc1.SaveAs("myHmassXreco_JPsi_KstarC_Pi.png");

  TCanvas cc1p("cc1p","",1);
  myHmassXreco3pre->Draw();
  cc1p.SaveAs("myHmassXrecoPrefit_JPsi_KstarC_Pi.png");

  TCanvas cc2("cc2","",1);
  cc2.Divide(3,1);
  cc2.cd(1); myHmassXgen3a->Draw();
  cc2.cd(2); myHmassXgen3b->Draw("same");
  cc2.cd(3); myHmassXgen3c->Draw("same");
  cc2.SaveAs("myHmassXgen_JPsi_KstarC_Pi.png");

  TCanvas cc2z("cc2z","",1);
  cc2z.Divide(3,1);
  cc2z.cd(1); myHmassXgen3az->Draw();
  cc2z.cd(2); myHmassXgen3bz->Draw("same");
  cc2z.cd(3); myHmassXgen3cz->Draw("same");
  cc2z.SaveAs("myHmassXgenZoom_JPsi_KstarC_Pi.png");

  TCanvas cc3("cc3","",1);
  myHmassXgen3c->Draw("same");
  cc3.SaveAs("myHmassXgen_JPsi_KstarC_Pi___WithPi.png");  

  TCanvas cc3z("cc3z","",1);
  myHmassXgen3cz->Draw("same");
  cc3z.SaveAs("myHmassXgenZoom_JPsi_KstarC_Pi___WithPi.png");  

  TCanvas cd1("cd1","",1);
  myHmassXreco4->Draw();
  cd1.SaveAs("myHmassXreco_JPsi_K1.png");

  TCanvas cd1p("cd1p","",1);
  myHmassXreco4pre->Draw();
  cd1p.SaveAs("myHmassXrecoPrefit_JPsi_K1.png");

  TCanvas cd2("cd2","",1);
  cd2.Divide(3,1);
  cd2.cd(1); myHmassXgen4a->Draw();
  cd2.cd(2); myHmassXgen4b->Draw("same");
  cd2.cd(3); myHmassXgen4c->Draw("same");
  cd2.SaveAs("myHmassXgen_JPsi_K1.png");

  TCanvas cd2z("cd2z","",1);
  cd2z.Divide(3,1);
  cd2z.cd(1); myHmassXgen4az->Draw();
  cd2z.cd(2); myHmassXgen4bz->Draw("same");
  cd2z.cd(3); myHmassXgen4cz->Draw("same");
  cd2z.SaveAs("myHmassXgenZoom_JPsi_K1.png");

  TCanvas cd3("cd3","",1);
  myHmassXgen4cz->Draw("same");
  cd3.SaveAs("myHmassXgenZoom_JPsi_K1___WithPi.png");  

  TCanvas ce1("ce1","",1);
  myHmassXreco5->Draw();
  ce1.SaveAs("myHmassXreco_Xc0_JPsi_KstarC.png");

  TCanvas ce1p("ce1p","",1);
  myHmassXreco5pre->Draw();
  ce1p.SaveAs("myHmassXrecoPrefit_Xc0_JPsi_KstarC.png");

  TCanvas ce2("ce2","",1);
  myHmassXgen5->Draw();
  ce2.SaveAs("myHmassXgen_Xc0_JPsi_KstarC.png");

  TCanvas ce2z("ce2z","",1);
  myHmassXgen5z->Draw();
  ce2z.SaveAs("myHmassXgenZoom_Xc0_JPsi_KstarC.png");

  TCanvas cf1("cf1","",1);
  myHmassXreco6->Draw();
  cf1.SaveAs("myHmassXreco_Xc1_JPsi_KstarC.png");

  TCanvas cf1p("cf1p","",1);
  myHmassXreco6pre->Draw();
  cf1p.SaveAs("myHmassXrecoPrefit_Xc1_JPsi_KstarC.png");

  TCanvas cf2("cf2","",1);
  myHmassXgen6->Draw();
  cf2.SaveAs("myHmassXgen_Xc1_JPsi_KstarC.png");

  TCanvas cf2z("cf2z","",1);
  myHmassXgen6z->Draw();
  cf2z.SaveAs("myHmassXgenZoom_Xc1_JPsi_KstarC.png");

  TCanvas cg1("cg1","",1);
  myHmassXreco7->Draw();
  cg1.SaveAs("myHmassXreco_Xc1_JPsi_KoP.png");

  TCanvas cg1p("cg1p","",1);
  myHmassXreco7pre->Draw();
  cg1p.SaveAs("myHmassXrecoPrefit_Xc1_JPsi_KoP.png");

  TCanvas cg2("cg2","",1);
  cg2.Divide(2,1);
  cg2.cd(1); myHmassXgen7a->Draw();
  cg2.cd(2); myHmassXgen7b->Draw();
  cg2.SaveAs("myHmassXgen_Xc1_JPsi_KoP.png");

  TCanvas cg2z("cg2z","",1);
  cg2z.Divide(2,1);
  cg2z.cd(1); myHmassXgen7az->Draw();
  cg2z.cd(2); myHmassXgen7bz->Draw();
  cg2z.SaveAs("myHmassXgenZoom_Xc1_JPsi_KoP.png");


  // ------------------------------------------------
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


  // summary plots
  myHmassXrecoB0        -> SetTitle("");
  myHmassXrecoB0        -> SetLineColor(1);
  myHmassXrecoB0Pre     -> SetTitle("");
  myHmassXrecoB0Pre     -> SetLineColor(1);
  myHmassXrecoBPlus     -> SetLineColor(2);
  myHmassXrecoBPlusPre  -> SetLineColor(2);
  myHmassXrecoChic0     -> SetLineColor(3);  
  myHmassXrecoChic0Pre  -> SetLineColor(3);
  myHmassXrecoChic1     -> SetLineColor(4);
  myHmassXrecoChic1Pre  -> SetLineColor(4);
  myHmassXrecoChic2     -> SetLineColor(6);
  myHmassXrecoChic2Pre  -> SetLineColor(6);

  // Comparison between decays from different hadrons
  TCanvas c9a("c9a","",1);
  myHmassXrecoB0->SetTitle("");
  myHmassXrecoB0->GetXaxis()->SetTitle("m('X')");
  myHmassXrecoB0        -> DrawNormalized();
  myHmassXrecoBPlus     -> DrawNormalized("same");
  //myHmassXrecoChic0     -> DrawNormalized("same");  
  myHmassXrecoChic1     -> DrawNormalized("same");
  myHmassXrecoChic2     -> DrawNormalized("same");
  c9a.SaveAs("SummaryFromBHadMassX.png");

  TCanvas c9b("c9b","",1);
  myHmassXrecoB0Pre->SetTitle("");
  myHmassXrecoB0Pre->GetXaxis()->SetTitle("m('X')");
  myHmassXrecoB0Pre     -> DrawNormalized();    
  myHmassXrecoBPlusPre  -> DrawNormalized("same");
  //myHmassXrecoChic0Pre  -> DrawNormalized("same");
  myHmassXrecoChic1Pre  -> DrawNormalized("same");
  myHmassXrecoChic2Pre  -> DrawNormalized("same");
  c9b.SaveAs("SummaryFromBHadMassXprefit.png");

  // Comparison between different decays from B0
  myHmassXreco2pre  -> SetTitle("");
  myHmassXreco2pre  -> GetXaxis()->SetTitle("m('X')");
  myHmassXreco1pre  -> SetLineColor(1);
  myHmassXreco2pre  -> SetLineColor(2);
  myHmassXreco3pre  -> SetLineColor(3);
  myHmassXreco4pre  -> SetLineColor(4);
  myHmassXreco2     -> SetTitle("");
  myHmassXreco2     -> GetXaxis()->SetTitle("m('X')");
  myHmassXreco1     -> SetLineColor(1);
  myHmassXreco2     -> SetLineColor(2);
  myHmassXreco3     -> SetLineColor(3);
  myHmassXreco4     -> SetLineColor(4);

  TCanvas c10a("c10a","",1);
  myHmassXreco1->Rebin();
  myHmassXreco2->Rebin();
  myHmassXreco3->Rebin();
  myHmassXreco4->Rebin();
  myHmassXreco2->DrawNormalized();
  myHmassXreco1->DrawNormalized("same");
  myHmassXreco3->DrawNormalized("same");
  myHmassXreco4->DrawNormalized("same");
  c10a.SaveAs("SummaryFromB0MassX.png");

  TCanvas c10b("c10b","",1);
  myHmassXreco1pre->Rebin();
  myHmassXreco2pre->Rebin();
  myHmassXreco3pre->Rebin();
  myHmassXreco4pre->Rebin();
  myHmassXreco2pre->DrawNormalized();
  myHmassXreco1pre->DrawNormalized("same");
  myHmassXreco3pre->DrawNormalized("same");
  myHmassXreco4pre->DrawNormalized("same");
  c10b.SaveAs("SummaryFromB0MassXprefit.png");
}

