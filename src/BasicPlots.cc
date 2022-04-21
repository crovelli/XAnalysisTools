#include "TMath.h"
#include "../include/BasicPlots.hh"    

using namespace std;

BasicPlots::BasicPlots(TTree *tree)     
    : XBase(tree) { 
  
}

BasicPlots::~BasicPlots() {

  // output
  outFile_->cd();
  outTree_->Write();
  outFile_->Close();
}     

void BasicPlots::Loop() {

  if (fChain == 0) return;

  // -------------------------------------------------------------------
  // Histos: mass checks
  TH1F *H_mJPsi_true   = new TH1F("H_mJPsi_true",   "H_mJPsi_true",   100, 2.7, 3.5);
  TH1F *H_mJPsi_comb   = new TH1F("H_mJPsi_comb",   "H_mJPsi_comb",   100, 2.7, 3.5);
  //
  TH1F *H_mRho_true    = new TH1F("H_mRho_true",   "H_mRho_true",   100, 0.3, 1.0);
  TH1F *H_mRho_truePF  = new TH1F("H_mRho_truePF", "H_mRho_truePF", 100, 0.3, 1.0);
  TH1F *H_mRho_comb    = new TH1F("H_mRho_comb",   "H_mRho_comb",   100, 0.3, 1.0);
  //
  TH1F *H_mKs_true     = new TH1F("H_mKs_true",   "H_mKs_true",   100, 0.40, 0.60);
  TH1F *H_mKs_truePF   = new TH1F("H_mKs_truePF", "H_mKs_truePF", 100, 0.40, 0.60);
  TH1F *H_mKs_comb     = new TH1F("H_mKs_comb",   "H_mKs_comb",   100, 0.40, 0.60);
  // 
  TH1F *H_mX_true      = new TH1F("H_mX_true",   "H_mX_true",   100, 3.6, 4.2);
  TH1F *H_mX_truePF    = new TH1F("H_mX_truePF", "H_mX_truePF", 100, 3.6, 4.2);
  TH1F *H_mX_comb      = new TH1F("H_mX_comb",   "H_mX_comb",   100, 3.6, 4.2);
  //
  TH1F *H_mB_true      = new TH1F("H_mB_true",   "H_mB_true",   100, 4.5, 6.);
  TH1F *H_mB_truePF    = new TH1F("H_mB_truePF", "H_mB_truePF", 100, 4.5, 6.);
  TH1F *H_mB_comb      = new TH1F("H_mB_comb",   "H_mB_comb",   100, 4.5, 6.);


  // -------------------------------------------------------------------  
  // Histos: trigger emulation
  TH1F *H_matchingTrackPt   = new TH1F("H_matchingTrackPt",  "H_matchingTrackPt",  50, 0., 10.);
  TH1F *H_matchingTrackD0   = new TH1F("H_matchingTrackD0",  "H_matchingTrackD0",  10, 0., 10.);
  TH1F *H_matchingMumudr    = new TH1F("H_matchingMumudr",   "H_matchingMumudr",   10, 0., 0.5); 
  TH1F *H_matchingJPsiPt    = new TH1F("H_matchingJPsiPt",   "H_matchingJPsiPt",   30, 0., 30.); 
  TH1F *H_matchingJPsiDca   = new TH1F("H_matchingJPsiDca",  "H_matchingJPsiDca",  20, 0., 0.5);  
  TH1F *H_matchingJPsiLxyS  = new TH1F("H_matchingJPsiLxyS", "H_matchingJPsiLxyS", 30, 0., 30.); 
  TH1F *H_matchingJPsiCosA  = new TH1F("H_matchingJPsiCosA", "H_matchingJPsiCosA", 20, -1., 1.);    
  TH1F *H_matchingJPsiProb  = new TH1F("H_matchingJPsiProb", "H_matchingJPsiProb", 10,  0., 1.);  

  TH1F *H_Bmass_trigger0_trueB   = new TH1F("H_Bmass_trigger0_trueB",   "H_Bmass_trigger0_trueB",   100, 4.5, 6.); 
  TH1F *H_Bmass_triggermu1_trueB = new TH1F("H_Bmass_triggermu1_trueB", "H_Bmass_triggermu1_trueB", 100, 4.5, 6.); 
  TH1F *H_Bmass_triggermu2_trueB = new TH1F("H_Bmass_triggermu2_trueB", "H_Bmass_triggermu2_trueB", 100, 4.5, 6.); 
  TH1F *H_Bmass_triggermu3_trueB = new TH1F("H_Bmass_triggermu3_trueB", "H_Bmass_triggermu3_trueB", 100, 4.5, 6.); 
  TH1F *H_Bmass_triggerpi1_trueB = new TH1F("H_Bmass_triggerpi1_trueB", "H_Bmass_triggerpi1_trueB", 100, 4.5, 6.); 
  TH1F *H_Bmass_triggerpi2_trueB = new TH1F("H_Bmass_triggerpi2_trueB", "H_Bmass_triggerpi2_trueB", 100, 4.5, 6.); 
  // 
  // ------------------------------------------------------------------- 


  // Loop over events
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;

  // for (Long64_t jentry=0; jentry<nentries;jentry++) {
  for (Long64_t jentry=0; jentry<5000;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%500==0) cout << jentry << endl;
    


    // ------------------------------------------------
    // Accessing gen level info
    genMuPIdx = -1;
    genMuMIdx = -1;
    genPiPIdx = -1;
    genPiMIdx = -1;
    genKsIdx  = -1;
    genLevelInfo();
    if (genMuPIdx<0 || genMuMIdx<0) cout << "problem with mc-truth, muons" << endl;
    if (genPiPIdx<0 || genPiMIdx<0) cout << "problem with mc-truth, pions" << endl;
    if (genKsIdx<0)                 cout << "problem with mc-truth, Ks"    << endl;

    TLorentzVector p4genmup, p4genmum;
    p4genmup.SetPtEtaPhiM(genMuPPt, genMuPEta, genMuPPhi, 0.1056);
    p4genmum.SetPtEtaPhiM(genMuMPt, genMuMEta, genMuMPhi, 0.1056);

    TLorentzVector p4genpip, p4genpim;
    p4genpip.SetPtEtaPhiM(genPiPPt, genPiPEta, genPiPPhi, 0.1396);
    p4genpim.SetPtEtaPhiM(genPiMPt, genPiMEta, genPiMPhi, 0.1396);

    TLorentzVector p4genks;
    p4genks.SetPtEtaPhiM(genKsPt, genKsEta, genKsPhi, 0.497611);

    // Reco - Gen match
    recoMatchGenmuM =-1;
    recoMatchGenmuP =-1;
    recoMatchGenpiM =-1;
    recoMatchGenpiP =-1;
    recoMatchGenks  =-1;
    genMatch(p4genmup, p4genmum, p4genpip, p4genpim, p4genks);


    // ------------------------------------------
    // Basic plots

    bool rhoTBF  = true;  // to fill just once per event the true distribution
    bool jpsiTBF = true;  
    bool ksTBF   = true;
    bool xTBF    = true;         
    bool bTBF    = true;

    for (int ii=0; ii<nB0; ii++) {

      int pi1idx = B0_pi1_idx[ii];
      int pi2idx = B0_pi2_idx[ii];
      int mu1idx = B0_mu1_idx[ii];
      int mu2idx = B0_mu2_idx[ii];
      int ksidx  = B0_k0short_idx[ii];

      bool trueRho  = (pi1idx==recoMatchGenpiM && pi2idx==recoMatchGenpiP) || (pi2idx==recoMatchGenpiM && pi1idx==recoMatchGenpiP);
      bool trueJPsi = (mu1idx==recoMatchGenmuM && mu2idx==recoMatchGenmuP) || (mu2idx==recoMatchGenmuM && mu1idx==recoMatchGenmuP);
      bool trueKs   = ksidx==recoMatchGenks;
      bool trueX    = trueRho && trueJPsi;
      bool trueB    = trueRho && trueJPsi && trueKs;
	
      TLorentzVector p4mu1;
      TLorentzVector p4mu2;
      p4mu1.SetPtEtaPhiM(B0_MuMu_prefit_mu1_pt[ii], B0_MuMu_prefit_mu1_eta[ii], B0_MuMu_prefit_mu1_phi[ii], 0.1056);
      p4mu2.SetPtEtaPhiM(B0_MuMu_prefit_mu2_pt[ii], B0_MuMu_prefit_mu2_eta[ii], B0_MuMu_prefit_mu2_phi[ii], 0.1056);

      TLorentzVector p4pi1;
      TLorentzVector p4pi2;
      p4pi1.SetPtEtaPhiM(B0_PiPi_prefit_pi1_pt[ii], B0_PiPi_prefit_pi1_eta[ii], B0_PiPi_prefit_pi1_phi[ii], 0.1396);
      p4pi2.SetPtEtaPhiM(B0_PiPi_prefit_pi2_pt[ii], B0_PiPi_prefit_pi2_eta[ii], B0_PiPi_prefit_pi2_phi[ii], 0.1396);

      TLorentzVector p4ks;
      p4ks.SetPtEtaPhiM(B0_K0s_mcFitted_pt[ii], B0_K0s_mcFitted_eta[ii], B0_K0s_mcFitted_phi[ii], 0.497611);

      float massJPsi = (p4mu1+p4mu2).M();
      float ptJPsi   = (p4mu1+p4mu2).Pt();
      float massRho  = (p4pi1+p4pi2).M();
      float massX    = (p4mu1+p4mu2+p4pi1+p4pi2).M();
      float massB    = (p4mu1+p4mu2+p4pi1+p4pi2+p4ks).M();

      if (trueRho) {
	if (rhoTBF) {
	  H_mRho_true->Fill(massRho);
	  H_mRho_truePF->Fill(B0_finalFit_Rho_mass[ii]);
	  rhoTBF=false;
	}
      }
      else H_mRho_comb->Fill(massRho);  

      if (trueJPsi) {
	if (jpsiTBF) {
	  H_mJPsi_true->Fill(massJPsi);
	  jpsiTBF=false;
	}
      }
      else { 
	H_mJPsi_comb->Fill(massJPsi);  
      }

      if (trueKs) {
	if (ksTBF) {
	  H_mKs_true->Fill(B0_K0s_prefit_mass[ii]);
	  H_mKs_truePF->Fill(B0_K0s_nmcFitted_mass[ii]);
	  ksTBF=false;
	}
      }
      else H_mKs_comb->Fill(B0_K0s_prefit_mass[ii]);  

      if (trueX) {
	if (xTBF) {  
	  H_mX_true->Fill(massX);
	  H_mX_truePF->Fill(B0_finalFit_X_mass[ii]);
	  xTBF=false;  
	}
      }
      else H_mX_comb->Fill(massX);  

      if (trueB) {
	if (bTBF) { 

	  // mass
	  H_mB_true->Fill(massB);
	  H_mB_truePF->Fill(B0_finalFit_mass[ii]);

	  // HLT checks for true candidates
	  int iHLT_DoubleMu4_JpsiTrk_Displaced= (int)HLT_DoubleMu4_JpsiTrk_Displaced;
	  bool triggerok0 = (bool)iHLT_DoubleMu4_JpsiTrk_Displaced;
	  if (triggerok0) H_Bmass_trigger0_trueB -> Fill(massB);

	  // Variables used at HLT level: muons
	  if (B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[ii]) H_matchingMumudr -> Fill(B0_MuMu_mu1_dr[ii]);
	  if (B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[ii]) H_matchingMumudr -> Fill(B0_MuMu_mu2_dr[ii]);
	  if (B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[ii] && B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[ii]) {
	    H_matchingJPsiPt   -> Fill(ptJPsi);
	    H_matchingJPsiDca  -> Fill(B0_MuMu_DCA[ii]);
	    H_matchingJPsiLxyS -> Fill(B0_MuMu_LxySign[ii]);
	    H_matchingJPsiCosA -> Fill(B0_MuMu_cosAlpha[ii]);
	    H_matchingJPsiProb -> Fill(B0_MuMu_sv_prob[ii]);
	  }

	  // HLT checks: muons part
	  bool mumuok1 = false;
	  if (B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[ii] && B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[ii]) mumuok1 = true;
	  if (triggerok0 && mumuok1) H_Bmass_triggermu1_trueB -> Fill(massB);
	  
	  bool mumuok2 = true;
	  if (B0_MuMu_prefit_mu1_pt[ii]<4 || fabs(B0_MuMu_prefit_mu1_eta[ii])>2.5) mumuok2 = false;
	  if (B0_MuMu_prefit_mu2_pt[ii]<4 || fabs(B0_MuMu_prefit_mu2_eta[ii])>2.5) mumuok2 = false;
	  if (B0_MuMu_mu1_dr[ii]>2) mumuok2 = false; 
	  if (B0_MuMu_mu2_dr[ii]>2) mumuok2 = false; 
	  if (massJPsi>3.3 || massJPsi<2.9) mumuok2 = false; 
	  if (ptJPsi<6.9) mumuok2 = false; 
	  if (B0_MuMu_DCA[ii]>0.5) mumuok2 = false; 
	  if (triggerok0 && mumuok1 && mumuok2) H_Bmass_triggermu2_trueB -> Fill(massB);  

	  bool mumuok3 = true;
	  if (B0_MuMu_LxySign[ii]<3)    mumuok3 = false;
	  if (B0_MuMu_cosAlpha[ii]<0.9) mumuok3 = false;  
	  if (B0_MuMu_sv_prob[ii]<0.1)  mumuok3 = false;  
	  if (triggerok0 && mumuok1 && mumuok2 && mumuok3) H_Bmass_triggermu3_trueB -> Fill(massB); 

	  // HLT checks: tracks part
	  bool pipiok1 = false;
	  if (B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[ii] || B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[ii] || B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[ii] || B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[ii]) pipiok1 = true;
	  if (triggerok0 && pipiok1) H_Bmass_triggerpi1_trueB -> Fill(massB); 

	  bool pipiok2 = false;
	  if (B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[ii] &&
	      B0_PiPi_prefit_pi1_pt[ii]>1.2 &&
	      fabs(B0_PiPi_prefit_pi1_eta[ii])<2.5 &&
	      B0_PiPi_pi1_d0sig[ii]>2) pipiok2 = true; 
	  if (B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[ii] &&
	      B0_PiPi_prefit_pi2_pt[ii]>1.2 &&
	      fabs(B0_PiPi_prefit_pi2_eta[ii])<2.5 &&
	      B0_PiPi_pi2_d0sig[ii]>2) pipiok2 = true; 
	  if (B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[ii] &&
	      B0_K0s_matchTrack1_pt[ii]>1.2 &&
	      fabs(B0_K0s_matchTrack1_eta[ii])<2.5 &&
	      B0_K0s_matchTrack1_D0sign[ii]>2) pipiok2 = true; 
	  if (B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[ii] &&
	      B0_K0s_matchTrack2_pt[ii]>1.2 &&
	      fabs(B0_K0s_matchTrack2_eta[ii])<2.5 &&
	      B0_K0s_matchTrack2_D0sign[ii]>2) pipiok2 = true; 
	  if (triggerok0 && pipiok1 && pipiok2) H_Bmass_triggerpi2_trueB -> Fill(massB); 

	  // Variables used at HLT level: tracks
	  if (B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[ii]) {
	    H_matchingTrackPt->Fill(B0_PiPi_prefit_pi1_pt[ii]);
	    H_matchingTrackD0->Fill(B0_PiPi_pi1_d0sig[ii]); 
	  }
	  if (B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[ii]) {
	    H_matchingTrackPt->Fill(B0_PiPi_prefit_pi2_pt[ii]);
	    H_matchingTrackD0->Fill(B0_PiPi_pi2_d0sig[ii]);  
	  }
	  if (B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[ii]) {
	    H_matchingTrackPt->Fill(B0_K0s_matchTrack1_pt[ii]);
	    H_matchingTrackD0->Fill(B0_K0s_matchTrack1_D0sign[ii]);
	  }
	  if (B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[ii]) {
	    H_matchingTrackPt->Fill(B0_K0s_matchTrack2_pt[ii]);  
	    H_matchingTrackD0->Fill(B0_K0s_matchTrack2_D0sign[ii]);
	  }
 
	  bTBF=false;  
	}
      }
      else { 
	H_mB_comb->Fill(massB);  
      }

    } // Loop over Bs

    // Filling the output tree
    outTree_->Fill();

  } // Loop over events

  
  // Plots
  gStyle->SetOptStat(0000);

  TCanvas c1a("c1a","",1);
  H_mJPsi_true -> SetLineWidth(2);
  H_mJPsi_true -> SetLineColor(2);
  H_mJPsi_comb -> SetLineWidth(2);
  H_mJPsi_comb -> SetLineColor(4);
  H_mJPsi_true -> Draw(); 
  H_mJPsi_comb -> Draw("same"); 
  c1a.SaveAs("jpsi_trueVsComb.png");

  TCanvas c1b("c1b","",1);
  H_mRho_true -> SetLineWidth(2);
  H_mRho_true -> SetLineColor(2);
  H_mRho_comb -> SetLineWidth(2);
  H_mRho_comb -> SetLineColor(4);
  H_mRho_true -> Draw(); 
  H_mRho_comb -> Draw("same"); 
  c1b.SaveAs("rho_trueVsComb.png");

  TCanvas c1c("c1c","",1);
  H_mKs_true -> SetLineWidth(2);
  H_mKs_true -> SetLineColor(2);
  H_mKs_comb -> SetLineWidth(2);
  H_mKs_comb -> SetLineColor(4);
  H_mKs_true -> Draw(); 
  H_mKs_comb -> Draw("same"); 
  c1c.SaveAs("ks_trueVsComb.png");

  TCanvas c1d("c1d","",1);
  H_mX_true -> SetLineWidth(2);
  H_mX_true -> SetLineColor(2);
  H_mX_comb -> SetLineWidth(2);
  H_mX_comb -> SetLineColor(4);
  H_mX_true -> Draw(); 
  H_mX_comb -> Draw("same"); 
  c1d.SaveAs("X_trueVsComb.png");

  TCanvas c1e("c1e","",1);
  H_mB_true -> SetLineWidth(2);
  H_mB_true -> SetLineColor(2);
  H_mB_comb -> SetLineWidth(2);
  H_mB_comb -> SetLineColor(4);
  H_mB_true -> Draw(); 
  H_mB_comb -> Draw("same"); 
  c1e.SaveAs("B_trueVsComb.png");

  // ----------------------------------------------

  TCanvas c2a("c2a","",1);
  H_mKs_true   -> SetLineWidth(2);
  H_mKs_true   -> SetLineColor(2);
  H_mKs_truePF -> SetLineWidth(2);
  H_mKs_truePF -> SetLineColor(7);
  H_mKs_truePF -> DrawNormalized(); 
  H_mKs_true   -> DrawNormalized("same"); 
  c2a.SaveAs("ksTrue_fitEffect.png");

  TCanvas c2b("c2b","",1);
  H_mRho_true   -> SetLineWidth(2);
  H_mRho_true   -> SetLineColor(2);
  H_mRho_truePF -> SetLineWidth(2);
  H_mRho_truePF -> SetLineColor(7);
  H_mRho_truePF -> DrawNormalized(); 
  H_mRho_true   -> DrawNormalized("same"); 
  c2b.SaveAs("rhoTrue_fitEffect.png");

  TCanvas c2d("c2d","",1);
  H_mX_true   -> SetLineWidth(2);
  H_mX_true   -> SetLineColor(2);
  H_mX_truePF -> SetLineWidth(2);
  H_mX_truePF -> SetLineColor(7);
  H_mX_truePF -> DrawNormalized(); 
  H_mX_true   -> DrawNormalized("same"); 
  c2d.SaveAs("xTrue_fitEffect.png");

  TCanvas c2e("c2e","",1);
  H_mB_true   -> SetLineWidth(2);
  H_mB_true   -> SetLineColor(2);
  H_mB_truePF -> SetLineWidth(2);
  H_mB_truePF -> SetLineColor(7);
  H_mB_truePF -> DrawNormalized(); 
  H_mB_true   -> DrawNormalized("same"); 
  c2e.SaveAs("bTrue_fitEffect.png");

  // ----------------------------------------
  TCanvas c4a("c4a","",1);
  H_Bmass_trigger0_trueB   -> SetLineColor(1);
  H_Bmass_triggermu1_trueB -> SetLineColor(2);
  H_Bmass_triggermu2_trueB -> SetLineColor(3);
  H_Bmass_triggermu3_trueB -> SetLineColor(7);
  H_Bmass_trigger0_trueB   -> Draw();
  H_Bmass_triggermu1_trueB -> Draw("same"); 
  H_Bmass_triggermu2_trueB -> Draw("same"); 
  H_Bmass_triggermu3_trueB -> Draw("same"); 
  c4a.SaveAs("hltEmulation_muons.png");

  TCanvas c4ab("c4ab","",1);
  H_Bmass_trigger0_trueB   -> SetLineColor(1);
  H_Bmass_triggerpi1_trueB -> SetLineColor(2);
  H_Bmass_triggerpi2_trueB -> SetLineColor(7);
  H_Bmass_trigger0_trueB   -> Draw();
  H_Bmass_triggerpi1_trueB -> Draw("same"); 
  H_Bmass_triggerpi2_trueB -> Draw("same"); 
  c4ab.SaveAs("hltEmulation_pions.png");

  TCanvas c4ac("c4ac","",1);
  H_matchingTrackPt -> SetLineColor(1);
  H_matchingTrackPt -> SetLineColor(2);
  H_matchingTrackPt -> Draw();
  c4ac.SaveAs("hltEmulation_matchingPionsPt.png");

  TCanvas c4ad("c4ad","",1);
  H_matchingTrackD0 -> SetLineColor(1);
  H_matchingTrackD0 -> SetLineColor(2);
  H_matchingTrackD0 -> Draw();
  c4ad.SaveAs("hltEmulation_matchingPionsD0.png");

  TCanvas c4b("c4b","",1);
  H_matchingMumudr -> SetLineWidth(2);
  H_matchingMumudr -> SetLineColor(2);
  H_matchingMumudr -> DrawNormalized(); 
  c4b.SaveAs("hltEmulation_matchingMuonsDr.png");

  TCanvas c4c("c4c","",1);
  H_matchingJPsiPt -> SetLineWidth(2);
  H_matchingJPsiPt -> SetLineColor(2);
  H_matchingJPsiPt -> DrawNormalized(); 
  c4c.SaveAs("hltEmulation_matchingJPsiPt.png");

  TCanvas c4d("c4d","",1);
  H_matchingJPsiDca -> SetLineWidth(2);
  H_matchingJPsiDca -> SetLineColor(2);
  H_matchingJPsiDca -> DrawNormalized(); 
  c4d.SaveAs("hltEmulation_matchingJPsiDca.png");

  TCanvas c4e("c4e","",1);
  H_matchingJPsiLxyS -> SetLineWidth(2);
  H_matchingJPsiLxyS -> SetLineColor(2);
  H_matchingJPsiLxyS -> DrawNormalized("same"); 
  c4e.SaveAs("hltEmulation_matchingJPsiLxySign.png");

  TCanvas c4f("c4f","",1);
  H_matchingJPsiCosA -> SetLineWidth(2);
  H_matchingJPsiCosA -> SetLineColor(2);
  H_matchingJPsiCosA -> DrawNormalized(); 
  c4f.SetLogy();
  c4f.SaveAs("hltEmulation_matchingJPsiCosAlpha.png");   

  TCanvas c4g("c4g","",1);
  H_matchingJPsiProb -> SetLineWidth(2); 
  H_matchingJPsiProb -> SetLineColor(2); 
  H_matchingJPsiProb -> DrawNormalized(); 
  c4g.SaveAs("hltEmulation_matchingJPsiSvProb.png");
}


void BasicPlots::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");
  
  bookOutputTree();
};


void BasicPlots::bookOutputTree() 
{
  outTree_ = new TTree("Xtree", "Xtree");
  
  cout << "Booking tree" << endl;
  outTree_->Branch("genMuPEta", &genMuPEta, "genMuPEta/F");    
  outTree_->Branch("genMuPPhi", &genMuPPhi, "genMuPPhi/F");    
  outTree_->Branch("genMuPPt",  &genMuPPt,  "genMuPPt/F");    
  //
  outTree_->Branch("genMuMEta", &genMuMEta, "genMuMEta/F");    
  outTree_->Branch("genMuMPhi", &genMuMPhi, "genMuMPhi/F");    
  outTree_->Branch("genMuMPt",  &genMuMPt,  "genMuMPt/F");    
  // 
  outTree_->Branch("genPiPEta", &genPiPEta, "genPiPEta/F");    
  outTree_->Branch("genPiPPhi", &genPiPPhi, "genPiPPhi/F");    
  outTree_->Branch("genPiPPt",  &genPiPPt,  "genPiPPt/F");    
  //
  outTree_->Branch("genPiMEta", &genPiMEta, "genPiMEta/F");    
  outTree_->Branch("genPiMPhi", &genPiMPhi, "genPiMPhi/F");    
  outTree_->Branch("genPiMPt",  &genPiMPt,  "genPiMPt/F");    
}

void BasicPlots::genLevelInfo()
{
  for (int igen=0; igen<nGenPart; igen++) { 
    if (GenPart_genPartIdxMother[igen]<0) continue;
    int pdgId    = GenPart_pdgId[igen];
    int mumIdx   = GenPart_genPartIdxMother[igen];
    int pdgIdMum = GenPart_pdgId[mumIdx];
    if (pdgId==13   && pdgIdMum==443) genMuPIdx = igen;
    if (pdgId==-13  && pdgIdMum==443) genMuMIdx = igen;
    if (pdgId==211  && pdgIdMum==113) genPiPIdx = igen;  
    if (pdgId==-211 && pdgIdMum==113) genPiMIdx = igen;  
    if (abs(pdgId)==310  && abs(pdgIdMum)==511) genKsIdx = igen;  
  }
  if (genMuPIdx<0 || genMuMIdx<0 || genPiPIdx<0 || genPiMIdx<0) {
    cout << "There is a problem with gen level" << endl;
    exit; 
  }

  genMuPEta = GenPart_eta[genMuPIdx];
  genMuPPhi = GenPart_phi[genMuPIdx];
  genMuPPt  = GenPart_pt[genMuPIdx];
  genMuMEta = GenPart_eta[genMuMIdx];
  genMuMPhi = GenPart_phi[genMuMIdx];
  genMuMPt  = GenPart_pt[genMuMIdx];
  genPiPEta = GenPart_eta[genPiPIdx];
  genPiPPhi = GenPart_phi[genPiPIdx];
  genPiPPt  = GenPart_pt[genPiPIdx];
  genPiMEta = GenPart_eta[genPiMIdx];
  genPiMPhi = GenPart_phi[genPiMIdx];
  genPiMPt  = GenPart_pt[genPiMIdx];
  genKsEta  = GenPart_eta[genKsIdx];
  genKsPhi  = GenPart_phi[genKsIdx];
  genKsPt   = GenPart_pt[genKsIdx];
} 

void BasicPlots::genMatch(TLorentzVector p4genmup, TLorentzVector p4genmum, TLorentzVector p4genpip, TLorentzVector p4genpim, TLorentzVector p4genks)
{
  float minDrGenmuM = 100;
  float minDrGenmuP = 100;
  float minDrGenpiM = 100; 
  float minDrGenpiP = 100; 
  float minDrGenks  = 100;

  for (int ii=0; ii<nB0; ii++) {

    TLorentzVector p4recomu1, p4recomu2;
    p4recomu1.SetPtEtaPhiM(B0_MuMu_prefit_mu1_pt[ii], B0_MuMu_prefit_mu1_eta[ii], B0_MuMu_prefit_mu1_phi[ii], 0.1056);
    p4recomu2.SetPtEtaPhiM(B0_MuMu_prefit_mu2_pt[ii], B0_MuMu_prefit_mu2_eta[ii], B0_MuMu_prefit_mu2_phi[ii], 0.1056);

    TLorentzVector p4recopi1, p4recopi2;
    p4recopi1.SetPtEtaPhiM(B0_PiPi_prefit_pi1_pt[ii], B0_PiPi_prefit_pi1_eta[ii], B0_PiPi_prefit_pi1_phi[ii], 0.1396);
    p4recopi2.SetPtEtaPhiM(B0_PiPi_prefit_pi2_pt[ii], B0_PiPi_prefit_pi2_eta[ii], B0_PiPi_prefit_pi2_phi[ii], 0.1396);

    TLorentzVector p4recoks;
    p4recoks.SetPtEtaPhiM(B0_K0s_mcFitted_pt[ii], B0_K0s_mcFitted_eta[ii], B0_K0s_mcFitted_phi[ii], 0.497611);


    // ------------------------------------------------
    // muons
    float dRmu_genp_reco1 = p4recomu1.DeltaR(p4genmup);
    float dRmu_genp_reco2 = p4recomu2.DeltaR(p4genmup);
    float dRmu_genm_reco1 = p4recomu1.DeltaR(p4genmum);
    float dRmu_genm_reco2 = p4recomu2.DeltaR(p4genmum);

    if (dRmu_genp_reco1<minDrGenmuP && dRmu_genp_reco1<dRmu_genp_reco2 && dRmu_genp_reco1<0.03) {
      minDrGenmuP = dRmu_genp_reco1;
      recoMatchGenmuP = B0_mu1_idx[ii];
    }
    if (dRmu_genp_reco2<minDrGenmuP && dRmu_genp_reco2<dRmu_genp_reco1 && dRmu_genp_reco2<0.03) {
      minDrGenmuP = dRmu_genp_reco2;
      recoMatchGenmuP = B0_mu2_idx[ii];
    }

    if (dRmu_genm_reco1<minDrGenmuM && dRmu_genm_reco1<dRmu_genm_reco2 && dRmu_genm_reco1<0.03) {
      minDrGenmuM = dRmu_genm_reco1;
      recoMatchGenmuM = B0_mu1_idx[ii];
    }
    if (dRmu_genm_reco2<minDrGenmuM && dRmu_genm_reco2<dRmu_genm_reco1 && dRmu_genm_reco2<0.03) {
      minDrGenmuM = dRmu_genm_reco2;
      recoMatchGenmuM = B0_mu2_idx[ii];
    }


    // ------------------------------------------------
    // pions
    float dRpi_genp_reco1 = p4recopi1.DeltaR(p4genpip);
    float dRpi_genp_reco2 = p4recopi2.DeltaR(p4genpip);
    float dRpi_genm_reco1 = p4recopi1.DeltaR(p4genpim);
    float dRpi_genm_reco2 = p4recopi2.DeltaR(p4genpim);

    if (dRpi_genp_reco1<minDrGenpiP && dRpi_genp_reco1<dRpi_genp_reco2 && dRpi_genp_reco1<0.03) {
      minDrGenpiP = dRpi_genp_reco1;
      recoMatchGenpiP = B0_pi1_idx[ii];
    }
    if (dRpi_genp_reco2<minDrGenpiP && dRpi_genp_reco2<dRpi_genp_reco1 && dRpi_genp_reco2<0.03) {
      minDrGenpiP = dRpi_genp_reco2;
      recoMatchGenpiP = B0_pi2_idx[ii];
    }

    if (dRpi_genm_reco1<minDrGenpiM && dRpi_genm_reco1<dRpi_genm_reco2 && dRpi_genm_reco1<0.03) {
      minDrGenpiM = dRpi_genm_reco1;
      recoMatchGenpiM = B0_pi1_idx[ii];
    }
    if (dRpi_genm_reco2<minDrGenpiM && dRpi_genm_reco2<dRpi_genm_reco1 && dRpi_genm_reco2<0.03) {
      minDrGenpiM = dRpi_genm_reco2;
      recoMatchGenpiM = B0_pi2_idx[ii];
    }


    // ------------------------------------------------
    // ks
    float dRks_gen_reco = p4recoks.DeltaR(p4genks);
    if (dRks_gen_reco<minDrGenks && dRks_gen_reco<0.03) {
      minDrGenks = dRks_gen_reco;
      recoMatchGenks = B0_k0short_idx[ii];
    }

  }  // Loop over Bs
}
