#include "TMath.h"
#include "TTree.h"

#include <iostream> 
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "../include/Trigger.hh"    

using namespace std;

Trigger::Trigger(TTree *tree)     
    : XBase(tree) { 
}

Trigger::~Trigger() { }     

void Trigger::Loop() {

  if (fChain == 0) return;

  // Histos
  TH1F *H_triggers = new TH1F("H_triggers","H_trigger", 16, 0.5, 16.5);
  //
  TH1F *H_LMuPt_mybit_on  = new TH1F("H_LMuPt_mybit_on",  "H_LMuPt_mybit_on",  20, 0., 20.);
  TH1F *H_LMuPt_mybit_off = new TH1F("H_LMuPt_mybit_off", "H_LMuPt_mybit_off", 20, 0., 20.);
  TH1F *H_SMuPt_mybit_on  = new TH1F("H_SMuPt_mybit_on",  "H_SMuPt_mybit_on",  15, 0., 15.);
  TH1F *H_SMuPt_mybit_off = new TH1F("H_SMuPt_mybit_off", "H_SMuPt_mybit_off", 15, 0., 15.);
  TH1F *H_LPiPt_mybit_on  = new TH1F("H_LPiPt_mybit_on",  "H_LPiPt_mybit_on",  30, 0., 6.);
  TH1F *H_LPiPt_mybit_off = new TH1F("H_LPiPt_mybit_off", "H_LPiPt_mybit_off", 30, 0., 6.);
  TH1F *H_SPiPt_mybit_on  = new TH1F("H_SPiPt_mybit_on",  "H_SPiPt_mybit_on",  20, 0., 4.);
  TH1F *H_SPiPt_mybit_off = new TH1F("H_SPiPt_mybit_off", "H_SPiPt_mybit_off", 20, 0., 4.);




  // Loop over events
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%500==0) cout << jentry << endl;

    // Accessing gen level info
    genLevelInfo();    
    float leadingMuPt    = genMuPPt;
    float subleadingMuPt = genMuMPt;
    if (genMuPPt<genMuMPt) {
      leadingMuPt    = genMuMPt;
      subleadingMuPt = genMuPPt;
    }
    float leadingPiPt    = genPiPPt;
    float subleadingPiPt = genPiMPt;
    if (genPiPPt<genPiMPt) {
      leadingPiPt    = genPiMPt;
      subleadingPiPt = genPiPPt;
    }

    // Fired triggers
    int iHLT_Dimuon25_Jpsi = (int)HLT_Dimuon25_Jpsi;
    int iHLT_Dimuon20_Jpsi_Barrel_Seagulls= (int)HLT_Dimuon20_Jpsi_Barrel_Seagulls;
    int iHLT_DoubleMu4_JpsiTrk_Displaced= (int)HLT_DoubleMu4_JpsiTrk_Displaced;
    int iHLT_DoubleMu4_JpsiTrkTrk_Displaced= (int)HLT_DoubleMu4_JpsiTrkTrk_Displaced;
    int iHLT_DoubleMu4_3_Jpsi_Displaced= (int)HLT_DoubleMu4_3_Jpsi_Displaced;
    int iHLT_DoubleMu4_3_Jpsi= (int)HLT_DoubleMu4_3_Jpsi;
    int iHLT_DoubleMu4_Jpsi_Displaced= (int)HLT_DoubleMu4_Jpsi_Displaced;
    int iHLT_Dimuon18_PsiPrime= (int)HLT_Dimuon18_PsiPrime;
    int iHLT_Dimuon10_PsiPrime_Barrel_Seagulls= (int)HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
    int iHLT_DoubleMu4_PsiPrimeTrk_Displaced= (int)HLT_DoubleMu4_PsiPrimeTrk_Displaced;
    int iHLT_Dimuon0_Jpsi3p5_Muon2= (int)HLT_Dimuon0_Jpsi3p5_Muon2;
    int iHLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi= (int)HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;
    int iHLT_DoubleMu2_Jpsi_DoubleTrk1_Phi= (int)HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;
    int iHLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05 = (int)HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
    int iHLT_DoubleMu4_3_Bs = (int)HLT_DoubleMu4_3_Bs;
    //
    H_triggers -> Fill(1); 
    if (iHLT_Dimuon25_Jpsi==1) H_triggers -> Fill(2);
    if (iHLT_Dimuon20_Jpsi_Barrel_Seagulls==1) H_triggers -> Fill(3);
    if (iHLT_DoubleMu4_JpsiTrk_Displaced==1) H_triggers -> Fill(4);                //
    if (iHLT_DoubleMu4_JpsiTrkTrk_Displaced==1) H_triggers -> Fill(5);
    if (iHLT_DoubleMu4_3_Jpsi_Displaced==1) H_triggers -> Fill(6);                 //
    if (iHLT_DoubleMu4_3_Jpsi==1) H_triggers -> Fill(7);                           //
    if (iHLT_DoubleMu4_Jpsi_Displaced==1) H_triggers -> Fill(8);
    if (iHLT_Dimuon18_PsiPrime==1) H_triggers -> Fill(9);
    if (iHLT_Dimuon10_PsiPrime_Barrel_Seagulls==1) H_triggers -> Fill(10);
    if (iHLT_DoubleMu4_PsiPrimeTrk_Displaced==1) H_triggers -> Fill(11);
    if (iHLT_Dimuon0_Jpsi3p5_Muon2==1) H_triggers -> Fill(12);
    if (iHLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi==1) H_triggers -> Fill(13);
    if (iHLT_DoubleMu2_Jpsi_DoubleTrk1_Phi==1) H_triggers -> Fill(14);
    if (iHLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05==1) H_triggers -> Fill(15);
    if (iHLT_DoubleMu4_3_Bs==1) H_triggers -> Fill(16);


    if (iHLT_DoubleMu4_JpsiTrk_Displaced==1) {
      H_LMuPt_mybit_on -> Fill(leadingMuPt);
      H_SMuPt_mybit_on -> Fill(subleadingMuPt);
      H_LPiPt_mybit_on -> Fill(leadingPiPt);
      H_SPiPt_mybit_on -> Fill(subleadingPiPt);
    } else {
      H_LMuPt_mybit_off -> Fill(leadingMuPt);
      H_SMuPt_mybit_off -> Fill(subleadingMuPt);
      H_LPiPt_mybit_off -> Fill(leadingPiPt);
      H_SPiPt_mybit_off -> Fill(subleadingPiPt);
    }
    
  } // Loop over events

  
  // Plots
  gStyle->SetOptStat(0000);

  TLegend leg(0.65,0.70,0.95,0.85);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(H_LMuPt_mybit_on, "Fired");
  leg.AddEntry(H_LMuPt_mybit_off,"Not-Fired");

  TCanvas c1("c1","",1);
  H_triggers -> SetLineWidth(2);
  H_triggers -> SetLineColor(2);
  H_triggers -> GetXaxis()->SetTitle("bit");
  H_triggers -> Draw(); 
  c1.SaveAs("triggers.png");

  TCanvas c2a("c2a","DoubleMu4_JpsiTrk_Displaced",1);
  H_LMuPt_mybit_on  -> GetXaxis()->SetTitle("Leading #mu p_{T} [GeV]"); 
  H_LMuPt_mybit_off -> GetXaxis()->SetTitle("Leading #mu p_{T} [GeV]"); 
  H_LMuPt_mybit_on  -> SetTitle("");
  H_LMuPt_mybit_off -> SetTitle(""); 
  H_LMuPt_mybit_on  -> SetLineWidth(2);
  H_LMuPt_mybit_on  -> SetLineColor(2);
  H_LMuPt_mybit_off -> SetLineWidth(2);
  H_LMuPt_mybit_off -> SetLineColor(4);
  H_LMuPt_mybit_off -> DrawNormalized();
  H_LMuPt_mybit_on  -> DrawNormalized("same"); 
  leg.Draw();
  c2a.SaveAs("LeadingGenMuPt_DoubleMu4_JpsiTrk_Displaced.png");

  TCanvas c2b("c2b","DoubleMu4_JpsiTrk_Displaced",1);
  H_SMuPt_mybit_on  -> GetXaxis()->SetTitle("Subleading #mu p_{T} [GeV]"); 
  H_SMuPt_mybit_off -> GetXaxis()->SetTitle("Subleading #mu p_{T} [GeV]"); 
  H_SMuPt_mybit_on  -> SetTitle("");
  H_SMuPt_mybit_off -> SetTitle(""); 
  H_SMuPt_mybit_on  -> SetLineWidth(2);
  H_SMuPt_mybit_on  -> SetLineColor(2);
  H_SMuPt_mybit_off -> SetLineWidth(2);
  H_SMuPt_mybit_off -> SetLineColor(4);
  H_SMuPt_mybit_off -> DrawNormalized();
  H_SMuPt_mybit_on  -> DrawNormalized("same"); 
  leg.Draw();
  c2b.SaveAs("SubLeadingGenMuPt_DoubleMu4_JpsiTrk_Displaced.png");

  // --------------------------------------------

  TCanvas c12a("c12a","DoubleMu4_JpsiTrk_Displaced",1);
  H_LPiPt_mybit_on  -> GetXaxis()->SetTitle("Leading #pi p_{T} [GeV]"); 
  H_LPiPt_mybit_off -> GetXaxis()->SetTitle("Leading #pi p_{T} [GeV]"); 
  H_LPiPt_mybit_on  -> SetTitle("");
  H_LPiPt_mybit_off -> SetTitle(""); 
  H_LPiPt_mybit_on  -> SetLineWidth(2);
  H_LPiPt_mybit_on  -> SetLineColor(2);
  H_LPiPt_mybit_off -> SetLineWidth(2);
  H_LPiPt_mybit_off -> SetLineColor(4);
  H_LPiPt_mybit_off -> DrawNormalized();
  H_LPiPt_mybit_on  -> DrawNormalized("same"); 
  leg.Draw();
  c12a.SaveAs("LeadingGenPiPt_DoubleMu4_JpsiTrk_Displaced.png");

  TCanvas c12b("c12b","DoubleMu4_JpsiTrk_Displaced",1);
  H_SPiPt_mybit_on  -> GetXaxis()->SetTitle("Subleading #pi p_{T} [GeV]"); 
  H_SPiPt_mybit_off -> GetXaxis()->SetTitle("Subleading #pi p_{T} [GeV]"); 
  H_SPiPt_mybit_on  -> SetTitle("");
  H_SPiPt_mybit_off -> SetTitle(""); 
  H_SPiPt_mybit_on  -> SetLineWidth(2);
  H_SPiPt_mybit_on  -> SetLineColor(2);
  H_SPiPt_mybit_off -> SetLineWidth(2);
  H_SPiPt_mybit_off -> SetLineColor(4);
  H_SPiPt_mybit_off -> DrawNormalized();
  H_SPiPt_mybit_on  -> DrawNormalized("same"); 
  leg.Draw();
  c12b.SaveAs("SubLeadingGenPiPt_DoubleMu4_JpsiTrk_Displaced.png");
}



void Trigger::genLevelInfo()
{
  for (int igen=0; igen<nGenPart; igen++) { 
    if (GenPart_genPartIdxMother[igen]<0) continue;
    int pdgId = GenPart_pdgId[igen];
    int mumIdx = GenPart_genPartIdxMother[igen];
    int pdgIdMum = GenPart_pdgId[mumIdx];
    if (pdgId==13   && pdgIdMum==443) genMuPIdx = igen;
    if (pdgId==-13  && pdgIdMum==443) genMuMIdx = igen;
    if (pdgId==211  && pdgIdMum==113) genPiPIdx = igen;  
    if (pdgId==-211 && pdgIdMum==113) genPiMIdx = igen;  
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
} 
