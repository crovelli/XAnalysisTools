#ifndef BasicPlots_h
#define BasicPlots_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "TTree.h"  
#include "TCanvas.h"  
#include "TStyle.h"  
#include <iostream>   
#include <vector>
#include <string>
#include "TLorentzVector.h"  

#include "./XBase.h"

using namespace std;

class BasicPlots : public XBase{
public:

  //! constructor
  BasicPlots(TTree *tree=0);
  //! destructor
  virtual ~BasicPlots();
  //! loop over events
  void Loop();
  //! extra functions
  void PrepareOutputs(std::string filename);             
  void bookOutputTree();
  void genLevelInfo();
  void genMatch(TLorentzVector p4genmup, TLorentzVector p4genmum, TLorentzVector p4genpip, TLorentzVector p4genpim, TLorentzVector p4genks);

private:
  
  // Analysis methods
  // 
 
  // ---- outputs
  TFile* outFile_;
  TTree* outTree_;

  // dataset name
  std::string _datasetName;      

  //---common variables
  int genKsIdx;
  int genMuPIdx;
  int genMuMIdx;
  int genPiPIdx;
  int genPiMIdx;
  // 
  int recoMatchGenmuM;
  int recoMatchGenmuP;
  int recoMatchGenpiM;
  int recoMatchGenpiP;
  int recoMatchGenks;

  //---output tree branches variables
  float genMuPEta, genMuPPhi, genMuPPt;
  float genMuMEta, genMuMPhi, genMuMPt;
  float genPiPEta, genPiPPhi, genPiPPt;
  float genPiMEta, genPiMPhi, genPiMPt;
  float genKsEta,  genKsPhi,  genKsPt;
};

#endif
