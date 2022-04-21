#ifndef Trigger_h
#define Trigger_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>

#include "./XBase.h"

using namespace std;

class Trigger : public XBase{
public:

  //! constructor
  Trigger(TTree *tree=0);
  //! destructor
  virtual ~Trigger();
  //! loop over events
  void Loop();
  //! extra functions
  void genLevelInfo();

private:
  
  // Analysis methods
  // 

  //---common variables
  int genMuPIdx;
  int genMuMIdx;
  int genPiPIdx;
  int genPiMIdx;

  //---output tree branches variables
  float genMuPEta, genMuPPhi, genMuPPt;
  float genMuMEta, genMuMPhi, genMuMPt;
  float genPiPEta, genPiPPhi, genPiPPt;
  float genPiMEta, genPiMPhi, genPiMPt;
};

#endif
