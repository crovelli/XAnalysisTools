#ifndef GenLevelStudy_h
#define GenLevelStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>

#include "./XBase.h"

using namespace std;

class GenLevelStudy : public XBase{
public:

  //! constructor
  GenLevelStudy(TTree *tree=0);
  //! destructor
  virtual ~GenLevelStudy();
  //! loop over events
  void Loop();
  //! extra functions
  void genLevelInfo();

private:
  
  // Analysis methods
  bool contains(vector<int> vec, const int & elem);

  //---common variables

  //---output tree branches variables
};

#endif
