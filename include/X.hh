/// This X class is an auxiliary class which contains basic
/// functionality useful for any analysis of B-Parked events.
/// It derives from XBase.

#ifndef X_h
#define X_h

#include "./XBase.h"

#include <TLorentzVector.h>
#include <string>
#include <vector>
#include <map>

class X : public XBase{

public:

  X(TTree *tree=0);
  virtual ~X();

private:

protected:

};

#endif
