// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "../include/Application.hh"
#if Application == 1
#include "./BasicPlots.cc"
#endif
#if Application == 2
#include "./Trigger.cc"
#endif
#if Application == 3
#include "./Selection.cc"
#endif
#if Application == 4
#include "./GenLevelStudy.cc"
#endif

int main(int argc, char* argv[]) {

  char inputFileName[500];
  char outputFileName[500];
  char dataset[150];
  if ( argc < 2 ){
    std::cout << "missing argument: insert at least inputFile with list of root files" << std::endl; 
    std::cout << "XApp inputFile [outputFile] [1=MC,0=data] [dataset]" << std::endl;
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  if (argc <= 2 ) strcpy(outputFileName,argv[1]);
  else strcpy(outputFileName,argv[2]);
  int isMC=1;
  if(argc==5) {
    isMC=atoi(argv[3]);
    strcpy(dataset,argv[4]);
  }
  
  // -------------------------
  // Loading the file
  TChain *theChain = new TChain("Events");
  char Buffer[5000];
  char MyRootFile[10000];
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);
  char tmpFileName[256];
  vector<string> filesToRemove;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
        sscanf(Buffer,"%s",MyRootFile);
        theChain->Add(TString(MyRootFile));
        std::cout << "chaining " << MyRootFile << std::endl;
      }
  }
  inputFile->close();
  delete inputFile;
  
  
#if Application == 1
  
  BasicPlots plots(theChain);
  
  plots.PrepareOutputs(outputFileName);   
  
  plots.Loop();

#endif

#if Application == 2

  Trigger plots(theChain);

  plots.Loop();

#endif

#if Application == 3

  Selection plots(theChain);

  plots.Loop();

#endif

#if Application == 4

  GenLevelStudy plots(theChain);

  plots.Loop();

#endif

  return 0;

}
