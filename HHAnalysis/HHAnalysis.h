#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <list>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TROOT.h"

#include "PlotFunctions/MapBranches.h"

class HHAnalysis
{
 public:
  HHAnalysis();
  ~HHAnalysis();


  float GetExpectedSignificance(TH1 *histBkg, TH1 *histSignal, int mode, int method, double maxUnc, double *error);

  double getExpectedSignificance (TH1 *hBkg, TH1 *hSig, int mode,int method, double maxUnc,double *error);

  /**\brief Creates an histogram for each category of btag, each sample and each variable of interest. Saves the TH1D* of output in a root file.
   */
  void CreateSaveDistri(std::vector<std::string> vectInFiles, int selectionType, std::vector<std::string> vectVariables, std::vector<int> vectCategories, std::string outputFileName);

  /**\brief Returns the weight according to the weightType:
    - 0: for no selection weightMC*weightvertex*weightpileup*Lumi/LumiMC (default)
    - 1: isPassed selection weightinit*Lumi/LumiMC
    - 2: low mass selection weightlowmass
    - 3: high mass selection weighthighmass
   */
  double GetWeight(int weightType, ChrisLib::MapBranches mapBranches);

  /**\brief Returns true if the event is selected and false if not\n
     selectionType%10 gives the basic selection (linked to the weight):
    - 0: no selection
    - 1: isPassed (skimmed files)
    - X2: low mass
    - X3: high mass
    selectionType/10 gives the extra cut to be applied if needed:
    - 1X: no extra cut (so basic low/high mass selection)
    - 2X: + mgamgam cut (mh=125.09 +-4.7(4.3) GeV for low (high) mass)
    - 3X: + mhh cut (mhh < (>) 350 GeV for low (high) mass)
  */
  bool IsEventSelected(int selectionType, ChrisLib::MapBranches mapBranches);

  bool IsLowMass(ChrisLib::MapBranches mapBranches);
  bool IsHighMass(ChrisLib::MapBranches mapBranches);

  TH1D* InitialiseHist(std::string histName, std::string variable); 

  void DrawDistriForLambdas(TFile *inFile, std::vector<std::string> listVariables, std::vector<int> vectCategories, std::vector<std::string> vectSamples, std::string path, std::string extension="eps");

  std::vector<double> ReturnExtremalBins(TH1* hist);
  
  void DrawCompLONLOForLambdas(TFile *LOFile, TTree *LOTree, TFile *NLOFile, TTree *NLOTree, std::list<std::string> vectVariables, std::string savePath);

};
