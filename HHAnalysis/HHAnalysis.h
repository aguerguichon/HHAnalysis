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

  /* Reads the configuration file. Different supported options are:
    - inFile (vector<string>): root files to be used for the study.
    - variable (vector<string>): variables used to create the histograms and to be drawn.
    - sample (vector<string>): samples used. Keyword 'ALL' to used them all
    - category (vector<int>)b-tag categories used.
    - savePathPlot (string): path were the plots will be saved. Can be the absolute path or the absolute path+beginning of the name (name of the plot is automatically of the form tagcatX_varY).
    - outFileName (string): name of the root file where histograms will be stored.
    - selectionType (int): type of selection basic (selectionType%10) or extra (selectionType/10). More details in HHAnalysis::IsEventSelected.
    - extraInfo (string): extra information to be written on plots (default value "").
    - infoForWorkspace (int): information to be stored and used for workspace (default value 0). Possible values are 0 (none), 1 (mass), 2 (yields), 12 (yields +fit yields).
   */
  HHAnalysis(std::string configFileName);
  ~HHAnalysis();


  float GetExpectedSignificance(TH1 *histBkg, TH1 *histSignal, int mode, int method, double maxUnc, double *error);

  double getExpectedSignificance (TH1 *hBkg, TH1 *hSig, int mode,int method, double maxUnc,double *error);

  /**\brief Creates an histogram for each category of btag, each sample and each variable of interest. Saves the TH1D* of output in a root file.
   */
  //  void CreateSaveDistri(std::vector<std::string> vectInFiles, int selectionType, std::vector<std::string> vectVariables, std::vector<int> vectCategories, std::string outputFileName);
  void CreateSaveDistri(bool saveMassForWorkspace=0);

  /**\brief Returns the weight according to the weightType:
    - 0: for no selection weightMC*weightvertex*weightpileup*Lumi/LumiMC (default)
    - 1: isPassed selection weightinit*Lumi/LumiMC
    - 2: low mass selection weightlowmass
    - 3: high mass selection weighthighmass
   */
  double GetWeight(int weightType, ChrisLib::MapBranches mapBranches);

  /**\brief Returns true if the event is selected and false if not\n
     selectionType%10 gives the basic selection (linked to the weight):
    - X0: no selection
    - X1: isPassed (skimmed files)
    - X2: low mass
    - X3: high mass

    selectionType/10 gives the extra cut to be applied if needed:
    - X: no extra cut (so basic low/high mass selection)
    - 1X: + mgamgam cut (mh=125.09 +-4.7(4.3) GeV for low (high) mass)
    - 2X: + mhh cut (mhh < (>) 350 GeV for low (high) mass)
    - 3X: + mhh cut "reversed" (mhh > (<) 350 GeV for low (high) mass)
  */
  bool IsEventSelected(int selectionType, ChrisLib::MapBranches mapBranches);

  bool IsLowMass(ChrisLib::MapBranches mapBranches);
  bool IsHighMass(ChrisLib::MapBranches mapBranches);

  void InitialiseHist(TH1D* &hist, std::string histName, std::string variable); 

  void DrawDistriForLambdas(std::string extension="eps");

  void MakePdf( std::string latexFileName, std::vector<std::string> vectHistNames, std::string comment );

  void SaveYields(bool fitYield=0);
  void FitYields(TH1D *histYield, std::string name);
  std::vector<double> ReturnExtremalBins(TH1* hist);
  
  std::string GetOutFileName();
  int GetTypeInfoForWorkspace();

 private:
  std::map <std::string, TH1D*> m_mapHist;

  std::vector <std::string> m_vectInFiles;
  std::vector <std::string> m_vectVariables;
  std::vector <std::string> m_vectSamples;
  std::vector <int> m_vectCategories;

  std::string m_savePathPlot;
  std::string m_outFileName;
  std::string m_extraInfo;

  int m_selectionType;
  int m_infoForWorkspace;
};
