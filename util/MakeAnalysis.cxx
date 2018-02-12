#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TROOT.h"

#include "HHAnalysis/HHAnalysis.h"

using namespace std;

int main(int argc, char *argv[])
{

  HHAnalysis hh;
  //  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/LOSignal.root"}, 0, {"mbbgamgam", "pTb1", "pTb2", "pTbb", "pTgam1", "pTgam2", "pTgamgam", "dRbb", "dRgamgam", "dRgamb"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/DistriLO/UnskimmedHighLevel.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/LOSignal.root"}, 0, {"mbbgamgam", "pTb1"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection0.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/skim_file/LOSignal.root"}, 1, {"mbbgamgam", "pTb1"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection1.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/skim_file/LOSignal.root"}, 12, {"mbbgamgam", "pTb1"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection12.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/skim_file/LOSignal.root"}, 13, {"mbbgamgam", "pTb1"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection13.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/skim_file/LOSignal.root"}, 22, {"mbbgamgam", "pTb1", "mgamgam"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection22.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/skim_file/LOSignal.root"}, 23, {"mbbgamgam", "pTb1", "mgamgam"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection23.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/skim_file/LOSignal.root"}, 32, {"mbbgamgam", "pTb1"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection32.root");
  hh.CreateSaveDistri({"/sps/atlas/d/delgove/private/HH/tree_v19/skim_file/LOSignal.root"}, 33, {"mbbgamgam", "pTb1"}, {2, 1} ,"/sps/atlas/a/aguerguichon/HHAnalysis/Test/TestSelection33.root");

  TFile *inFile=TFile::Open("/sps/atlas/a/aguerguichon/HHAnalysis/DistriLO/UnskimmedHighLevel.root");

  //hh.DrawDistriForLambdas(inFile, {"mbbgamgam", "pTb1", "pTb2", "pTbb", "pTgam1", "pTgam2", "pTgamgam", "dRbb", "dRgamgam", "dRgamb" }, {2, 1}, {"ALL"}, "/sps/atlas/a/aguerguichon/HHAnalysis/DistriLO/UnskimmedHighLevel_", "pdf");

  inFile->Close();
  cout<<"End of program.\n";
}
