#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TError.h"

#include "HHAnalysis/HHAnalysis.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description desc("Configuration files.");
  vector<string> configFiles;

  //define all options in the program 
  desc.add_options()
    ( "help", "Display this help message")
    ( "configFiles", po::value<vector <string> >(&configFiles), "" )
    ;

  //Define options gathered by position 
  po::positional_options_description p;
  p.add("configFiles", -1);

  // create a map vm that contains options and all arguments of options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {cout << desc; return 0;}
  //=========================================================

  gErrorIgnoreLevel =kError; //removes warnings from root
  int infoForWorkspace{0};

  for (unsigned int iFile=0; iFile <configFiles.size(); iFile++){
    HHAnalysis hh ( configFiles[iFile] );
    infoForWorkspace = hh.GetTypeInfoForWorkspace();
    hh.CreateSaveDistri(infoForWorkspace);
    if (infoForWorkspace%10 == 2) hh.SaveYields(infoForWorkspace/10);
    hh.DrawDistriForLambdas("pdf");
  }

  cout<<"End of program.\n";
}
