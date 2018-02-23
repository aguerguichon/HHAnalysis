#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include "boost/program_options.hpp"

#include "TFile.h"                                                                      
#include "TF1.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TKey.h"
#include "TString.h"

#include "HHAnalysis/HHAnalysis.h"
#include "PlotFunctions/DrawPlot.h"
#include "PlotFunctions/DrawOptions.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "PlotFunctions/SideFunctions.h"
#include "PlotFunctions/MapBranches.h"

using namespace ChrisLib;
using namespace std;

namespace po = boost::program_options;


//==========================================================================================
HHAnalysis::HHAnalysis(string configFileName){
  po::options_description configOptions("Configuration options");
  configOptions.add_options()
    ("help", "Help")
    ("inFile", po::value<vector<string>>(&m_vectInFiles), "Root files to be used for the study.")
    ("variable", po::value<vector<string>>(&m_vectVariables), "Variables used to create the histograms and to be drawn.")
    ("sample", po::value<vector<string>>(&m_vectSamples), "Samples used. Keyword 'ALL' to used them all. ")
    ("category", po::value<vector<int>>(&m_vectCategories), "b-tag categories used.")
    ("savePathPlot", po::value<string>(&m_savePathPlot), "Path were the plots will be saved. Can be the absolute path or the absolute path+beginning of the name (name is of the form tagcatX_varY).")
    ("outFileName", po::value<string>(&m_outFileName), "Name of the root file where histograms will be stored.")
    ("selectionType", po::value<int>(&m_selectionType), "Type of selection basic (selectionType%10) or extra (selectionType/10)")
    ("extraInfo", po::value<string>(&m_extraInfo)->default_value(""), "Extra information to be written on plots.")
    ("infoForWorkspace", po::value<int>(&m_infoForWorkspace)->default_value(0), "Information to be stored and used for workspace.")
    ;
  
  po::variables_map vm;
  ifstream ifs( configFileName, ifstream::in );
  po::store(po::parse_config_file(ifs, configOptions), vm);
  po::notify(vm);

  cout<<"Configuration file "<<configFileName<<" loaded."<<endl;
}


HHAnalysis::~HHAnalysis(){
  m_mapHist.clear();
  m_vectInFiles.clear();
  m_vectVariables.clear();
  m_vectSamples.clear();
  m_vectCategories.clear();
}


//==========================================================================================
float HHAnalysis::GetExpectedSignificance(TH1 *histBkg, TH1 *histSignal, int mode, int method, double maxUnc, double *error) {
  float significance{0}; 
  double iS{0};
  double iB{0};
  double errorS{0};
  double errorB{0};
  double dS{0};
  double dB{0};
  double deltaSigma{0};

  if (!histBkg || !histSignal) throw invalid_argument( "HHAnalysis::GetExpectedSignificance: no histograms provided." );

  // if (method >=0 && mode==2){
  //   vector<int > *newbins = new vector<int>;
  //   *newbins = getRebinBins(histBkg,histSignal,method,maxUnc);
  //   rebinHisto(histSignal,newbins, false,false);   
  //   rebinHisto(histBkg,newbins, false,false);   
  // }
  if (histBkg->GetNbinsX() != histSignal->GetNbinsX()) throw invalid_argument( "HHAnalysis::GetExpectedSignificance: signal and bkg histograms do not have the same number of bins.");

  if (mode==1){
    if (histBkg->Integral(0, histBkg->GetNbinsX()+2)!=0){
      iS = histSignal->IntegralAndError(0, histSignal->GetNbinsX()+2, errorS);
      iB = histBkg->IntegralAndError(0, histBkg->GetNbinsX()+2, errorB);
      //      sig = pow(histSignal->Integral(0,histSignal->GetNbinsX()+2) /sqrt(histBkg->Integral(0,histBkg->GetNbinsX()+2)),2) ;  
      significance =  2 * ( (iS+iB) * log(1 + iS/iB) - iS ) ;
      dS = log(1+iS/iB);
      dB = log(1+iS/iB)-iS/iB;
      deltaSigma = sqrt( dS*dS*errorS*errorS + dB*dB*errorB*errorB );
      *error = deltaSigma/sqrt(significance);
    }
    else {
      significance = 0;
      *error = 0;
    }
  }

  if (mode==2){ 
    *error = 0;     
    for (int iBin=0; iBin<= histSignal->GetNbinsX(); iBin++) {
      errorS=histSignal->GetBinError(iBin);
      errorB=histBkg->GetBinError(iBin);
      iS = histSignal->GetBinContent(iBin);
      iB = histBkg->GetBinContent(iBin);
      if ( iB > 0 ) {
	significance += 2 * ( (iS+iB) * log(1 + iS/iB) - iS );
      	dS = log(1+iS/iB);
	dB = log(1+iS/iB)-iS/iB;
	deltaSigma = sqrt( dS*dS*errorS*errorS + dB*dB*errorB*errorB);
	*error += deltaSigma*deltaSigma;
      }
    }
    //cout << "start significance fonction " << endl;
    *error = sqrt(*error/significance);
  }
   
  return sqrt(significance);

}


//==================================
//HHAnalysis::getExpectedSignificance(TH1 *hBkg, TH1 *hSig, int mode, int method, double maxUnc,double *error) {  
//From David's code: /afs/in2p3.fr/home/d/delgove/public/HH/Analysis/Tools/ToolsForOptimisation.C
double HHAnalysis::getExpectedSignificance(TH1 *hBkg, TH1 *hSig, int mode, int method, double maxUnc,double *error) {  
  double sig = 0; 
  if (!hBkg || !hSig) {
    cout << "no histogram " << endl;
    return sig;
  }  
  // if (method>=0 && mode==2){
  //   vector<int > *newbins = new vector<int>;
  //   *newbins = getRebinBins(hBkg,hSig,method,maxUnc);
  //   rebinHisto(hSig,newbins, false,false);   
  //   rebinHisto(hBkg,newbins, false,false);   
  // }
  if (hBkg->GetNbinsX() != hSig->GetNbinsX()) {
    cout << "not the same number of bin " << endl;
    return sig;
  }  
  if (mode==1){
    if (hBkg->Integral(0,hBkg->GetNbinsX()+2)!=0){
      double error_s;
      double error_b;
      double iS = hSig->IntegralAndError(0,hSig->GetNbinsX()+2,error_s);
      double iB = hBkg->IntegralAndError(0,hBkg->GetNbinsX()+2,error_b);
      //      sig = pow(hSig->Integral(0,hSig->GetNbinsX()+2) /sqrt(hBkg->Integral(0,hBkg->GetNbinsX()+2)),2) ;  
      sig =  2 * ( (iS+iB) * log( 1 + iS / iB ) - iS ) ;
      double ds = log(1+iS/iB);
      double db = log(1+iS/iB)-iS/iB;
      double deltasigma = sqrt(ds*ds*error_s*error_s+db*db*error_b*error_b);
      *error = deltasigma/sqrt(sig);
      /*cout << "start significance fonction " << endl;
      cout << " Signal : " << iS << "+/- " << error_s << endl;                                                                                                          
      cout << "Background : " << iB << "+/- " << error_b << endl;     
      cout << sqrt(2 * ( (iS+iB) * log( 1 + iS / iB ) - iS ) )<< endl;
      cout << "start significance fonction " << endl;*/
    }
    else {
      sig = 0;
      *error = 0;
    }
  }
  if (mode==2){ 
    *error = 0;     
    //    cout << "start significance fonction " << endl;
    int nbins = hSig->GetNbinsX();
    for (int bin=0; bin<=nbins; bin++) {
      double error_s=hSig->GetBinError(bin);
      double error_b=hBkg->GetBinError(bin);
      double iS = hSig->GetBinContent(bin);
      double iB = hBkg->GetBinContent(bin);
      /*cout << " Signal : " << iS << "+/- " << error_s << endl;
	cout << "Background : " << iB << "+/- " << error_b << endl;*/
      if ( iB > 0 ) {
	sig += 2 * ( (iS+iB) * log( 1 + iS / iB ) - iS );
      	double ds = log(1+iS/iB);
	double db = log(1+iS/iB)-iS/iB;
	double deltasigma = sqrt(ds*ds*error_s*error_s+db*db*error_b*error_b);
	*error += deltasigma*deltasigma;
	//	cout << sqrt(2 * ( (iS+iB) * log( 1 + iS / iB ) - iS ) )<< endl;
      }
    }
    //cout << "start significance fonction " << endl;
    *error = sqrt(*error/sig);
  }
   
  return sqrt(sig);
}

//==================================================================================
void HHAnalysis::CreateSaveDistri(bool saveMassForWorkspace){
  cout<<"HHAnalysis::CreateSaveDistri\n";

  string *sampleName=new string;
  string histName;
  TBranch *sampleBranch; //!
  double weight{0};
  TFile *inFile{0};
  TTree *inTree{0};
  list<string> listVariables(m_vectVariables.begin(),m_vectVariables.end());

  bool doSave=0; vector <ofstream*> vectOfstream(m_vectCategories.size()); 

  listVariables.merge({"tagcat"});
  listVariables.merge({"weightMC", "weightvertex", "weightpileup", "weightinit","Lumi", "LumiMC", "weightlowmass", "weighthighmass"});
  listVariables.merge({"pTb1corlow", "pTb2corlow", "pTb1corhigh", "pTb2corhigh", "mbbcorlow", "mbbcorhigh"});

  TH1::AddDirectory(kFALSE); 

  for (unsigned int iFile=0; iFile<m_vectInFiles.size(); iFile++){
    inFile=TFile::Open( (m_vectInFiles[iFile]).c_str() );
    if (!inFile) throw invalid_argument( "HHAnalysis::CreateSaveDistri: no file provided." );
    inTree=(TTree*)inFile->Get("ntuple");
    if (!inFile) throw invalid_argument( "HHAnalysis::CreateSaveDistri: no tree provided." );

    if (saveMassForWorkspace==1) {
      for (unsigned int iCat=0; iCat<m_vectCategories.size(); iCat++){
	vectOfstream[iCat] = new ofstream( (m_savePathPlot+"tagcat"+to_string(m_vectCategories[iCat])+"_listOfEvents.csv").c_str(), ios::out);
      }
      doSave=1;
    }

    MapBranches MB; 
    MB.LinkTreeBranches(inTree, 0, listVariables);
    inTree->SetBranchStatus("sample", 1);
    inTree->SetBranchAddress("sample", &sampleName, &sampleBranch);
  
    for (unsigned int iEntry=0; iEntry<inTree->GetEntries(); iEntry++){
      inTree->GetEntry(iEntry);
      if ( !IsEventSelected(m_selectionType, MB) || !IsSampleSelected(sampleName) ) continue;
      weight=GetWeight(m_selectionType%10, MB);
      if ( m_vectInFiles[iFile].find("LO")!=string::npos ) weight*=0.34*2.28;
	
      //Create entries of m_mapHist and fill hists
      for (unsigned int iCat=0; iCat<m_vectCategories.size(); iCat++){
	if ( MB.GetInt("tagcat")!=m_vectCategories[iCat] ) continue;
	if (doSave && m_selectionType%10==2) *vectOfstream[iCat]<<MB.GetDouble("mgamgam")<<"\t"<<MB.GetDouble("mbbcorlow")<<"\t"<<MB.GetDouble("mbbgamgamcorlow")<<"\t"<<MB.GetDouble("mbbgamgamcorlowmodified")<<"\t"<<"1.0"<<"\n";
	if (doSave && m_selectionType%10==3) *vectOfstream[iCat]<<MB.GetDouble("mgamgam")<<"\t"<<MB.GetDouble("mbbcorhigh")<<"\t"<<MB.GetDouble("mbbgamgamcorhigh")<<"\t"<<MB.GetDouble("mbbgamgamcorhighmodified")<<"\t"<<"1.0"<<"\n";

	for (unsigned int iVar=0; iVar<m_vectVariables.size(); iVar++ ){
	  histName="tagcat"+to_string( MB.GetInt("tagcat") )+"_var"+m_vectVariables[iVar]+"_sample"+*sampleName;

	  if ( !m_mapHist.count(histName) ) { InitialiseHist(m_mapHist[histName], histName, m_vectVariables[iVar] );} 
	  if ( MB.IsInt(m_vectVariables[iVar].c_str()) ) m_mapHist[histName]->Fill( MB.GetInt( m_vectVariables[iVar].c_str() ), weight );
	  else if ( m_vectVariables[iVar].find("phi")!=string::npos ) m_mapHist[histName]->Fill( MB.GetDouble( m_vectVariables[iVar].c_str() )*180/acos(-1), weight );
	  else m_mapHist[histName]->Fill( MB.GetDouble( m_vectVariables[iVar].c_str() ), weight );
	} // end iVar
      }//end iCat  
    } //end iEntry
    inFile->Close("R");
  }//end iFile

  //Saving histograms in a root file
  TFile *outputFile=new TFile (m_outFileName.c_str(), "RECREATE");
  if (!outputFile) throw invalid_argument("HHAnalysis::CreatSaveDistri outputFile does not exist.");
  for(auto it : m_mapHist) it.second->Write("", TObject::kWriteDelete);
  outputFile->Close("R");
  cout<<"Histograms saved in "<<m_outFileName<<endl;
  
  if (doSave) {for (unsigned int iStream=0; iStream<vectOfstream.size(); iStream++){vectOfstream[iStream]->close();} }
  vectOfstream.clear();
  delete outputFile; outputFile=0;  
  delete sampleName; sampleName=0;
  delete inFile; inFile=0; 
  cout<<"HHAnalysis::CreateSaveDistri done.\n";
}

//=================================================================================
bool HHAnalysis::IsSampleSelected(string *sampleName){
  for(unsigned int iSample=0; iSample<m_vectSamples.size(); iSample++){
    if (m_vectSamples[iSample].compare("ALL")==0) return true; //take all samples
    else if (m_vectSamples[iSample].compare(*sampleName)==0) return true;
  }
  return false;
}

//==================================================================================
bool HHAnalysis::IsEventSelected(int selectionType, MapBranches MB){

  if (selectionType%10 == 0 || selectionType%10 == 1){ //unskimmed or isPassed (ie skimmed)
    switch (selectionType/10){
    case 0: {return true; break; } //no extra selection
    case 2: { if (MB.GetDouble("mbbgamgam")<350  ) return true; else return false; break;}
    case 3: { if (MB.GetDouble("mbbgamgam")>350  ) return true; else return false; break;} 
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for unskimmed of skimmed selection: "+selectionType);
    }
  }

  else if (selectionType%10 == 2){ //lowmass
    if ( !IsLowMass(MB) ) return false;
    switch(selectionType/10){
    case 0: { return true; break;} //basic low mass selection
    case 1: { if (MB.GetDouble("mgamgam")>120.39 && MB.GetDouble("mgamgam")<129.79 ) return true; else return false; break;} //adding mgamgam cut
    case 2: { if (MB.GetDouble("mbbgamgam")<350  ) return true; else return false; break;} //adding mhh cut
    case 3: { if (MB.GetDouble("mbbgamgam")>350  ) return true; else return false; break;} //reversed mhh cut
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for low mass selection: "+selectionType);
    }
  }

  else if (selectionType%10 == 3){ //highmass
    if ( !IsHighMass(MB) ) return false;
    switch(selectionType/10){
    case 0: { return true; break;} //basic high mass selection
    case 1: { if (MB.GetDouble("mgamgam")>120.79 && MB.GetDouble("mgamgam")<129.39  ) return true; else return false; break;} //adding mgamgam cut
    case 2: { if (MB.GetDouble("mbbgamgam")>350  ) return true; else return false; break;}
    case 3: { if (MB.GetDouble("mbbgamgam")<350  ) return true; else return false; break;} //reversed mhh cut
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for high mass selection: "+selectionType);
    }
  }
  return false;
}

//==================================================================================
double HHAnalysis::GetWeight(int weightType, MapBranches MB){
  double weight{0};
  switch (weightType){
  case 0: {weight=MB.GetDouble("weightMC")*MB.GetDouble("weightvertex")*MB.GetDouble("weightpileup")*(MB.GetDouble("Lumi")/MB.GetDouble("LumiMC") ); break;} //for no selection cut at all
  case 1:{ weight=MB.GetDouble("weightinit")*(MB.GetDouble("Lumi")/MB.GetDouble("LumiMC") ); break;} //for selection cut at isPassed (equivalent to skimmed files)
  case 2:{weight=MB.GetDouble("weightlowmass"); break;}
  case 3:{weight=MB.GetDouble("weighthighmass"); break;}

  default: weight=MB.GetDouble("weightMC")*MB.GetDouble("weightvertex")*MB.GetDouble("weightpileup")*(MB.GetDouble("Lumi")/MB.GetDouble("LumiMC") ); //no selection at all
  }
  
  return weight;
}


//=================================================================================
bool HHAnalysis::IsLowMass(MapBranches MB){
  if ( MB.GetDouble("pTb1corlow")>40 && MB.GetDouble("pTb2corlow")>25 && MB.GetDouble("mbbcorlow")>80 && MB.GetDouble("mbbcorlow")<140 ) return true;
  else return false;
}


//=================================================================================
bool HHAnalysis::IsHighMass(MapBranches MB){
  if ( MB.GetDouble("pTb1corhigh")>100 && MB.GetDouble("pTb2corhigh")>30 && MB.GetDouble("mbbcorhigh")>90 && MB.GetDouble("mbbcorhigh")<140 ) return true;
  else return false;
}


//==================================================================================
void HHAnalysis::InitialiseHist(TH1D* &hist, string histName, string strVariable){
  TString var=strVariable;
  if (var.Contains("mv2c")) hist= new TH1D (histName.c_str(), "", 200, -1, 1);
  else if (var.BeginsWith("mgamgam")) hist= new TH1D (histName.c_str(), "", 60, 110, 140);
  else if (var=="mbb") hist= new TH1D (histName.c_str(), "", 40, 0, 200);
  else if (var.BeginsWith("mhhtruth")) hist= new TH1D (histName.c_str(), "", 48, 100e3, 1500e3);
  else if (var.BeginsWith("met")) hist= new TH1D (histName.c_str(), "", 40, 0, 200);
  else if (var.BeginsWith("m")) hist= new TH1D (histName.c_str(), "", 100, 0, 1000);
  else if (var.Contains("pT")) hist= new TH1D (histName.c_str(), "", 60, 0, 600);
  else if (var.BeginsWith("d")) hist= new TH1D (histName.c_str(), "", 25, 0, 5);
  else if (var.Contains("eta")) hist= new TH1D (histName.c_str(), "", 50, -5, 5);
  else if (var.Contains("phi")) hist= new TH1D (histName.c_str(), "", 18, -180, 180);
  else if (var.BeginsWith("n") && var.Contains("jet")) hist= new TH1D (histName.c_str(), "", 10, 0, 10);
  else hist= new TH1D (histName.c_str(), "", 2000, -1000, 1000);
  return;
}


//==================================================================================
void HHAnalysis::DrawDistriForLambdas(string extension){
  cout<<"HHAnalysis::DrawDistriForLambdas done.\n";

  TObjArray *objArrayString{0};
  string histName, plotName, catName, sampleName, varName, legLatex;
  vector <string> vectOpt, vectHistNames;
  vector <double> vectExtremalBins;
  vector <TH1*> vectHistTmp;
  TString tmp, name;

  for (unsigned int iCat=0; iCat<m_vectCategories.size(); iCat++){
    for (unsigned int iVar=0; iVar<m_vectVariables.size(); iVar++){
      vectExtremalBins.push_back(15e10); //arbitrary value
      vectExtremalBins.push_back(0);
      plotName="tagcat"+to_string(m_vectCategories[iCat])+"_var"+m_vectVariables[iVar];
      for(auto &it : m_mapHist){ //it.first =histName, it.second=hist
	name=it.first;
	if ( !name.BeginsWith((plotName+"_").c_str()) ) continue;
	objArrayString = name.Tokenize("_");
	for(int iString=0; iString<objArrayString->GetEntriesFast(); iString++){
	  tmp = ((TObjString*)objArrayString->At(iString))->GetString();
	  if ( tmp.Contains("tagcat") ){ catName=tmp.ReplaceAll("tagcat","");}
	  if ( tmp.Contains("var") ) varName=tmp.ReplaceAll("var","");
	  if ( tmp.Contains("sample") ) sampleName=tmp.ReplaceAll("sample","");
	}
	vectOpt.push_back( ("legend="+ sampleName ).c_str() );
	
	if ( ReturnExtremalBins(it.second)[0]<vectExtremalBins[0] ) vectExtremalBins[0]=ReturnExtremalBins(it.second)[0];
	if ( ReturnExtremalBins(it.second)[1]>vectExtremalBins[1] ) vectExtremalBins[1]=ReturnExtremalBins(it.second)[1];
	vectHistTmp.push_back(it.second);
      } // end it m_mapHist ie loop over hist

      DrawOptions drawOpt;
      drawOpt.AddOption(vectOpt);
      drawOpt.AddOption("rangeUserX", (to_string(vectExtremalBins[0])+" "+to_string(vectExtremalBins[1])).c_str() );
      drawOpt.AddOption("outName", (m_savePathPlot+plotName).c_str() );
      drawOpt.AddOption("legendPos", "0.8 0.9");
      drawOpt.AddOption("latex", ("Category "+ catName +" b-tag").c_str() );
      drawOpt.AddOption("latexOpt", "0.45 0.9");
      if (m_extraInfo!="") {drawOpt.AddOption("latex", m_extraInfo.c_str()); drawOpt.AddOption("latexOpt", "0.45 0.85");}
      drawOpt.AddOption("xTitle", varName);
      drawOpt.AddOption("normalize", "1");
      drawOpt.AddOption("extension", extension);
      drawOpt.Draw( vectHistTmp );
      
      vectHistNames.push_back(m_savePathPlot+plotName);
      vectOpt.clear();
      vectExtremalBins.clear();
      vectHistTmp.clear();
    } //end iVar
  }//end iCat
  
  plotName=m_savePathPlot+"Summary";
  MakePdf( plotName, vectHistNames, m_extraInfo );
  system( ("mv *Summary.pdf "+plotName+".pdf").c_str() );
  system( "rm *Summary*"); 
  
  delete objArrayString; objArrayString=0;
  cout<<"HHAnalysis::DrawDistriForLambdas done.\n";
}




//==================================================================================
void HHAnalysis::MakePdf( string latexFileName, vector<string> vectHistNames, string comment ){
  fstream stream;
  string var;
  latexFileName+=".tex";
  stream.open( latexFileName.c_str(), fstream::out | fstream::trunc );
  WriteLatexHeader( stream, "HH study" , "Antinea Guerguichon" );

  stream <<comment <<"\\newline"<<endl;
  stream << "\\indent Variables: ";
  for (unsigned int iVar=0; iVar< m_vectVariables.size(); iVar++)    {
    if (m_vectVariables[iVar].find("_")!=string::npos) var=ReplaceString("_", "\\_")(m_vectVariables[iVar]);
    else var=m_vectVariables[iVar];

    if (iVar == m_vectVariables.size()-1) stream<< var <<"\\newline  ";
    else  stream  << var <<", ";
  }

  WriteLatexMinipage( stream, vectHistNames, 2);
  stream << "\\end{document}" << endl;
  string commandLine = "pdflatex  -interaction=batchmode " + latexFileName;
  system( commandLine.c_str() );
  system( commandLine.c_str() );
  system( commandLine.c_str() );
  stream.close();

  cout<<"HHAnalysis::MakePdf "+latexFileName+"  Done"<<endl;
  return;
}

//==================================================================================
void HHAnalysis::SaveYields(bool fitYield){
  cout<<"HHAnalysis::SaveYields\n";
  TString histName, tmp;
  TFile *outRootFile{0};
  if (fitYield) outRootFile=new TFile((m_savePathPlot+"fit.root").c_str(), "RECREATE");  
  TH1D *histYield{0};
  double error{0};

  for (unsigned int iCat=0; iCat<m_vectCategories.size(); iCat++){
    histYield=new TH1D(("yields_tagcat"+to_string(m_vectCategories[iCat])).c_str(), "", 24, -12.5, 11.5);
    ofstream outFile( (m_savePathPlot+"tagcat"+to_string(m_vectCategories[iCat])+"_yields.csv").c_str(), ios::out );
    for (auto &it : m_mapHist){
      histName=it.first;
      if (!histName.Contains(("tagcat"+to_string(m_vectCategories[iCat])+"_varmbbgamgam_").c_str() ) ) continue;
      TObjArray *objArrayString = histName.Tokenize("_");
      for(int iString=0; iString<objArrayString->GetEntriesFast(); iString++){
	tmp = ((TObjString*)objArrayString->At(iString))->GetString();
	if ( tmp.Contains("sample") ) {tmp=tmp.ReplaceAll("sample",""); tmp=tmp.ReplaceAll("hh",""); tmp=tmp.ReplaceAll("minus", "-"); tmp=tmp.ReplaceAll("plus", ""); }
      }
      outFile<<tmp<<"\t"<<it.second->Integral(0, it.second->GetNbinsX()+1)<<endl;
      if (fitYield) {
	histYield->SetBinContent( histYield->GetXaxis()->FindFixBin(tmp.Atof()), it.second->Integral(0, it.second->GetNbinsX()+1) );
	it.second->IntegralAndError(0, it.second->GetNbinsX()+1, error);
	histYield->SetBinError( histYield->GetXaxis()->FindFixBin(tmp.Atof()), error );
      }
    }
    outFile.close();
    if (fitYield) FitYields( histYield, to_string(m_vectCategories[iCat]));
  }

  if (fitYield) outRootFile->Close();
  delete outRootFile; outRootFile=0;
  cout<<"HHAnalysis::SaveYields done.\n";
  return;
}
//=================================================================================
void HHAnalysis::FitYields(TH1D* histYield, string name){
  cout<<"HHAnalysis::FitYields\n";
  ofstream outFile( m_savePathPlot+"tagcat"+name+"_lambdaParabola.csv", ios::out);
  string latex;
  TF1 *parabola=new TF1(("parabola"+name).c_str(), "[0]+[1]*x+[2]*x*x", -10, 10);
  histYield->SetMarkerStyle(20);
  histYield->SetMarkerColor(kBlack);
  histYield->Fit( parabola, "R", "P");
  histYield->Write();
  DrawOptions drawOpt;
  drawOpt.AddOption("outName", (m_savePathPlot+"tagcat"+name+"_yieldFit").c_str() );
  drawOpt.AddOption("latex", ("Category "+ name +" b-tag").c_str() );
  drawOpt.AddOption("latexOpt", "0.45 0.9");
  if (m_extraInfo!="") {drawOpt.AddOption("latex", m_extraInfo.c_str()); drawOpt.AddOption("latexOpt", "0.45 0.85");}
  drawOpt.AddOption("latex", "Fitting function p0+p1*x+p2*x*x") ;
  drawOpt.AddOption("latexOpt", "0.6 0.8");
  latex="p0: "+to_string(histYield->GetFunction(("parabola"+name).c_str())->GetParameter(0));
  drawOpt.AddOption("latex", latex.c_str()) ;
  drawOpt.AddOption("latexOpt", "0.6 0.77");
  latex="p1: "+to_string(histYield->GetFunction(("parabola"+name).c_str())->GetParameter(1));
  drawOpt.AddOption("latex", latex.c_str()) ;
  drawOpt.AddOption("latexOpt", "0.6 0.74");
  latex="p2: "+to_string(histYield->GetFunction(("parabola"+name).c_str())->GetParameter(2));
  drawOpt.AddOption("latex", latex.c_str()) ;
  drawOpt.AddOption("latexOpt", "0.6 0.71");
  drawOpt.AddOption("xTitle", "__HASHTAGlambda");
  drawOpt.AddOption("drawStyle", "3");
  drawOpt.AddOption("yTitle", "Number of expected events");
  drawOpt.AddOption("rangeUserY", "0 0.99");
  drawOpt.Draw( {histYield} );
  
  outFile<<"Fittingfunction"<<"\t"<<"p0+p1*x+p2*x*x\n"<<"p0\t"<<histYield->GetFunction(("parabola"+name).c_str())->GetParameter(0)<<"\n"<<"p1\t"<<histYield->GetFunction(("parabola"+name).c_str())->GetParameter(1)<<"\n"<<"p2\t"<<histYield->GetFunction(("parabola"+name).c_str())->GetParameter(2)<<"\n";
  outFile.close();
  cout<<"Fit of yields saved in "<<m_savePathPlot<<"tagcat"<<name<<"_yieldFit.pdf and parameters in "<<m_savePathPlot<<"tagcat"<<name<<"_lambdaParabola.csv\n";
  cout<<"HHAnalysis::FitYields done.\n";
  return;
}

//==================================================================================
vector<double> HHAnalysis::ReturnExtremalBins(TH1* hist){
  unsigned int lowBin = 1;
  unsigned int upBin = hist->GetNbinsX();
  double minX=hist->GetXaxis()->GetXmin();
  double maxX=hist->GetXaxis()->GetXmax();
  vector <double> vectExtremalBins;

  while ( hist->GetBinContent( lowBin ) == 0 && lowBin!=upBin ) lowBin++;
  while ( hist->GetBinContent( upBin ) ==0 && lowBin!=upBin ) upBin--;
  if ( lowBin != upBin ) {
    minX = hist->GetXaxis()->GetBinLowEdge( lowBin );
    maxX = hist->GetXaxis()->GetBinUpEdge( upBin );
  }
  vectExtremalBins.push_back(minX);   vectExtremalBins.push_back(maxX);
  return vectExtremalBins;
}

//==================================================================================
string HHAnalysis::GetOutFileName(){
  return m_outFileName;
}

//==================================================================================
int HHAnalysis::GetTypeInfoForWorkspace(){
  return m_infoForWorkspace;
}
