#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include "boost/program_options.hpp"

#include "TFile.h"                                                                          
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
void HHAnalysis::CreateSaveDistri(){
  cout<<"HHAnalysis::CreateSaveDistri\n";

  string *sampleName=new string;
  string histName;
  TBranch *sampleBranch; //!
  double weight{0};
  TFile *inFile{0};
  TTree *inTree{0};
  list<string> listVariables(m_vectVariables.begin(),m_vectVariables.end());
  
  listVariables.merge({"tagcat"});
  listVariables.merge({"weightMC", "weightvertex", "weightpileup", "weightinit","Lumi", "LumiMC", "weightlowmass", "weighthighmass"});
  listVariables.merge({"pTb1corlow", "pTb2corlow", "pTb1corhigh", "pTb2corhigh", "mbbcorlow", "mbbcorhigh"});

  TH1::AddDirectory(kFALSE); 

  for (unsigned int iFile=0; iFile<m_vectInFiles.size(); iFile++){
    inFile=TFile::Open( (m_vectInFiles[iFile]).c_str() );
    if (!inFile) throw invalid_argument( "HHAnalysis::CreateSaveDistri: no file provided." );
    inTree=(TTree*)inFile->Get("ntuple");
    if (!inFile) throw invalid_argument( "HHAnalysis::CreateSaveDistri: no tree provided." );

    MapBranches mapBranches; 
    mapBranches.LinkTreeBranches(inTree, 0, listVariables);
    inTree->SetBranchStatus("sample", 1);
    inTree->SetBranchAddress("sample", &sampleName, &sampleBranch);
  
    for (unsigned int iEntry=0; iEntry<inTree->GetEntries(); iEntry++){
      inTree->GetEntry(iEntry);
      if ( !IsEventSelected(m_selectionType, mapBranches) ) continue;
      weight=GetWeight(m_selectionType%10, mapBranches);
      if ( m_vectInFiles[iFile].find("LO")!=string::npos ) weight*=0.34*2.28;
	
      //Create entries of m_mapHist and fill hists
      for (unsigned int iCat=0; iCat<m_vectCategories.size(); iCat++){
	if ( mapBranches.GetInt("tagcat")!=m_vectCategories[iCat] ) continue;
	for (unsigned int iVar=0; iVar<m_vectVariables.size(); iVar++ ){
	  histName="tagcat"+to_string( mapBranches.GetInt("tagcat") )+"_var"+m_vectVariables[iVar]+"_sample"+*sampleName;

	  if ( !m_mapHist.count(histName) ) { InitialiseHist(m_mapHist[histName], histName, m_vectVariables[iVar] );} 
	  if ( mapBranches.IsInt(m_vectVariables[iVar].c_str()) ) m_mapHist[histName]->Fill( mapBranches.GetInt( m_vectVariables[iVar].c_str() ), weight );
	  else if ( m_vectVariables[iVar].find("phi")!=string::npos ) m_mapHist[histName]->Fill( mapBranches.GetDouble( m_vectVariables[iVar].c_str() )*180/acos(-1), weight );
	  else m_mapHist[histName]->Fill( mapBranches.GetDouble( m_vectVariables[iVar].c_str() ), weight );
	} // end iVar
      }//end iCat  
    } //end iEntry
    inFile->Close("R");
  }//end iFile

  //Saving histograms in a root file to be post treated.
  TFile *outputFile=new TFile (m_outFileName.c_str(), "RECREATE");
  if (!outputFile) throw invalid_argument("HHAnalysis::CreatSaveDistri outputFile does not exist.");
  for(auto it : m_mapHist) it.second->Write("", TObject::kWriteDelete);
  outputFile->Close("R");
  cout<<"Histograms saved in "<<m_outFileName<<endl;
  
  delete outputFile; outputFile=0;  
  delete sampleName; sampleName=0;
  delete inFile; inFile=0; 
  cout<<"HHAnalysis::CreateSaveDistri done.\n";
}


//==================================================================================
bool HHAnalysis::IsEventSelected(int selectionType, MapBranches mapBranches){

  if (selectionType%10 == 0 || selectionType%10 == 1){ //unskimmed or isPassed (ie skimmed)
    switch (selectionType/10){
    case 0: {return true; break; } //no extra selection
    case 2: { if (mapBranches.GetDouble("mbbgamgam")<350  ) return true; else return false; break;}
    case 3: { if (mapBranches.GetDouble("mbbgamgam")>350  ) return true; else return false; break;} 
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for unskimmed of skimmed selection: "+selectionType);
    }
  }

  else if (selectionType%10 == 2){ //lowmass
    if ( !IsLowMass(mapBranches) ) return false;
    switch(selectionType/10){
    case 0: { return true; break;} //basic low mass selection
    case 1: { if (mapBranches.GetDouble("mgamgam")>120.39 && mapBranches.GetDouble("mgamgam")<129.79 ) return true; else return false; break;} //adding mgamgam cut
    case 2: { if (mapBranches.GetDouble("mbbgamgam")<350  ) return true; else return false; break;} //adding mhh cut
    case 3: { if (mapBranches.GetDouble("mbbgamgam")>350  ) return true; else return false; break;} //reversed mhh cut
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for low mass selection: "+selectionType);
    }
  }

  else if (selectionType%10 == 3){ //highmass
    if ( !IsHighMass(mapBranches) ) return false;
    switch(selectionType/10){
    case 0: { return true; break;} //basic high mass selection
    case 1: { if (mapBranches.GetDouble("mgamgam")>120.79 && mapBranches.GetDouble("mgamgam")<129.39  ) return true; else return false; break;} //adding mgamgam cut
    case 2: { if (mapBranches.GetDouble("mbbgamgam")>350  ) return true; else return false; break;}
    case 3: { if (mapBranches.GetDouble("mbbgamgam")<350  ) return true; else return false; break;} //reversed mhh cut
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for high mass selection: "+selectionType);
    }
  }
  return false;
}

//==================================================================================
double HHAnalysis::GetWeight(int weightType, MapBranches mapBranches){
  double weight{0};
  switch (weightType){
  case 0: {weight=mapBranches.GetDouble("weightMC")*mapBranches.GetDouble("weightvertex")*mapBranches.GetDouble("weightpileup")*(mapBranches.GetDouble("Lumi")/mapBranches.GetDouble("LumiMC") ); break;} //for no selection cut at all
  case 1:{ weight=mapBranches.GetDouble("weightinit")*(mapBranches.GetDouble("Lumi")/mapBranches.GetDouble("LumiMC") ); break;} //for selection cut at isPassed (equivalent to skimmed files)
  case 2:{weight=mapBranches.GetDouble("weightlowmass"); break;}
  case 3:{weight=mapBranches.GetDouble("weighthighmass"); break;}

  default: weight=mapBranches.GetDouble("weightMC")*mapBranches.GetDouble("weightvertex")*mapBranches.GetDouble("weightpileup")*(mapBranches.GetDouble("Lumi")/mapBranches.GetDouble("LumiMC") ); //no selection at all
  }
  
  return weight;
}


//=================================================================================
bool HHAnalysis::IsLowMass(MapBranches mapBranches){
  if ( mapBranches.GetDouble("pTb1corlow")>40 && mapBranches.GetDouble("pTb2corlow")>25 && mapBranches.GetDouble("mbbcorlow")>80 && mapBranches.GetDouble("mbbcorlow")<140 ) return true;
  else return false;
}


//=================================================================================
bool HHAnalysis::IsHighMass(MapBranches mapBranches){
  if ( mapBranches.GetDouble("pTb1corhigh")>100 && mapBranches.GetDouble("pTb2corhigh")>30 && mapBranches.GetDouble("mbbcorhigh")>90 && mapBranches.GetDouble("mbbcorhigh")<140 ) return true;
  else return false;
}


//==================================================================================
void HHAnalysis::InitialiseHist(TH1D* &hist, string histName, string strVariable){
  TString var=strVariable;
  if (var.Contains("mv2c")) hist= new TH1D (histName.c_str(), "", 200, -1, 1);
  else if (var.BeginsWith("m")) hist= new TH1D (histName.c_str(), "", 100, 0, 1000);
  else if (var.Contains("pT")) hist= new TH1D (histName.c_str(), "", 60, 0, 600);
  else if (var.BeginsWith("d")) hist= new TH1D (histName.c_str(), "", 25, 0, 5);
  else if (var.Contains("eta")) hist= new TH1D (histName.c_str(), "", 1000, -5, 5);
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
  vector <string> vectOpt;
  vector <double> vectExtremalBins;
  vector <TH1*> vectHistTmp;
  TString tmp, name;

  for (unsigned int iCat=0; iCat<m_vectCategories.size(); iCat++){
    for (unsigned int iVar=0; iVar<m_vectVariables.size(); iVar++){
      vectExtremalBins.push_back(15e10); //arbitrary value
      vectExtremalBins.push_back(0);
        plotName="tagcat"+to_string(m_vectCategories[iCat])+"_var"+m_vectVariables[iVar];
      for(auto &it : m_mapHist){ //it.first =histName, it.second=hist
	if ( it.first.find(plotName.c_str())==string::npos ) continue;
	if (m_vectSamples[0]!="ALL"){
	  for (unsigned int iSample=0; iSample<m_vectSamples.size(); iSample++){
	    if ( it.first.find(m_vectSamples[iSample].c_str())==string::npos ) continue;
	  } //end iSample
	}// end if m_vectSamples[0]!="ALL"

	name=it.second->GetName();
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
      drawOpt.AddOption("latexOpt", "0.45 0.85");
      if (m_extraInfo!="") {drawOpt.AddOption("latex", m_extraInfo.c_str()); drawOpt.AddOption("latexOpt", "0.45 0.8");}
      drawOpt.AddOption("xTitle", varName);
      drawOpt.AddOption("normalize", "1");
      drawOpt.AddOption("extension", extension);
      drawOpt.Draw( vectHistTmp );
      
      cout<<m_savePathPlot<<plotName<<"."<<extension<<" has been created.\n";
      vectOpt.clear();
      vectExtremalBins.clear();
      vectHistTmp.clear();
    } //end iVar
  }//end iCat
  
  delete objArrayString; objArrayString=0;
  cout<<"HHAnalysis::DrawDistriForLambdas done.\n";
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
/*void HHAnalysis::DrawCompLONLOForLambdas(TFile *LOFile, TTree *LOTree, TFile *NLOFile, TTree *NLOTree, list<string> vectVariables, string savePath){
}*/
