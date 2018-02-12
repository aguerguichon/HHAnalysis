#include <iostream>
#include <istream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include "boost/multi_array.hpp"

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

using boost::extents;


//==========================================================================================
HHAnalysis::HHAnalysis(){
}


HHAnalysis::~HHAnalysis(){
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
void HHAnalysis::CreateSaveDistri(vector<string> vectInFiles, int selectionType, vector<string> vectVariables, vector<int> vectCategories, string outputFileName){

  cout<<"HHAnalysis::CreateSaveDistri\n";

  map< string, TH1D* > mapHist;
  string *sampleName=new string;
  string histName;
  TBranch *sampleBranch; //!
  double weight{0};
  TFile *inFile{0};
  TTree *inTree{0};
  list<string> listVariables(vectVariables.begin(),vectVariables.end());
  
  listVariables.push_back("tagcat");
  listVariables.push_back("weightMC"); listVariables.push_back("weightvertex"); listVariables.push_back("weightpileup"); listVariables.push_back("weightinit"); listVariables.push_back("Lumi"); listVariables.push_back("LumiMC"); listVariables.push_back("weightlowmass"); listVariables.push_back("weighthighmass");
  listVariables.push_back("pTb1corlow"); listVariables.push_back("pTb2corlow"); listVariables.push_back("pTb1corhigh"); listVariables.push_back("pTb2corhigh"); listVariables.push_back("mbbcorlow"); listVariables.push_back("mbbcorhigh");

  TH1::AddDirectory(kFALSE); 

  for (unsigned int iFile=0; iFile<vectInFiles.size(); iFile++){
    inFile=TFile::Open( (vectInFiles[iFile]).c_str() );
    if (!inFile) throw invalid_argument( "HHAnalysis::CreateSaveDistri: no file provided." );
    inTree=(TTree*)inFile->Get("ntuple");
    if (!inFile) throw invalid_argument( "HHAnalysis::CreateSaveDistri: no tree provided." );

    MapBranches mapBranches; 
    mapBranches.LinkTreeBranches(inTree, 0, listVariables);
    inTree->SetBranchStatus("sample", 1);
    inTree->SetBranchAddress("sample", &sampleName, &sampleBranch);
  
    for (unsigned int iEntry=0; iEntry<inTree->GetEntries(); iEntry++){
      inTree->GetEntry(iEntry);
      if ( !IsEventSelected(selectionType, mapBranches) ) continue;
      weight=GetWeight(selectionType%10, mapBranches);
      if ( vectInFiles[iFile].find("LO")!=string::npos ) weight*=0.34*2.28;
	
      //Create entries of mapHist and fill hists
      for (unsigned int iCat=0; iCat<vectCategories.size(); iCat++){
	if ( mapBranches.GetInt("tagcat")!=vectCategories[iCat] ) continue;
	for (unsigned int iVar=0; iVar<vectVariables.size(); iVar++ ){
	  histName="tagcat"+to_string( mapBranches.GetInt("tagcat") )+"_var"+vectVariables[iVar]+"_sample"+*sampleName;

	  if ( !mapHist.count(histName) ) mapHist[histName]= InitialiseHist(histName, vectVariables[iVar]); 
	  mapHist[histName]->Fill( mapBranches.GetDouble( vectVariables[iVar].c_str() ), weight );
	} // end iVar
      }//end iCat  
    } //end iEntry
    inFile->Close();
  }//end iFile

  //Saving histograms in a root file to be post treated.
  TFile *outputFile=new TFile (outputFileName.c_str(), "RECREATE");
  for(auto &it : mapHist) it.second->Write();
  outputFile->Close();
  
  delete sampleName; sampleName=0;
  delete inFile; inFile=0; 
  delete outputFile; outputFile=0;
  cout<<"HHAnalysis::CreateSaveDistri done.\n";
}


//==================================================================================
bool HHAnalysis::IsEventSelected(int selectionType, MapBranches mapBranches){

  if (selectionType%10 == 0 || selectionType%10 == 1) return true; //if no selection or isPassed (ie skimmed)

  else if (selectionType%10 == 2){ //lowmass
    if ( !IsLowMass(mapBranches) ) return false;
    switch(selectionType/10){
    case 1: { return true; break;} //basic low mass selection
    case 2: { if (mapBranches.GetDouble("mgamgam")>120.39 && mapBranches.GetDouble("mgamgam")<129.79 ) return true; else return false; break;} //adding mgamgam cut
    case 3: { if (mapBranches.GetDouble("mbbgamgam")<350  ) return true; else return false; break;} //adding mhh cut
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for low mass selection: "+selectionType);
    }
  }

  else if (selectionType%10 == 3){ //highmass
    if ( !IsHighMass(mapBranches) ) return false;
    switch(selectionType/10){
    case 1: { return true; break;} //basic high mass selection
    case 2: { if (mapBranches.GetDouble("mgamgam")>120.79 && mapBranches.GetDouble("mgamgam")<129.39  ) return true; else return false; break;} //adding mgamgam cut
    case 3: { if (mapBranches.GetDouble("mbbgamgam")>350  ) return true; else return false; break;}
    default: throw invalid_argument("HHAnalysis::IsEventSelected: the following selectionType does not exist for high mass selection: "+selectionType);
    }
  }
}

//==================================================================================
double HHAnalysis::GetWeight(int weightType, MapBranches mapBranches){
  double weight{0};
  switch (weightType){
  case 0: {weight=mapBranches.GetDouble("weightMC")*mapBranches.GetDouble("weightvertex")*mapBranches.GetDouble("weightpileup")*(mapBranches.GetDouble("Lumi")/mapBranches.GetDouble("LumiMC") ); break;} //for no selection cut at all
  case 1:{ weight=mapBranches.GetDouble("weightinit")*(mapBranches.GetDouble("Lumi")/mapBranches.GetDouble("LumiMC") ); break;} //for selection cut at isPassed (equivalent to skimmed files)
  case 2:{weight=mapBranches.GetDouble("weightlowmass"); break;}
  case 3:{weight=mapBranches.GetDouble("weighthighmass"); break;}

  default: weight=weight=mapBranches.GetDouble("weightMC")*mapBranches.GetDouble("weightvertex")*mapBranches.GetDouble("weightpileup")*(mapBranches.GetDouble("Lumi")/mapBranches.GetDouble("LumiMC") ); //no selection at all
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
TH1D* HHAnalysis::InitialiseHist(string histName, string strVariable){
  TString var=strVariable;
  TH1D* hist{0};
  if (var.Contains("mv2c")) hist= new TH1D (histName.c_str(), "", 200, -1, 1);
  else if (var.BeginsWith("m")) hist= new TH1D (histName.c_str(), "", 100, 0, 1000);
  else if (var.Contains("pT")) hist= new TH1D (histName.c_str(), "", 60, 0, 600);
  else if (var.BeginsWith("d")) hist= new TH1D (histName.c_str(), "", 25, 0, 5);
  else if (var.Contains("eta") || var.Contains("phi")) hist= new TH1D (histName.c_str(), "", 1000, -5, 5);
  else hist= new TH1D (histName.c_str(), "", 2000, -1000, 1000);
  return hist;
}

//==================================================================================
void HHAnalysis::DrawDistriForLambdas(TFile *inFile, vector<string> vectVariables, vector<int> vectCategories, vector<string> vectSamples, string path, string extension){

  cout<<"HHAnalysis::DrawDistriForLambdas done.\n";
  TKey *key{0};
  TH1D *hist{0};
  TClass *cl{0};
  TObjArray *objArrayString{0};
  string histName, plotName, catName, sampleName, varName, legLatex;
  map< string, vector<TH1*> > mapHist;
  vector <string> vectOpt;
  vector <double> vectExtremalBins;
  vector <TH1*> vectHistTmp;
  TString tmp, name;

  TIter next( inFile->GetListOfKeys() );
  while ( (key = (TKey *) next()) ) { //iteration on the keys of the root file
    cl = gROOT->GetClass(key->GetClassName()); //getting the class of each key
    if (cl->InheritsFrom("TH1D")){ //loop over all histograms
      hist = (TH1D*)key->ReadObj();
      histName=hist->GetName();
      for (unsigned int iCat=0; iCat<vectCategories.size(); iCat++){
	for (unsigned int iVar=0; iVar<vectVariables.size(); iVar++){
	  plotName="tagcat"+to_string(vectCategories[iCat])+"_var"+vectVariables[iVar];
	  if ( !mapHist.count(plotName) ) mapHist[plotName];
	  for (unsigned int iSample=0; iSample<vectSamples.size(); iSample++){
	    if ( histName.find(plotName.c_str())==string::npos ) continue;
	    if ( histName.find(vectSamples[iSample].c_str())==string::npos && vectSamples[iSample]!="ALL") continue;
	    mapHist[plotName].push_back(hist);	    
	  }//end iSample
	} //end it vectVariables
      }//end iCat
    } //if Inherits
  } //end while

  for(auto &it : mapHist) {
    vectExtremalBins.push_back(0);
    vectExtremalBins.push_back(0);
    for (auto &h : it.second) {
      name=h->GetName();
      objArrayString = name.Tokenize("_");
      for(int iString=0; iString<objArrayString->GetEntriesFast(); iString++){
	tmp = ((TObjString*)objArrayString->At(iString))->GetString();
	if ( tmp.Contains("tagcat") ){ catName=tmp.ReplaceAll("tagcat","");}
	if ( tmp.Contains("var") ) varName=tmp.ReplaceAll("var","");
	if ( tmp.Contains("sample") ) sampleName=tmp.ReplaceAll("sample","");
      }
      vectOpt.push_back( ("legend="+ sampleName ).c_str() );
      if ( ReturnExtremalBins(h)[0]<vectExtremalBins[0] ) vectExtremalBins[0]=ReturnExtremalBins(h)[0];
      if ( ReturnExtremalBins(h)[1]>vectExtremalBins[1] ) vectExtremalBins[1]=ReturnExtremalBins(h)[1];
    } //end loop over hist

    DrawOptions drawOpt;
    vectHistTmp=it.second;
    drawOpt.AddOption(vectOpt);
    drawOpt.AddOption("rangeUserX", (to_string(vectExtremalBins[0])+" "+to_string(vectExtremalBins[1])).c_str() );
    drawOpt.AddOption("outName", (path+it.first).c_str() );
    drawOpt.AddOption("legendPos", "0.75 0.9");
    drawOpt.AddOption("latex", ("Category "+ catName +" b-tag").c_str() );
    drawOpt.AddOption("latexOpt", "0.45 0.8");
    drawOpt.AddOption("xTitle", varName);
    drawOpt.AddOption("normalize", "1");
    drawOpt.AddOption("extension", extension);
    drawOpt.Draw( vectHistTmp );
    vectOpt.clear();
    vectHistTmp.clear();
    vectExtremalBins.clear();

  }//end over key of map ie plots

  delete key; key=0;
  delete cl; cl=0;
  delete objArrayString; objArrayString=0;
  delete hist; hist=0;
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
void HHAnalysis::DrawCompLONLOForLambdas(TFile *LOFile, TTree *LOTree, TFile *NLOFile, TTree *NLOTree, list<string> vectVariables, string savePath){
}