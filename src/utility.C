#include "../include/utility.h"
#include <iostream>
#include <string>
#include "TMath.h"

//Author Ezekiel Wertz
//Implementation for useful task functions

namespace utility{

  bool check_number(const char *myChar){
   bool mybool=false;
  	for(int i=0; i< strlen(myChar); i++){
		//So we ever find not a digit. Immediately return mybool=false. So we have a string
		//cout << mybool << endl;
		if(isdigit(myChar[i])==false){
		return mybool;
		}
  	}
  //If we get through the for loop we never found not a digit. Make it true and return
  mybool=true;
  return mybool;
  }//end check_number

  //convert from degrees to radians
  double DegToRad(double myNum){
  return (myNum * (TMath::Pi()/180.0));
  }
  double RadToDeg(double datRad){
  return ((datRad * 180.0)/(TMath::Pi()));
  }

  //convert an int to a TString. Just a helper function
  TString intToTString(int datInt){
  stringstream ss;
  ss << datInt;
  std::string datInt_temp;
  ss >> datInt_temp;
  TString datInt_string(datInt_temp.c_str());
  //cout << datInt_string << endl;
  return datInt_string;
  }

  //This is a helper function that will tell you the default directory for all analysis output
  TString getOutputDir(){
  return output_directory;
  }

  //Helper function to make output file name for parsed data root file
  TString makeOutputFileNameParse(TString exp,TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/Zeke_parsed_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //TODO:Will need a function to make file name for compared yield data and MC

  //Helper function to make output file name for yield or ratio information
  TString makeYieldReportFileName(TString exp,TString pass,  TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/yieldRep_Zeke_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(), pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }
  
  //Helper function to make output file name for HCal Efficiency from data
  TString makeOutputFileNameHCalEffParse(TString exp,TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/Zeke_HCalEff_parsed_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function to make output file name for HCal Efficiency from MC
  TString makeOutputFileName_HCalEffMC(TString exp, TString Kin){
  TString outfile = Form("%s/Zeke_HCalEff_MC_%s_%s",(getOutputDir()).Data(),exp.Data(),Kin.Data());
 
  return outfile;
  }

  //TODO: Will need a function to make file name for compare HCal Efficiency data with MC

  //Helper function to make output file for HCal Efficiency Bootstrap
  TString makeOutputFileName_HCalBS(TString exp,TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event){
  TString outfile = Form("%s/Zeke_HCalEff_bootstrap_%s_%s_%s_%s_%i_event_%i_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,first_event,last_event);
  //cout << outfile << endl;
  return outfile;
  }
  //Helper function for making file name for HCal Efficiency report
  TString makeReportFileName_HCalEff(TString exp,TString pass,  TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/efficiencyRep_Zeke_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(), pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function for making file name for HCal Efficiency report with boot strap
  TString makeReportFileName_HCalEffBS(TString exp,TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event){
  TString outfile = Form("%s/efficiencyRep_Zeke_%s_%s_%s_%s_%i_event_%i_%i.txt",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,first_event,last_event);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function for output file name for scale field study for MC
  TString makeOutputFileName_MCsf(TString exp, TString Kin, int SBS_field,double sf,TString target){
  TString outfile = Form("%s/MC_Zeke_%s_%s_%s_%i_sf%.1f.root",(getOutputDir()).Data(),exp.Data(),Kin.Data(),target.Data(),SBS_field,sf);
  //cout << outfile << endl;
  return outfile;
  }

} //end namespace
