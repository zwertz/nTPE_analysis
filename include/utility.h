#ifndef UTILITY_H
#define UTILITY_H


//Author Ezekiel Wertz
//A location to hold useful functions that do a task

#include "TString.h"

namespace utility{

  string out_dir_temp = "/work/halla/sbs/ewertz/nTPE_analysis/output";
  TString output_directory = TString(out_dir_temp);

  bool check_number(const char *myChar);

  double DegToRad(double myNum);

  double RadToDeg(double datRad);

  TString intToTString(int datInt);

  TString getOutputDir();

  TString makeOutputFileNameParse(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  TString makeOutputFileName_MCParse(TString exp, TString Kin, int SBS_field);

  TString makeYieldReportFileName(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  TString makeOutputFileNameHCalEffParse(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  TString makeOutputFileName_HCalBS(TString exp, TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event);

  TString makeReportFileName_HCalEff(TString exp, TString pass, TString Kin, int SBS_field,TString target);

  TString makeReportFileName_HCalEffBS(TString exp, TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event);

  TString makeOutputFileName_MCsf(TString exp, TString Kin, int SBS_field,double sf,TString target);
 
  TString makeOutputFileName_HCalEffMC(TString exp, TString Kin);

  vector<string> findHistFiles(TString replay_type,TString histDirectory,TString partialName);

  void matchMCFiles(TString replay_type,vector<string>& histFiles,vector<string>& rootFiles, TString rootDirectory);

  void syncJobNumbers(vector<string>& proton_vec,vector<string>& neutron_vec);

}
#endif
