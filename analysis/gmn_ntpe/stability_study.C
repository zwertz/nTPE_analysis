//Author: Ezekiel Wertz
//07/31/2024
//Purpose: A script to evaluate the stability of various analysis cuts. Ultimately needed to justify final cut decisions. Presumably will evolve to to justify systematic error analysis.

#include "TFile.h"
#include <vector>
#include "../../src/utility.C"
#include "../../src/exp_constants.C"
#include "../../src/fits.C"
#include "../../src/parse_config.C"
#include "../../src/plots.C"
#include "../../src/cutvar.C"
void stability_study(const char *setup_file_name){

  //Define a clock to check macro processing time
  TStopwatch *watch = new TStopwatch();
  watch->Start( kTRUE );

  //parse object to get in the information that The One Config file has and is manipulated
  parse_config mainConfig(setup_file_name);

  //store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
  TString exp = mainConfig.getExp();
  TString kin = mainConfig.getKin();
  int sbs_field = mainConfig.getSBSField();
  TString pass = mainConfig.getPass();
  TString target = mainConfig.getTarg();
  TString cutvar_string = mainConfig.getCutVar();
  TString Data_input_file_name = mainConfig.getDataFile();
  TString MC_input_file_name = mainConfig.getMCFileName();
  

  //setup input file info
  TFile *data_file = new TFile(Data_input_file_name.Data());
  TFile *mc_file = new TFile(MC_input_file_name.Data());

  //Setup output file
  TString outfile = utility::makeOutputFileName_Stability(exp,pass,kin,sbs_field,target,cutvar_string);
  TString reportfile = utility::makeReportFileName_Stability(exp,pass,kin,sbs_field,target,cutvar_string);
  TFile *fout = new TFile(outfile,"RECREATE");

  //setup cutvar object, central to this stability study
  cutvar datVar(cutvar_string);



  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
