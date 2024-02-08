//Author Ezekiel Wertz
//02/01/2024
//Purpose: Parsing Script for LD2 data to produce output histograms for later analysis


#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "../../src/parse_config.C"
#include "../../src/data_object.C"

//Main
void data_elastic_parse(const char *setup_file_name){

  //Define a clock to check macro processing time
  TStopwatch *watch = new TStopwatch(); 
  watch->Start( kTRUE );

  //parse object to get in the information that The One Config file has and is manipulated
  parse_config mainConfig(setup_file_name); 
  //mainConfig.printDataYields();
  
  //store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
  vector<int> runNums = mainConfig.getRunNums();
  TCut globalcut = mainConfig.getGlobalCut();
  TString exp = mainConfig.getExp();
  TString kin = mainConfig.getKin();
  TString kinematic_file = mainConfig.getKinFileName();
  TString data_map = mainConfig.getDataFileName();
  int sbs_field = mainConfig.getSBSField();
  TString target = mainConfig.getTarg();
  double W2_mean = mainConfig.getW2Mean();
  double W2_sigma = mainConfig.getW2Sigma();
  double W2_sigfac = mainConfig.getW2SigFac();
  double dxO_n = mainConfig.get_dxOn();
  double dy0_n = mainConfig.get_dyOn();

  int W2_low = W2_mean - W2_sigfac*W2_sigma;
  int W2_high = W2_mean + W2_sigfac*W2_sigma;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;	

}//end Main
