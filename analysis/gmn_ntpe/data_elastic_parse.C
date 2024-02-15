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
#include "TString.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "../../src/utility.C"
#include "../../include/physics_constants.h"
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/parse_config.C"
#include "../../src/data_object.C"
#include "../../src/cuts.C"

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
  TString pass = mainConfig.getPass();
  int sbs_field = mainConfig.getSBSField();
  int maxtracks = mainConfig.getMAXNTRACKS();
  TString target = mainConfig.getTarg();
  double W2_mean = mainConfig.getW2Mean();
  double W2_sigma = mainConfig.getW2Sigma();
  double W2_sigfac = mainConfig.getW2SigFac();
  double dxO_n = mainConfig.get_dxOn();
  double dyO_n = mainConfig.get_dyOn();
  double dxsig_n = mainConfig.get_dxsign();
  double dysig_n = mainConfig.get_dysign();
  double dxO_p = mainConfig.get_dxOp();
  double dyO_p = mainConfig.get_dyOp();
  double dxsig_p = mainConfig.get_dxsigp();
  double dysig_p = mainConfig.get_dysigp();
  double dx_pn = mainConfig.get_dxpn();
  double dx_low = mainConfig.get_dxLow();
  double dx_high = mainConfig.get_dxHigh();
  double dy_low = mainConfig.get_dyLow();
  double dy_high = mainConfig.get_dyHigh();
  int useAlshield = mainConfig.getAlshield();
  double dxsig_n_fac = mainConfig.get_dxSignFac();
  double dxsig_p_fac = mainConfig.get_dxSigpFac();
  double dysig_n_fac = mainConfig.get_dySignFac();
  double dysig_p_fac = mainConfig.get_dySigpFac();

  int W2_low = W2_mean - W2_sigfac*W2_sigma;
  int W2_high = W2_mean + W2_sigfac*W2_sigma;

  //store all important kinematic info in local variables
  kinematic_obj myKin(kinematic_file,kin);
  double EBeam = myKin.getBeamEnergy();
  double hcaldist = myKin.getHCalDist();
  double sbsdist = myKin.getSBSDist();
  double bbtheta = myKin.getBBAngle_Rad();
  double hcaltheta = myKin.getHCalAngle_Rad();
  
  //setup hcal active area with bounds that match database depending on pass
  vector<double> hcalaa = cuts::hcal_ActiveArea_data(1,1,pass);

  //proceed with active area/fid cuts. // Also put output file name makers somewhere

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;	

}//end Main
