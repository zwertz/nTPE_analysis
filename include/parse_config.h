#ifndef PARSE_CONFIG_H
#define PARSE_CONFIG_H

//class to handle parsing my config files for all types of analysis
//Author: Ezekiel Wertz

#include "TString.h"
#include "TCut.h"
#include <vector>

class parse_config{
  private:
  
  //common data analysis
  TString Exp,kin,data_file_name,kinematic_file_name,targ,pass,Data_file,fitopt;
  int SBS_field,useAlshield,MAXNTRACKS, e_method,hcalnclusmin;  
  double dxO_n,dyO_n,dxsig_n,dysig_n,dxO_p,dyO_p,dxsig_p,dysig_p,dx_pn,W2_low,W2_high,dx_low,dx_high,dy_low,dy_high,dxsig_n_fac,dxsig_p_fac,dysig_n_fac,dysig_p_fac,coin_mean,coin_sigma,coin_sig_fac,coin_profile_sig,dysig_cut_fac;
  TCut globalcut;
  vector<int> runnums;


  //common MC data analsyis
  TString proton_root_file, neutron_root_file, MC_file,rootfile_dir,histfile_dir,replay_type,partial_name_p,partial_name_n;
  double sf,Ntried_override,luminosity_override,genvol_override;
  bool sync_jobs,mc_override;

  //For MC HCal efficiency
  double hcalemin,proton_thresh_fac, neutron_thresh_fac,num_bin, pmin, pmax, Emin, Emax;

  //For data HCal Efficiency analysis
  double thetapq_low,thetapq_high,W2fitmax,W2fitmaxwide,binfac,hbinfac,fidx_min,fidx_max,fidy_min,fidy_max,spot_sig;

  //For stability studies
  TString EnergyCut, TrackQualityCut, TargetVertexCut, W2Cut, FidXCut, FidYCut, dyCut, eOverpCut, HCal_Shower_atime_Cut, HCal_Energy_Cut, OpticsCut, ProtonSpotCut, NeutronSpotCut, isProtonCut, isNeutronCut;


public:
  //Constructor
  parse_config(const char *setup_file_name);

  //destructor
  //no dynamic memory or pointers so no problem
  ~parse_config();

  //Get function for every private variable
  //Probably don't need set functions, once the constructor does the initial setup that should be only necessary

  TString getExp();
 
  TString getKin();

  TString getDataFileName();

  TString getKinFileName();

  TString getTarg();

  TString getProtonFileName();

  TString getNeutronFileName();

  TString getMCFileName();

  TString getPass();

  TString getRootFileDir();

  TString getHistFileDir();

  TString getReplayType();

  TString getPartialNameP();

  TString getPartialNameN();

  TString getDataFile();

  TString getFitOpt();

  TString getEnergyCut();

  TString getTrackQualityCut();

  TString getTargetVertexCut();

  TString getW2Cut();

  TString getFidXCut();

  TString getFidYCut();

  TString getdyCut();

  TString geteOverpCut();

  TString getHCal_Shower_atime_Cut();

  TString getHCal_Energy_Cut();

  TString getOpticsCut();

  TString getProtonSpotCut();

  TString getNeutronSpotCut();

  TString getisProtonCut();

  TString getisNeutronCut();

  int getSBSField();

  int getAlshield();

  int getMAXNTRACKS();

  int get_emethod();

  int get_HCalNclusMin();

  double get_dxOn();

  double get_dyOn();

  double get_dxsign();

  double get_dysign();

  double get_dxOp();

  double get_dyOp();

  double get_dxsigp();

  double get_dysigp();

  double getW2Low();

  double getW2High();

  double get_dxLow();

  double get_dxHigh();

  double get_dyLow();

  double get_dyHigh();

  double get_dxSignFac();

  double get_dxSigpFac();

  double get_dySignFac();

  double get_dySigpFac();

  double get_dySigCutFac();

  double get_dxpn();
  
  double getCoinMean();

  double getCoinSig();

  double getCoinSigFac();

  double getCoinProfSig();

  double getHCaleMin();

  double getProtonThreshFac();

  double getNeutronThreshFac();

  double getNumBin();

  double getPmax();

  double getPmin();

  double getEmin();

  double getEmax();

  double getThetapqLow();  

  double getThetapqHigh();

  double getW2FitMax();

  double getW2FitMaxWide();
 
  double getBinFac();

  double getHBinFac();

  double getNTriedOverride();

  double getLumiOverride();

  double getVolOverride();

  double get_sf();

  double get_fidxmin();

  double get_fidxmax();

  double get_fidymin();

  double get_fidymax();

  double get_spotsig();

  bool get_syncJobs();

  bool get_MCOverride();

  TCut getGlobalCut();

  vector<int> getRunNums();
  
  void printRunNums();

  //only work on LD2 data 
  void printDataYields();

  //only work on HCal Efficiency LH2 data
  void printDataHCalEff();

   //only work on MC LD2  
   void printMCYields();
   
   //only work on HCal Efficiency MC LH2 
   void printMCHCalEff();

};//end of class
#endif
