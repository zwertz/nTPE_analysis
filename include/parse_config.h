#ifndef PARSE_CONFIG_H
#define PARSE_CONFIG_H

//class to handle parsing my config files for all types of analysis
//Author: Ezekiel Wertz


class parse_config{
  private:
  
  //common data analysis
  TString Exp,kin,data_file_name,kinematic_file_name,targ,pass;
  int SBS_field,useAlshield,MAXNTRACKS;  
  double dxO_n,dyO_n,dxsig_n,dysig_n,dxO_p,dyO_p,dxsig_p,dysig_p,dx_pn,W2_mean,W2_sigma,W2_sigfac,dx_low,dx_high,dy_low,dy_high,dxsig_n_fac,dxsig_p_fac,dysig_n_fac,dysig_p_fac,coin_mean,coin_sigma,coin_sig_fac;
  TCut globalcut;
  vector<int> runnums;


  //common MC data analsyis
  TString proton_root_file, neutron_root_file, MC_file;
  double sf;

  //For MC HCal efficiency
  double hcalemin,proton_thresh_fac, neutron_thresh_fac,num_bin, pmin, pmax, Emin, Emax;

  //For data HCal Efficiency analysis
  double thetapq_low,thetapq_high,W2fitmax,W2fitmaxwide,binfac;

public:
  //Constructor
  parse_config(const char *setup_file_name);

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

  int getSBSField();

  int getAlshield();

  int getMAXNTRACKS();

  double get_dxOn();

  double get_dyOn();

  double get_dxsign();

  double get_dysign();

  double get_dxOp();

  double get_dyOp();

  double get_dxsigp();

  double get_dysigp();

  double getW2Mean();

  double getW2Sigma();

  double getW2SigFac();

  double get_dxLow();

  double get_dxHigh();

  double get_dyLow();

  double get_dyHigh();

  double get_dxSignFac();

  double get_dxSigpFac();

  double get_dySignFac();

  double get_dySigpFac();

  double get_dxpn();
  
  double getCoinMean();

  double getCoinSig();

  double getCoinSigFac();

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

  double get_sf();

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
