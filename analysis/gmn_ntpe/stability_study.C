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
  TString Data_input_file_name = mainConfig.getDataFile();
  TString MC_input_file_name = mainConfig.getMCFileName();
  TString EnergyCut = mainConfig.getEnergyCut();
  TString TrackQualityCut = mainConfig.getTrackQualityCut();
  TString TargetVertexCut = mainConfig.getTargetVertexCut();
  TString W2Cut = mainConfig.getW2Cut();
  TString FidXCut = mainConfig.getFidXCut();
  TString FidYCut = mainConfig.getFidYCut();
  TString dyCut = mainConfig.getdyCut();
  TString eOverpCut = mainConfig.geteOverpCut();
  TString HCal_Energy_Cut = mainConfig.getHCal_Energy_Cut();
  TString HCal_Shower_atime_Cut = mainConfig.getHCal_Shower_atime_Cut();
  TString OpticsCut = mainConfig.getOpticsCut();
  TString ProtonSpotCut = mainConfig.getProtonSpotCut();
  TString NeutronSpotCut = mainConfig.getNeutronSpotCut();
  TString isProtonCut = mainConfig.getisProtonCut();
  TString isNeutronCut = mainConfig.getisNeutronCut();
  double W2fitmax = mainConfig.getW2FitMax();
  double hbinfac = mainConfig.getHBinFac();
  double binfac = mainConfig.getBinFac();

  //Setup common histogram information
  double dx_hist_low = exp_constants::hcalposXi_mc;
  double dx_hist_high = exp_constants::hcalposXf_mc;
  double dx_hist_bin = exp_constants::hcal_vrange*hbinfac;
  double W2_hist_low = 0.0;
  double W2_hist_high = W2fitmax;
  double W2_hist_bin = W2fitmax*binfac;

  //Setup output file
  TString outfile = utility::makeOutputFileName_Stability(exp,pass,kin,sbs_field,target);
  TString reportfile = utility::makeReportFileName_Stability(exp,pass,kin,sbs_field,target);
  TFile *fout = new TFile(outfile,"RECREATE");

  //Intialize the TChain for data
  TChain *C_data = new TChain("Parse");
  C_data->Add(Data_input_file_name.Data());

  //Intialize the TChain for mc
  TChain *C_mc = new TChain("Parse");
  C_mc->Add(MC_input_file_name.Data());

  //Initialize data and MC string flags to separate the histograms
  TString data_flag_string = "data";
  TString mc_p_flag_string = "mc_p";
  TString mc_n_flag_string = "mc_n";

  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data ps
  TString ps_e_study_data = TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC ps 
  TString ps_e_study_mc = TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut;
  TString ps_e_study_mc_p = ps_e_study_mc + "&&" + isProtonCut;
  TString ps_e_study_mc_n = ps_e_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  cutvar psVar_data("BBps_e", ps_e_study_data, data_flag_string, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar psVar_mc_p("BBps_e", ps_e_study_mc_p, mc_p_flag_string, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar psVar_mc_n("BBps_e", ps_e_study_mc_n, mc_n_flag_string, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* ps_e_dxhist_data = psVar_data.get2DdxCutHisto();
  TH2D* ps_e_dxhist_mc_p = psVar_mc_p.get2DdxCutHisto();
  TH2D* ps_e_dxhist_mc_n = psVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_psVar_data = psVar_data.sliceAndProjectHisto_xMinxMax(ps_e_dxhist_data,psVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_psVar_mc_p = psVar_mc_p.sliceAndProjectHisto_xMinxMax(ps_e_dxhist_mc_p,psVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_psVar_mc_n = psVar_mc_n.sliceAndProjectHisto_xMinxMax(ps_e_dxhist_mc_n,psVar_mc_n.getAxisTitle().Data(),"dx");




  fout->Write();
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
