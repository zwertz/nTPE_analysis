//Author: Ezekiel Wertz
//07/31/2024
//Purpose: A script to evaluate the stability of various analysis cuts. Ultimately needed to justify final cut decisions. Presumably will evolve to to justify systematic error analysis.

//The exact ordering of this matters. ROOT for some reason cannot handle calls for files that have already been included.
#include "TFile.h"
#include <vector>
#include "../../src/utility.C"
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/fits.C"
#include "../../src/physics.C"
#include "../../src/parse_config.C"
#include "../../src/plots.C"
#include "../../src/calc_FFs_RCS_obj.C"
#include "../../src/cutvar.C"
#include "../../src/stability_analysis.C"

//This must live in the script so it has knowledge of the histograms for interpolation
//A fit function for dx slices. Primarily used in stability_analysis class
TH1D* histo_p;
TH1D* histo_n;

double mc_p_n_poly2_slice_fit(double *x, double *param){

        //MC float parameters
        double proton_scale = param[0];
        double R_sf = param[1];

        double dx_shift_p = param[2];
        double dx_shift_n = param[3];

        //Apply the shifts before any interpolation
        double proton = histo_p->Interpolate(x[0] - dx_shift_p);
        double neutron = histo_n->Interpolate(x[0] - dx_shift_n);

        //The total function is the proton and neutron peaks + a 2nd order polynomial for background
        double fit = proton_scale*( R_sf*neutron + proton) + fits::poly2_fit(x, &param[4]);
        return fit;
}

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
  TString psVar_string = "BBps_e";
  cutvar psVar_data(psVar_string, ps_e_study_data, data_flag_string, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar psVar_mc_p(psVar_string, ps_e_study_mc_p, mc_p_flag_string, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar psVar_mc_n(psVar_string, ps_e_study_mc_n, mc_n_flag_string, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
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

  //Fit type
  const char* FitType = "mc_p_n_poly2_slice_fit";

  //Maybe not best practice. Need histo_p and histo_n to be non null. So they get past the function call. Ultimately these will be manipulated in the function call
  histo_p = (TH1D*)slice_psVar_mc_p[0]->Clone("Test_histo_p");
  histo_n = (TH1D*)slice_psVar_mc_n[0]->Clone("Test_histo_n");
  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis psVar_stability(psVar_data,psVar_mc_p,psVar_mc_n,slice_psVar_data,slice_psVar_mc_p,slice_psVar_mc_n,FitType,histo_p,histo_n);


  fout->Write();
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
