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
#include "../../src/fit_histogram.C"
#include "../../src/cutvar.C"
#include "../../src/stability_analysis.C"

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
  TString TrackHitsCut = mainConfig.getTrackHitsCut();
  TString TrackQualityCut = mainConfig.getTrackQualityCut();
  TString TargetVertexCut = mainConfig.getTargetVertexCut();
  TString W2Cut = mainConfig.getW2Cut();
  TString FidXCut = mainConfig.getFidXCut();
  TString FidYCut = mainConfig.getFidYCut();
  TString dyCut = mainConfig.getdyCut();
  TString eOverpCut = mainConfig.geteOverpCut();
  TString HCal_Energy_Cut = mainConfig.getHCal_Energy_Cut();
  TString HCal_Shower_atime_Cut = mainConfig.getHCal_Shower_atime_Cut();
  TString OpticsCut_x = mainConfig.getOpticsCutX();
  TString OpticsCut_y = mainConfig.getOpticsCutY();
  TString ProtonSpotCut = mainConfig.getProtonSpotCut();
  TString NeutronSpotCut = mainConfig.getNeutronSpotCut();
  TString isProtonCut = mainConfig.getisProtonCut();
  TString isNeutronCut = mainConfig.getisNeutronCut();
  TString left_right = mainConfig.get_left_right();
  double W2fitmax = mainConfig.getW2FitMax();
  double hbinfac = mainConfig.getHBinFac();
  double binfac = mainConfig.getBinFac();
  int slice_mode = mainConfig.get_slice_mode();

  //Setup common histogram information
  double dx_hist_low = exp_constants::hcalposXi_mc;
  double dx_hist_high = exp_constants::hcalposXf_mc;
  double dx_hist_bin = exp_constants::hcal_vrange*hbinfac;
  double W2_hist_low = 0.0;
  double W2_hist_high = W2fitmax;
  double W2_hist_bin = W2fitmax*binfac;

  //Setup output file
  TString outfile = utility::makeOutputFileName_Stability(exp,pass,kin,sbs_field,target,slice_mode,left_right);
  //TString reportfile = utility::makeReportFileName_Stability(exp,pass,kin,sbs_field,target);
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


  ////////////////////////////////////////////////////////////////////////////////////////////
  //PS Energy Stability Analysis

  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data ps
  TString ps_e_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y + "&&"+ HCal_Shower_atime_Cut;
  //For MC ps 
  TString ps_e_study_mc = TrackHitsCut + "&&" + TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString ps_e_study_mc_p = ps_e_study_mc + "&&" + isProtonCut;
  TString ps_e_study_mc_n = ps_e_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString psVar_string = "BBps_e";
  cutvar psVar_data(psVar_string, ps_e_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar psVar_mc_p(psVar_string, ps_e_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar psVar_mc_n(psVar_string, ps_e_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
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
  const char* FitType = "fitFullShift_polyBG";
  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis psVar_stability(psVar_data,psVar_mc_p,psVar_mc_n,slice_psVar_data,slice_psVar_mc_p,slice_psVar_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //number of BBGEM hits on a track Stability Analysis

  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data BBGEM hits
  TString BBgem_nhits_study_data = EnergyCut + "&&" + TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y  + "&&" + HCal_Shower_atime_Cut;
  //For MC BBGEM hits
  TString BBgem_nhits_study_mc = EnergyCut + "&&" +TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString BBgem_nhits_study_mc_p = BBgem_nhits_study_mc + "&&" + isProtonCut;
  TString BBgem_nhits_study_mc_n = BBgem_nhits_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString BBgem_nhitsVar_string = "BBgem_nhits";
  cutvar BBgem_nhitsVar_data(BBgem_nhitsVar_string, BBgem_nhits_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar BBgem_nhitsVar_mc_p(BBgem_nhitsVar_string, BBgem_nhits_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar BBgem_nhitsVar_mc_n(BBgem_nhitsVar_string, BBgem_nhits_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* BBgem_nhits_dxhist_data = BBgem_nhitsVar_data.get2DdxCutHisto();
  TH2D* BBgem_nhits_dxhist_mc_p = BBgem_nhitsVar_mc_p.get2DdxCutHisto();
  TH2D* BBgem_nhits_dxhist_mc_n = BBgem_nhitsVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_BBgem_nhitsVar_data = BBgem_nhitsVar_data.sliceAndProjectHisto_xMinxMax(BBgem_nhits_dxhist_data,BBgem_nhitsVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_BBgem_nhitsVar_mc_p = BBgem_nhitsVar_mc_p.sliceAndProjectHisto_xMinxMax(BBgem_nhits_dxhist_mc_p,BBgem_nhitsVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_BBgem_nhitsVar_mc_n = BBgem_nhitsVar_mc_n.sliceAndProjectHisto_xMinxMax(BBgem_nhits_dxhist_mc_n,BBgem_nhitsVar_mc_n.getAxisTitle().Data(),"dx");

  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis BBgem_nhitsVar_stability(BBgem_nhitsVar_data,BBgem_nhitsVar_mc_p,BBgem_nhitsVar_mc_n,slice_BBgem_nhitsVar_data,slice_BBgem_nhitsVar_mc_p,slice_BBgem_nhitsVar_mc_n,FitType);


  ////////////////////////////////////////////////////////////////////////////////////////////
  //BBgem Track Chi2/ndf  Stability Analysis

  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data BBGEM Track Chi2
  TString BBgem_chi2ndf_study_data = EnergyCut + "&&" + TrackHitsCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC BBGEM Track Chi2
  TString BBgem_chi2ndf_study_mc = EnergyCut + "&&" +TrackHitsCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString BBgem_chi2ndf_study_mc_p = BBgem_chi2ndf_study_mc + "&&" + isProtonCut;
  TString BBgem_chi2ndf_study_mc_n = BBgem_chi2ndf_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString BBgem_chi2ndfVar_string = "BBgem_chi2ndf";
  cutvar BBgem_chi2ndfVar_data(BBgem_chi2ndfVar_string, BBgem_chi2ndf_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar BBgem_chi2ndfVar_mc_p(BBgem_chi2ndfVar_string, BBgem_chi2ndf_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar BBgem_chi2ndfVar_mc_n(BBgem_chi2ndfVar_string, BBgem_chi2ndf_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* BBgem_chi2ndf_dxhist_data = BBgem_chi2ndfVar_data.get2DdxCutHisto();
  TH2D* BBgem_chi2ndf_dxhist_mc_p = BBgem_chi2ndfVar_mc_p.get2DdxCutHisto();
  TH2D* BBgem_chi2ndf_dxhist_mc_n = BBgem_chi2ndfVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_BBgem_chi2ndfVar_data = BBgem_chi2ndfVar_data.sliceAndProjectHisto_xMinxMax(BBgem_chi2ndf_dxhist_data,BBgem_chi2ndfVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_BBgem_chi2ndfVar_mc_p = BBgem_chi2ndfVar_mc_p.sliceAndProjectHisto_xMinxMax(BBgem_chi2ndf_dxhist_mc_p,BBgem_chi2ndfVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_BBgem_chi2ndfVar_mc_n = BBgem_chi2ndfVar_mc_n.sliceAndProjectHisto_xMinxMax(BBgem_chi2ndf_dxhist_mc_n,BBgem_chi2ndfVar_mc_n.getAxisTitle().Data(),"dx");

  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis BBgem_chi2ndfVar_stability(BBgem_chi2ndfVar_data,BBgem_chi2ndfVar_mc_p,BBgem_chi2ndfVar_mc_n,slice_BBgem_chi2ndfVar_data,slice_BBgem_chi2ndfVar_mc_p,slice_BBgem_chi2ndfVar_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //Vertex Stability Analysis

  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data vertex
  TString vert_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC vertex 
  TString vert_study_mc = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString vert_study_mc_p = vert_study_mc + "&&" + isProtonCut;
  TString vert_study_mc_n = vert_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString vertVar_string = "BBtr_vz";
  cutvar vertVar_data(vertVar_string, vert_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar vertVar_mc_p(vertVar_string, vert_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar vertVar_mc_n(vertVar_string, vert_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  
  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* vert_dxhist_data = vertVar_data.get2DdxCutHisto();
  TH2D* vert_dxhist_mc_p = vertVar_mc_p.get2DdxCutHisto();
  TH2D* vert_dxhist_mc_n = vertVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_vertVar_data = vertVar_data.sliceAndProjectHisto_xMinxMax(vert_dxhist_data,vertVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_vertVar_mc_p = vertVar_mc_p.sliceAndProjectHisto_xMinxMax(vert_dxhist_mc_p,vertVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_vertVar_mc_n = vertVar_mc_n.sliceAndProjectHisto_xMinxMax(vert_dxhist_mc_n,vertVar_mc_n.getAxisTitle().Data(),"dx");

  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis vertVar_stability(vertVar_data,vertVar_mc_p,vertVar_mc_n,slice_vertVar_data,slice_vertVar_mc_p,slice_vertVar_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //W2 Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data W2
  TString W2_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC W2
  TString W2_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString W2_study_mc_p = W2_study_mc + "&&" + isProtonCut;
  TString W2_study_mc_n = W2_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString W2Var_string = "W2";
  cutvar W2Var_data(W2Var_string, W2_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar W2Var_mc_p(W2Var_string, W2_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar W2Var_mc_n(W2Var_string, W2_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* W2_dxhist_data = W2Var_data.get2DdxCutHisto();
  TH2D* W2_dxhist_mc_p = W2Var_mc_p.get2DdxCutHisto();
  TH2D* W2_dxhist_mc_n = W2Var_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_W2Var_data = W2Var_data.sliceAndProjectHisto_xMinxMax(W2_dxhist_data,W2Var_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_W2Var_mc_p = W2Var_mc_p.sliceAndProjectHisto_xMinxMax(W2_dxhist_mc_p,W2Var_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_W2Var_mc_n = W2Var_mc_n.sliceAndProjectHisto_xMinxMax(W2_dxhist_mc_n,W2Var_mc_n.getAxisTitle().Data(),"dx");

  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis W2Var_stability(W2Var_data,W2Var_mc_p,W2Var_mc_n,slice_W2Var_data,slice_W2Var_mc_p,slice_W2Var_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //BB E over p Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data BB E over p
  TString BB_E_over_p_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC BB E over p
  TString BB_E_over_p_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString BB_E_over_p_study_mc_p = BB_E_over_p_study_mc + "&&" + isProtonCut;
  TString BB_E_over_p_study_mc_n = BB_E_over_p_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString BB_E_over_pVar_string = "BB_E_over_p";
  cutvar BB_E_over_pVar_data(BB_E_over_pVar_string, BB_E_over_p_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar BB_E_over_pVar_mc_p(BB_E_over_pVar_string, BB_E_over_p_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar BB_E_over_pVar_mc_n(BB_E_over_pVar_string, BB_E_over_p_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* BB_E_over_p_dxhist_data = BB_E_over_pVar_data.get2DdxCutHisto();
  TH2D* BB_E_over_p_dxhist_mc_p = BB_E_over_pVar_mc_p.get2DdxCutHisto();
  TH2D* BB_E_over_p_dxhist_mc_n = BB_E_over_pVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_BB_E_over_pVar_data = BB_E_over_pVar_data.sliceAndProjectHisto_xMinxMax(BB_E_over_p_dxhist_data,BB_E_over_pVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_BB_E_over_pVar_mc_p = BB_E_over_pVar_mc_p.sliceAndProjectHisto_xMinxMax(BB_E_over_p_dxhist_mc_p,BB_E_over_pVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_BB_E_over_pVar_mc_n = BB_E_over_pVar_mc_n.sliceAndProjectHisto_xMinxMax(BB_E_over_p_dxhist_mc_n,BB_E_over_pVar_mc_n.getAxisTitle().Data(),"dx");

  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis BB_E_over_pVar_stability(BB_E_over_pVar_data,BB_E_over_pVar_mc_p,BB_E_over_pVar_mc_n,slice_BB_E_over_pVar_data,slice_BB_E_over_pVar_mc_p,slice_BB_E_over_pVar_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //HCal Energy Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data HCal Energy
  TString ehcal_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + dyCut + "&&" + OpticsCut_x + "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC ehcal 
  TString ehcal_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + dyCut+ "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString ehcal_study_mc_p = ehcal_study_mc + "&&" + isProtonCut;
  TString ehcal_study_mc_n = ehcal_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString ehcalVar_string = "ehcal";
  cutvar ehcalVar_data(ehcalVar_string, ehcal_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar ehcalVar_mc_p(ehcalVar_string, ehcal_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar ehcalVar_mc_n(ehcalVar_string, ehcal_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* ehcal_dxhist_data = ehcalVar_data.get2DdxCutHisto();
  TH2D* ehcal_dxhist_mc_p = ehcalVar_mc_p.get2DdxCutHisto();
  TH2D* ehcal_dxhist_mc_n = ehcalVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_ehcalVar_data = ehcalVar_data.sliceAndProjectHisto_xMinxMax(ehcal_dxhist_data,ehcalVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_ehcalVar_mc_p = ehcalVar_mc_p.sliceAndProjectHisto_xMinxMax(ehcal_dxhist_mc_p,ehcalVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_ehcalVar_mc_n = ehcalVar_mc_n.sliceAndProjectHisto_xMinxMax(ehcal_dxhist_mc_n,ehcalVar_mc_n.getAxisTitle().Data(),"dx");

  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis ehcalVar_stability(ehcalVar_data,ehcalVar_mc_p,ehcalVar_mc_n,slice_ehcalVar_data,slice_ehcalVar_mc_p,slice_ehcalVar_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //HCal Shower analog time diff Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data HCal Shower analog time diff
  TString hcal_sh_atime_diff_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&"  + dyCut+ "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  //For MC ps 
  TString hcal_sh_atime_diff_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString hcal_sh_atime_diff_study_mc_p = hcal_sh_atime_diff_study_mc + "&&" + isProtonCut;
  TString hcal_sh_atime_diff_study_mc_n = hcal_sh_atime_diff_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString hcal_sh_atime_diffVar_string = "hcal_sh_atime_diff";
  cutvar hcal_sh_atime_diffVar_data(hcal_sh_atime_diffVar_string, hcal_sh_atime_diff_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar hcal_sh_atime_diffVar_mc_p(hcal_sh_atime_diffVar_string, hcal_sh_atime_diff_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar hcal_sh_atime_diffVar_mc_n(hcal_sh_atime_diffVar_string, hcal_sh_atime_diff_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* hcal_sh_atime_diff_dxhist_data = hcal_sh_atime_diffVar_data.get2DdxCutHisto();
  TH2D* hcal_sh_atime_diff_dxhist_mc_p = hcal_sh_atime_diffVar_mc_p.get2DdxCutHisto();
  TH2D* hcal_sh_atime_diff_dxhist_mc_n = hcal_sh_atime_diffVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_hcal_sh_atime_diffVar_data = hcal_sh_atime_diffVar_data.sliceAndProjectHisto_xMinxMax(hcal_sh_atime_diff_dxhist_data,hcal_sh_atime_diffVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_hcal_sh_atime_diffVar_mc_p = hcal_sh_atime_diffVar_mc_p.sliceAndProjectHisto_xMinxMax(hcal_sh_atime_diff_dxhist_mc_p,hcal_sh_atime_diffVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_hcal_sh_atime_diffVar_mc_n = hcal_sh_atime_diffVar_mc_n.sliceAndProjectHisto_xMinxMax(hcal_sh_atime_diff_dxhist_mc_n,hcal_sh_atime_diffVar_mc_n.getAxisTitle().Data(),"dx");

  //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis hcal_sh_atime_diffVar_stability(hcal_sh_atime_diffVar_data,hcal_sh_atime_diffVar_mc_p,hcal_sh_atime_diffVar_mc_n,slice_hcal_sh_atime_diffVar_data,slice_hcal_sh_atime_diffVar_mc_p,slice_hcal_sh_atime_diffVar_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //dy Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data dy
  TString dy_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&"  + OpticsCut_x + "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC dy
  TString dy_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString dy_study_mc_p = dy_study_mc + "&&" + isProtonCut;
  TString dy_study_mc_n = dy_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString dyVar_string = "dy";
  cutvar dyVar_data(dyVar_string, dy_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar dyVar_mc_p(dyVar_string, dy_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar dyVar_mc_n(dyVar_string, dy_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* dy_dxhist_data = dyVar_data.get2DdxCutHisto();
  TH2D* dy_dxhist_mc_p = dyVar_mc_p.get2DdxCutHisto();
  TH2D* dy_dxhist_mc_n = dyVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_dyVar_data = dyVar_data.sliceAndProjectHisto_xMinxMax(dy_dxhist_data,dyVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_dyVar_mc_p = dyVar_mc_p.sliceAndProjectHisto_xMinxMax(dy_dxhist_mc_p,dyVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_dyVar_mc_n = dyVar_mc_n.sliceAndProjectHisto_xMinxMax(dy_dxhist_mc_n,dyVar_mc_n.getAxisTitle().Data(),"dx");

   //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis dyVar_stability(dyVar_data,dyVar_mc_p,dyVar_mc_n,slice_dyVar_data,slice_dyVar_mc_p,slice_dyVar_mc_n,FitType);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //Fid X, x expect Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data Fid X
  TString xexp_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&"  + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC Fid X
  TString xexp_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidYCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + dyCut+ "&&" + OpticsCut_y;
  TString xexp_study_mc_p = xexp_study_mc + "&&" + isProtonCut;
  TString xexp_study_mc_n = xexp_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString xexpVar_string = "xexp";
  cutvar xexpVar_data(xexpVar_string, xexp_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar xexpVar_mc_p(xexpVar_string, xexp_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar xexpVar_mc_n(xexpVar_string, xexp_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* xexp_dxhist_data = xexpVar_data.get2DdxCutHisto();
  TH2D* xexp_dxhist_mc_p = xexpVar_mc_p.get2DdxCutHisto();
  TH2D* xexp_dxhist_mc_n = xexpVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_xexpVar_data = xexpVar_data.sliceAndProjectHisto_xMinxMax(xexp_dxhist_data,xexpVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_xexpVar_mc_p = xexpVar_mc_p.sliceAndProjectHisto_xMinxMax(xexp_dxhist_mc_p,xexpVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_xexpVar_mc_n = xexpVar_mc_n.sliceAndProjectHisto_xMinxMax(xexp_dxhist_mc_n,xexpVar_mc_n.getAxisTitle().Data(),"dx");

   //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis xexpVar_stability(xexpVar_data,xexpVar_mc_p,xexpVar_mc_n,slice_xexpVar_data,slice_xexpVar_mc_p,slice_xexpVar_mc_n,FitType);

   ////////////////////////////////////////////////////////////////////////////////////////////
  //Fid Y, y expect Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data Fid Y
  TString yexp_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&"  + OpticsCut_x  + "&&" + HCal_Shower_atime_Cut;
  //For MC Fid Y
  TString yexp_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + dyCut+ "&&" + OpticsCut_x ;
  TString yexp_study_mc_p = yexp_study_mc + "&&" + isProtonCut;
  TString yexp_study_mc_n = yexp_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString yexpVar_string = "yexp";
  cutvar yexpVar_data(yexpVar_string, yexp_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar yexpVar_mc_p(yexpVar_string, yexp_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar yexpVar_mc_n(yexpVar_string, yexp_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* yexp_dxhist_data = yexpVar_data.get2DdxCutHisto();
  TH2D* yexp_dxhist_mc_p = yexpVar_mc_p.get2DdxCutHisto();
  TH2D* yexp_dxhist_mc_n = yexpVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_yexpVar_data = yexpVar_data.sliceAndProjectHisto_xMinxMax(yexp_dxhist_data,yexpVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_yexpVar_mc_p = yexpVar_mc_p.sliceAndProjectHisto_xMinxMax(yexp_dxhist_mc_p,yexpVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_yexpVar_mc_n = yexpVar_mc_n.sliceAndProjectHisto_xMinxMax(yexp_dxhist_mc_n,yexpVar_mc_n.getAxisTitle().Data(),"dx");

   //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis yexpVar_stability(yexpVar_data,yexpVar_mc_p,yexpVar_mc_n,slice_yexpVar_data,slice_yexpVar_mc_p,slice_yexpVar_mc_n,FitType);

   ////////////////////////////////////////////////////////////////////////////////////////////
  //Optics X validity Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data Optics X
  TString OpticsCut_x_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +   "&&"+ FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC OpticsCut_x
  TString OpticsCut_x_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  +  "&&" +FidYCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + dyCut+ "&&"  + OpticsCut_y;
  TString OpticsCut_x_study_mc_p = OpticsCut_x_study_mc + "&&" + isProtonCut;
  TString OpticsCut_x_study_mc_n = OpticsCut_x_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString OpticsCut_xVar_string = "BBtr_r_x-BBtr_r_th*0.9";
  cutvar OpticsCut_xVar_data(OpticsCut_xVar_string, OpticsCut_x_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar OpticsCut_xVar_mc_p(OpticsCut_xVar_string, OpticsCut_x_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar OpticsCut_xVar_mc_n(OpticsCut_xVar_string, OpticsCut_x_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* OpticsCut_x_dxhist_data = OpticsCut_xVar_data.get2DdxCutHisto();
  TH2D* OpticsCut_x_dxhist_mc_p = OpticsCut_xVar_mc_p.get2DdxCutHisto();
  TH2D* OpticsCut_x_dxhist_mc_n = OpticsCut_xVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_OpticsCut_xVar_data = OpticsCut_xVar_data.sliceAndProjectHisto_xMinxMax(OpticsCut_x_dxhist_data,OpticsCut_xVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_OpticsCut_xVar_mc_p = OpticsCut_xVar_mc_p.sliceAndProjectHisto_xMinxMax(OpticsCut_x_dxhist_mc_p,OpticsCut_xVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_OpticsCut_xVar_mc_n = OpticsCut_xVar_mc_n.sliceAndProjectHisto_xMinxMax(OpticsCut_x_dxhist_mc_n,OpticsCut_xVar_mc_n.getAxisTitle().Data(),"dx");

   //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis OpticsCut_xVar_stability(OpticsCut_xVar_data,OpticsCut_xVar_mc_p,OpticsCut_xVar_mc_n,slice_OpticsCut_xVar_data,slice_OpticsCut_xVar_mc_p,slice_OpticsCut_xVar_mc_n,FitType);

  //Check on correlation between W2 and Optics X
  TString W2_OpticsCut_x_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&"  + dyCut+ "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC W2 and dy study
  TString W2_OpticsCut_x_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&"+ FidYCut+ "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" +dyCut + "&&" + OpticsCut_y;
  TString W2_OpticsCut_x_study_mc_p = W2_OpticsCut_x_study_mc + "&&" + isProtonCut;
  TString W2_OpticsCut_x_study_mc_n = W2_OpticsCut_x_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString W2_OpticsCut_xVar_string = "BBtr_r_x-BBtr_r_th*0.9";
  cutvar W2_OpticsCut_xVar_data(W2_OpticsCut_xVar_string, W2_OpticsCut_x_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar W2_OpticsCut_xVar_mc_p(W2_OpticsCut_xVar_string, W2_OpticsCut_x_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar W2_OpticsCut_xVar_mc_n(W2_OpticsCut_xVar_string, W2_OpticsCut_x_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* W2_OpticsCut_x_W2hist_data = W2_OpticsCut_xVar_data.get2DW2CutHisto();
  TH2D* W2_OpticsCut_x_W2hist_mc_p = W2_OpticsCut_xVar_mc_p.get2DW2CutHisto();
  TH2D* W2_OpticsCut_x_W2hist_mc_n = W2_OpticsCut_xVar_mc_n.get2DW2CutHisto();


  //Optics Y validity Stability Analysis
  //First setup the cutvar string which should be every other cut but the one under consideration
  //For data Optics Y
  TString OpticsCut_y_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut +  "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut_x + "&&" + HCal_Shower_atime_Cut;
  //For MC OpticsCut_y
  TString OpticsCut_y_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + dyCut+ "&&"  + OpticsCut_x;
  TString OpticsCut_y_study_mc_p = OpticsCut_y_study_mc + "&&" + isProtonCut;
  TString OpticsCut_y_study_mc_n = OpticsCut_y_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString OpticsCut_yVar_string = "BBtr_r_y-0.9*BBtr_r_ph";
  cutvar OpticsCut_yVar_data(OpticsCut_yVar_string, OpticsCut_y_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar OpticsCut_yVar_mc_p(OpticsCut_yVar_string, OpticsCut_y_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar OpticsCut_yVar_mc_n(OpticsCut_yVar_string, OpticsCut_y_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* OpticsCut_y_dxhist_data = OpticsCut_yVar_data.get2DdxCutHisto();
  TH2D* OpticsCut_y_dxhist_mc_p = OpticsCut_yVar_mc_p.get2DdxCutHisto();
  TH2D* OpticsCut_y_dxhist_mc_n = OpticsCut_yVar_mc_n.get2DdxCutHisto();
  //Slice the histogram according to the given range of the cut var and make projections in based on the 2DHisto
  //For data  dx histo
  vector<TH1D*> slice_OpticsCut_yVar_data = OpticsCut_yVar_data.sliceAndProjectHisto_xMinxMax(OpticsCut_y_dxhist_data,OpticsCut_yVar_data.getAxisTitle().Data(),"dx");
  //For MC proton dx histo
  vector<TH1D*> slice_OpticsCut_yVar_mc_p = OpticsCut_yVar_mc_p.sliceAndProjectHisto_xMinxMax(OpticsCut_y_dxhist_mc_p,OpticsCut_yVar_mc_p.getAxisTitle().Data(),"dx");
  //For MC neutron dx histo
  vector<TH1D*> slice_OpticsCut_yVar_mc_n = OpticsCut_yVar_mc_n.sliceAndProjectHisto_xMinxMax(OpticsCut_y_dxhist_mc_n,OpticsCut_yVar_mc_n.getAxisTitle().Data(),"dx");

   //Setup a stability analysis object. This is useful for getting physics quantities for the overall dx slices
  stability_analysis OpticsCut_yVar_stability(OpticsCut_yVar_data,OpticsCut_yVar_mc_p,OpticsCut_yVar_mc_n,slice_OpticsCut_yVar_data,slice_OpticsCut_yVar_mc_p,slice_OpticsCut_yVar_mc_n,FitType);

  //Check on correlation between W2 and Optics Y
  TString W2_OpticsCut_y_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + FidXCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&"  + dyCut+ "&&" + OpticsCut_x + "&&" + HCal_Shower_atime_Cut;
  //For MC W2 and dy study
  TString W2_OpticsCut_y_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + FidXCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" +dyCut + "&&" + OpticsCut_x;
  TString W2_OpticsCut_y_study_mc_p = W2_OpticsCut_y_study_mc + "&&" + isProtonCut;
  TString W2_OpticsCut_y_study_mc_n = W2_OpticsCut_y_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString W2_OpticsCut_yVar_string = "BBtr_r_y-0.9*BBtr_r_ph";
  cutvar W2_OpticsCut_yVar_data(W2_OpticsCut_yVar_string, W2_OpticsCut_y_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar W2_OpticsCut_yVar_mc_p(W2_OpticsCut_yVar_string, W2_OpticsCut_y_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar W2_OpticsCut_yVar_mc_n(W2_OpticsCut_yVar_string, W2_OpticsCut_y_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* W2_OpticsCut_y_W2hist_data = W2_OpticsCut_yVar_data.get2DW2CutHisto();
  TH2D* W2_OpticsCut_y_W2hist_mc_p = W2_OpticsCut_yVar_mc_p.get2DW2CutHisto();
  TH2D* W2_OpticsCut_y_W2hist_mc_n = W2_OpticsCut_yVar_mc_n.get2DW2CutHisto();
  
  //Check on correlation between W2 and dy
  TString W2_dy_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&"  + OpticsCut_x + "&&" + OpticsCut_y + "&&" + HCal_Shower_atime_Cut;
  //For MC W2 and dy study
  TString W2_dy_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + FidXCut + "&&"+ FidYCut+ "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + OpticsCut_x + "&&" + OpticsCut_y;
  TString W2_dy_study_mc_p = W2_dy_study_mc + "&&" + isProtonCut;
  TString W2_dy_study_mc_n = W2_dy_study_mc + "&&" + isNeutronCut;
  //setup cutvar objects, central to this stability study, for data and mc
  TString W2_dyVar_string = "dy";
  cutvar W2_dyVar_data(W2_dyVar_string, W2_dy_study_data, data_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_data);
  //For MC proton
  cutvar W2_dyVar_mc_p(W2_dyVar_string, W2_dy_study_mc_p, mc_p_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);
  //For MC neutron
  cutvar W2_dyVar_mc_n(W2_dyVar_string, W2_dy_study_mc_n, mc_n_flag_string, slice_mode, left_right, kin, sbs_field, dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin,C_mc);

  //Get the 2D histgroms associated with the cut var for both data and mc
  TH2D* W2_dy_W2hist_data = W2_dyVar_data.get2DW2CutHisto();
  TH2D* W2_dy_W2hist_mc_p = W2_dyVar_mc_p.get2DW2CutHisto();
  TH2D* W2_dy_W2hist_mc_n = W2_dyVar_mc_n.get2DW2CutHisto();

  //draw the histos for pre shower plots
  TCanvas* c0 = psVar_stability.plotRsfTGraphError();
  TCanvas* c2 = psVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c3 = psVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c4 = psVar_stability.plotChi2_NDFGraph();
  TCanvas* c5 = psVar_stability.plotdxSliceHistos();
  TCanvas* c6 = psVar_stability.plotdxSliceResidual();
  TCanvas* c7 = psVar_stability.plot2DCutOverlay();
  TCanvas* c8 = psVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c9 = psVar_stability.plot1DCutOverlay();
  TCanvas* c10 = psVar_stability.plotNEntries();

  //draw thistos for the BBgem_nhits plots
  TCanvas* c11 = BBgem_nhitsVar_stability.plotRsfTGraphError();
  TCanvas* c12 = BBgem_nhitsVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c13 = BBgem_nhitsVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c14 = BBgem_nhitsVar_stability.plotChi2_NDFGraph();
  TCanvas* c15 = BBgem_nhitsVar_stability.plotdxSliceHistos();
  TCanvas* c16 = BBgem_nhitsVar_stability.plotdxSliceResidual();
  TCanvas* c17 = BBgem_nhitsVar_stability.plot2DCutOverlay();
  TCanvas* c18 = BBgem_nhitsVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c19 = BBgem_nhitsVar_stability.plot1DCutOverlay();
  TCanvas* c20 = BBgem_nhitsVar_stability.plotNEntries();

  //draw the histos for vertex plots
  TCanvas* c21 = BBgem_chi2ndfVar_stability.plotRsfTGraphError();
  TCanvas* c22 = BBgem_chi2ndfVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c23 = BBgem_chi2ndfVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c24 = BBgem_chi2ndfVar_stability.plotChi2_NDFGraph();
  TCanvas* c25 = BBgem_chi2ndfVar_stability.plotdxSliceHistos();
  TCanvas* c26 = BBgem_chi2ndfVar_stability.plotdxSliceResidual();
  TCanvas* c27 = BBgem_chi2ndfVar_stability.plot2DCutOverlay();
  TCanvas* c28 = BBgem_chi2ndfVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c29 = BBgem_chi2ndfVar_stability.plot1DCutOverlay();
  TCanvas* c30 = BBgem_chi2ndfVar_stability.plotNEntries();

  //draw the histos for vertex plots
  TCanvas* c31 = vertVar_stability.plotRsfTGraphError();
  TCanvas* c32 = vertVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c33 = vertVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c34 = vertVar_stability.plotChi2_NDFGraph();
  TCanvas* c35 = vertVar_stability.plotdxSliceHistos();
  TCanvas* c36 = vertVar_stability.plotdxSliceResidual();
  TCanvas* c37 = vertVar_stability.plot2DCutOverlay();
  TCanvas* c38 = vertVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c39 = vertVar_stability.plot1DCutOverlay();
  TCanvas* c40 = vertVar_stability.plotNEntries();
  
  //draw the histos for W2 plots
  TCanvas* c41 = W2Var_stability.plotRsfTGraphError();
  TCanvas* c42 = W2Var_stability.plotRsfTGraphError_xMin();
  TCanvas* c43 = W2Var_stability.plotRsfTGraphError_xMax();
  TCanvas* c44 = W2Var_stability.plotChi2_NDFGraph();
  TCanvas* c45 = W2Var_stability.plotdxSliceHistos();
  TCanvas* c46 = W2Var_stability.plotdxSliceResidual();
  TCanvas* c47 = W2Var_stability.plot2DCutOverlay();
  TCanvas* c48 = W2Var_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c49 = W2Var_stability.plot1DCutOverlay();
  TCanvas* c50 = W2Var_stability.plotNEntries();

  //draw the histos for E over p plots
  TCanvas* c51 = BB_E_over_pVar_stability.plotRsfTGraphError();
  TCanvas* c52 = BB_E_over_pVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c53 = BB_E_over_pVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c54 = BB_E_over_pVar_stability.plotChi2_NDFGraph();
  TCanvas* c55 = BB_E_over_pVar_stability.plotdxSliceHistos();
  TCanvas* c56 = BB_E_over_pVar_stability.plotdxSliceResidual();
  TCanvas* c57 = BB_E_over_pVar_stability.plot2DCutOverlay();
  TCanvas* c58 = BB_E_over_pVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c59 = BB_E_over_pVar_stability.plot1DCutOverlay();
  TCanvas* c60 = BB_E_over_pVar_stability.plotNEntries();

  //draw the histos for HCal E
  TCanvas* c61 = ehcalVar_stability.plotRsfTGraphError();
  TCanvas* c62 = ehcalVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c63 = ehcalVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c64 = ehcalVar_stability.plotChi2_NDFGraph();
  TCanvas* c65 = ehcalVar_stability.plotdxSliceHistos();
  TCanvas* c66 = ehcalVar_stability.plotdxSliceResidual();
  TCanvas* c67 = ehcalVar_stability.plot2DCutOverlay();
  TCanvas* c68 = ehcalVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c69 = ehcalVar_stability.plot1DCutOverlay();
  TCanvas* c70 = ehcalVar_stability.plotNEntries();

  //draw the histos for HCal sh analog time diff
  TCanvas* c71 = hcal_sh_atime_diffVar_stability.plotRsfTGraphError();
  TCanvas* c72 = hcal_sh_atime_diffVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c73 = hcal_sh_atime_diffVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c74 = hcal_sh_atime_diffVar_stability.plotChi2_NDFGraph();
  TCanvas* c75 = hcal_sh_atime_diffVar_stability.plotdxSliceHistos();
  TCanvas* c76 = hcal_sh_atime_diffVar_stability.plotdxSliceResidual();
  TCanvas* c77 = hcal_sh_atime_diffVar_stability.plot2DCutOverlay();
  TCanvas* c78 = hcal_sh_atime_diffVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c79 = hcal_sh_atime_diffVar_stability.plot1DCutOverlay();
  TCanvas* c80 = hcal_sh_atime_diffVar_stability.plotNEntries();

  //draw the histos for dy
  TCanvas* c81 = dyVar_stability.plotRsfTGraphError();
  TCanvas* c82 = dyVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c83 = dyVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c84 = dyVar_stability.plotChi2_NDFGraph();
  TCanvas* c85 = dyVar_stability.plotdxSliceHistos();
  TCanvas* c86 = dyVar_stability.plotdxSliceResidual();
  TCanvas* c87 = dyVar_stability.plot2DCutOverlay();
  TCanvas* c88 = dyVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c89 = dyVar_stability.plot1DCutOverlay();
  TCanvas* c90 = dyVar_stability.plotNEntries();

   //draw the histos for Fid X, x expect
  TCanvas* c91 = xexpVar_stability.plotRsfTGraphError();
  TCanvas* c92 = xexpVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c93 = xexpVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c94 = xexpVar_stability.plotChi2_NDFGraph();
  TCanvas* c95 = xexpVar_stability.plotdxSliceHistos();
  TCanvas* c96 = xexpVar_stability.plotdxSliceResidual();
  TCanvas* c97 = xexpVar_stability.plot2DCutOverlay();
  TCanvas* c98 = xexpVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c99 = xexpVar_stability.plot1DCutOverlay();
  TCanvas* c100 = xexpVar_stability.plotNEntries();

  //draw the histos for Fid Y, Y expect
  TCanvas* c101 = yexpVar_stability.plotRsfTGraphError();
  TCanvas* c102 = yexpVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c103 = yexpVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c104 = yexpVar_stability.plotChi2_NDFGraph();
  TCanvas* c105 = yexpVar_stability.plotdxSliceHistos();
  TCanvas* c106 = yexpVar_stability.plotdxSliceResidual();
  TCanvas* c107 = yexpVar_stability.plot2DCutOverlay();
  TCanvas* c108 = yexpVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c109 = yexpVar_stability.plot1DCutOverlay();
  TCanvas* c110 = yexpVar_stability.plotNEntries();

  //draw the histos for the Optics X cut
  TCanvas* c111 = OpticsCut_xVar_stability.plotRsfTGraphError();
  TCanvas* c112 = OpticsCut_xVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c113 = OpticsCut_xVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c114 = OpticsCut_xVar_stability.plotChi2_NDFGraph();
  TCanvas* c115 = OpticsCut_xVar_stability.plotdxSliceHistos();
  TCanvas* c116 = OpticsCut_xVar_stability.plotdxSliceResidual();
  TCanvas* c117 = OpticsCut_xVar_stability.plot2DCutOverlay();
  TCanvas* c118 = OpticsCut_xVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c119 = OpticsCut_xVar_stability.plot1DCutOverlay();
  TCanvas* c120 = OpticsCut_xVar_stability.plotNEntries();

  //draw the histos for the Optics Y cut
  TCanvas* c121 = OpticsCut_yVar_stability.plotRsfTGraphError();
  TCanvas* c122 = OpticsCut_yVar_stability.plotRsfTGraphError_xMin();
  TCanvas* c123 = OpticsCut_yVar_stability.plotRsfTGraphError_xMax();
  TCanvas* c124 = OpticsCut_yVar_stability.plotChi2_NDFGraph();
  TCanvas* c125 = OpticsCut_yVar_stability.plotdxSliceHistos();
  TCanvas* c126 = OpticsCut_yVar_stability.plotdxSliceResidual();
  TCanvas* c127 = OpticsCut_yVar_stability.plot2DCutOverlay();
  TCanvas* c128 = OpticsCut_yVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c129 = OpticsCut_yVar_stability.plot1DCutOverlay();
  TCanvas* c130 = OpticsCut_yVar_stability.plotNEntries();

  //Draw histos for W2 and dy study
  TCanvas* c131 = utility::printTH2D(W2_dy_W2hist_data,W2_dyVar_data.getDataMCFlag());
  TCanvas* c132 = utility::printTH2D(W2_dy_W2hist_mc_p,W2_dyVar_mc_p.getDataMCFlag());
  TCanvas* c133 = utility::printTH2D(W2_dy_W2hist_mc_n,W2_dyVar_mc_n.getDataMCFlag());

  //Draw histos for W2 and OpticsCut_x study
  TCanvas* c134 = utility::printTH2D(W2_OpticsCut_x_W2hist_data,W2_OpticsCut_xVar_data.getDataMCFlag());
  TCanvas* c135 = utility::printTH2D(W2_OpticsCut_x_W2hist_mc_p,W2_OpticsCut_xVar_mc_p.getDataMCFlag());
  TCanvas* c136 = utility::printTH2D(W2_OpticsCut_x_W2hist_mc_n,W2_OpticsCut_xVar_mc_n.getDataMCFlag());

  //Draw histos for W2 and OpticsCut_y study
  TCanvas* c137 = utility::printTH2D(W2_OpticsCut_y_W2hist_data,W2_OpticsCut_yVar_data.getDataMCFlag());
  TCanvas* c138 = utility::printTH2D(W2_OpticsCut_y_W2hist_mc_p,W2_OpticsCut_yVar_mc_p.getDataMCFlag());
  TCanvas* c139 = utility::printTH2D(W2_OpticsCut_y_W2hist_mc_n,W2_OpticsCut_yVar_mc_n.getDataMCFlag());

  //Write stuff to a pdf
  TString plotname = outfile;
  plotname.ReplaceAll(".root",".pdf");
  TString start = Form("%s%s",plotname.Data(),"(");
  //middle is the same as the name
  TString end = Form("%s%s",plotname.Data(),")");

  c0->Print(start.Data(),"pdf");
  c2->Print(plotname.Data(),"pdf");
  c3->Print(plotname.Data(),"pdf");
  c4->Print(plotname.Data(),"pdf");
  c5->Print(plotname.Data(),"pdf");
  c6->Print(plotname.Data(),"pdf");
  c7->Print(plotname.Data(),"pdf");
  c8->Print(plotname.Data(),"pdf");
  c9->Print(plotname.Data(),"pdf");
  c10->Print(plotname.Data(),"pdf");
  c11->Print(plotname.Data(),"pdf");
  c12->Print(plotname.Data(),"pdf");
  c13->Print(plotname.Data(),"pdf");
  c14->Print(plotname.Data(),"pdf");
  c15->Print(plotname.Data(),"pdf");
  c16->Print(plotname.Data(),"pdf");
  c17->Print(plotname.Data(),"pdf");
  c18->Print(plotname.Data(),"pdf");
  c19->Print(plotname.Data(),"pdf");
  c20->Print(plotname.Data(),"pdf");
  c21->Print(plotname.Data(),"pdf");
  c22->Print(plotname.Data(),"pdf");
  c23->Print(plotname.Data(),"pdf");
  c24->Print(plotname.Data(),"pdf");
  c25->Print(plotname.Data(),"pdf");
  c26->Print(plotname.Data(),"pdf");
  c27->Print(plotname.Data(),"pdf");
  c28->Print(plotname.Data(),"pdf");
  c29->Print(plotname.Data(),"pdf");
  c30->Print(plotname.Data(),"pdf");
  c31->Print(plotname.Data(),"pdf");
  c32->Print(plotname.Data(),"pdf");
  c33->Print(plotname.Data(),"pdf");
  c34->Print(plotname.Data(),"pdf");
  c35->Print(plotname.Data(),"pdf");
  c36->Print(plotname.Data(),"pdf");
  c37->Print(plotname.Data(),"pdf");
  c38->Print(plotname.Data(),"pdf");
  c39->Print(plotname.Data(),"pdf");
  c40->Print(plotname.Data(),"pdf");
  c41->Print(plotname.Data(),"pdf");
  c42->Print(plotname.Data(),"pdf");
  c43->Print(plotname.Data(),"pdf");
  c44->Print(plotname.Data(),"pdf");
  c45->Print(plotname.Data(),"pdf");
  c46->Print(plotname.Data(),"pdf");
  c47->Print(plotname.Data(),"pdf");
  c48->Print(plotname.Data(),"pdf");
  c49->Print(plotname.Data(),"pdf");
  c50->Print(plotname.Data(),"pdf");
  c51->Print(plotname.Data(),"pdf");
  c52->Print(plotname.Data(),"pdf");
  c53->Print(plotname.Data(),"pdf");
  c54->Print(plotname.Data(),"pdf");
  c55->Print(plotname.Data(),"pdf");
  c56->Print(plotname.Data(),"pdf");
  c57->Print(plotname.Data(),"pdf");
  c58->Print(plotname.Data(),"pdf");
  c59->Print(plotname.Data(),"pdf");
  c60->Print(plotname.Data(),"pdf");
  c61->Print(plotname.Data(),"pdf");
  c62->Print(plotname.Data(),"pdf");
  c63->Print(plotname.Data(),"pdf");
  c64->Print(plotname.Data(),"pdf");
  c65->Print(plotname.Data(),"pdf");
  c66->Print(plotname.Data(),"pdf");
  c67->Print(plotname.Data(),"pdf");
  c68->Print(plotname.Data(),"pdf");
  c69->Print(plotname.Data(),"pdf");
  c70->Print(plotname.Data(),"pdf");
  c71->Print(plotname.Data(),"pdf");
  c72->Print(plotname.Data(),"pdf");
  c73->Print(plotname.Data(),"pdf");
  c74->Print(plotname.Data(),"pdf");
  c75->Print(plotname.Data(),"pdf");
  c76->Print(plotname.Data(),"pdf");
  c77->Print(plotname.Data(),"pdf");
  c78->Print(plotname.Data(),"pdf");
  c79->Print(plotname.Data(),"pdf");
  c80->Print(plotname.Data(),"pdf");
  c81->Print(plotname.Data(),"pdf");
  c82->Print(plotname.Data(),"pdf");
  c83->Print(plotname.Data(),"pdf");
  c84->Print(plotname.Data(),"pdf");
  c85->Print(plotname.Data(),"pdf");
  c86->Print(plotname.Data(),"pdf");
  c87->Print(plotname.Data(),"pdf");
  c88->Print(plotname.Data(),"pdf");
  c89->Print(plotname.Data(),"pdf");
  c90->Print(plotname.Data(),"pdf");
  c91->Print(plotname.Data(),"pdf");
  c92->Print(plotname.Data(),"pdf");
  c93->Print(plotname.Data(),"pdf");
  c94->Print(plotname.Data(),"pdf");
  c95->Print(plotname.Data(),"pdf");
  c96->Print(plotname.Data(),"pdf");
  c97->Print(plotname.Data(),"pdf");
  c98->Print(plotname.Data(),"pdf");
  c99->Print(plotname.Data(),"pdf");
  c100->Print(plotname.Data(),"pdf");
  c101->Print(plotname.Data(),"pdf");
  c102->Print(plotname.Data(),"pdf");
  c103->Print(plotname.Data(),"pdf");
  c104->Print(plotname.Data(),"pdf");
  c105->Print(plotname.Data(),"pdf");
  c106->Print(plotname.Data(),"pdf");
  c107->Print(plotname.Data(),"pdf");
  c108->Print(plotname.Data(),"pdf");
  c109->Print(plotname.Data(),"pdf");
  c110->Print(plotname.Data(),"pdf");
  c111->Print(plotname.Data(),"pdf");
  c112->Print(plotname.Data(),"pdf");
  c113->Print(plotname.Data(),"pdf");
  c114->Print(plotname.Data(),"pdf");
  c115->Print(plotname.Data(),"pdf");
  c116->Print(plotname.Data(),"pdf");
  c117->Print(plotname.Data(),"pdf");
  c118->Print(plotname.Data(),"pdf");
  c119->Print(plotname.Data(),"pdf");
  c120->Print(plotname.Data(),"pdf");
  c121->Print(plotname.Data(),"pdf");
  c122->Print(plotname.Data(),"pdf");
  c123->Print(plotname.Data(),"pdf");
  c124->Print(plotname.Data(),"pdf");
  c125->Print(plotname.Data(),"pdf");
  c126->Print(plotname.Data(),"pdf");
  c127->Print(plotname.Data(),"pdf");
  c128->Print(plotname.Data(),"pdf");
  c129->Print(plotname.Data(),"pdf");
  c130->Print(plotname.Data(),"pdf");
  c131->Print(plotname.Data(),"pdf");
  c132->Print(plotname.Data(),"pdf");
  c133->Print(plotname.Data(),"pdf");
  c134->Print(plotname.Data(),"pdf");
  c135->Print(plotname.Data(),"pdf");
  c136->Print(plotname.Data(),"pdf");
  c137->Print(plotname.Data(),"pdf");
  c138->Print(plotname.Data(),"pdf");
  c139->Print(end.Data(),"pdf");
 
  fout->Write();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
