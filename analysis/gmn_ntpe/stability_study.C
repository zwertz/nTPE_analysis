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
  TString OpticsCut = mainConfig.getOpticsCut();
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
  TString ps_e_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC ps 
  TString ps_e_study_mc = TrackHitsCut + "&&" + TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut;
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
  TString BBgem_nhits_study_data = EnergyCut + "&&" + TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC BBGEM hits
  TString BBgem_nhits_study_mc = EnergyCut + "&&" +TrackQualityCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut;
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
  TString BBgem_chi2ndf_study_data = EnergyCut + "&&" + TrackHitsCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC BBGEM Track Chi2
  TString BBgem_chi2ndf_study_mc = EnergyCut + "&&" +TrackHitsCut + "&&" + TargetVertexCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut;
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
  TString vert_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC vertex 
  TString vert_study_mc = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + W2Cut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut;
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
  TString W2_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC W2
  TString W2_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + FidXCut + "&&" + FidYCut + "&&" + eOverpCut + "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut;
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
  TString BB_E_over_p_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" + HCal_Energy_Cut + "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC BB E over p
  TString BB_E_over_p_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut;
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
  TString ehcal_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + dyCut + "&&" + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC ehcal 
  TString ehcal_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + dyCut+ "&&" + OpticsCut;
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
  TString hcal_sh_atime_diff_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&"  + dyCut+ "&&" + OpticsCut;
  //For MC ps 
  TString hcal_sh_atime_diff_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&" + OpticsCut;
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
  TString dy_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&"  + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC dy
  TString dy_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" + FidYCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + OpticsCut;
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
  TString xexp_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidYCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&"  + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC Fid X
  TString xexp_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidYCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + dyCut+ "&&" + OpticsCut;
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
  TString yexp_study_data = TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut  + "&&" + W2Cut  +  "&&" + FidXCut + "&&" +eOverpCut+ "&&" + HCal_Energy_Cut + "&&" + dyCut+ "&&"  + OpticsCut + "&&" + HCal_Shower_atime_Cut;
  //For MC Fid Y
  TString yexp_study_mc =TrackHitsCut + "&&" + TrackQualityCut + "&&" + EnergyCut + "&&" + TargetVertexCut + "&&" + W2Cut  + "&&" + FidXCut + "&&" +eOverpCut+"&&" + HCal_Energy_Cut +  "&&" + dyCut+ "&&" + OpticsCut;
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

  //draw the histos for pre shower plots
  TCanvas* c0 = psVar_stability.plotRsfTGraphError();
  TCanvas* c2 = psVar_stability.plotdxSliceHistos();
  TCanvas* c3 = psVar_stability.plotdxSliceResidual();
  TCanvas* c4 = psVar_stability.plot2DCutOverlay();
  TCanvas* c5 = psVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c6 = psVar_stability.plot1DCutOverlay();
  TCanvas* c7 = psVar_stability.plotNEntries();

  //draw thistos for the BBgem_nhits plots
  TCanvas* c8 = BBgem_nhitsVar_stability.plotRsfTGraphError();
  TCanvas* c9 = BBgem_nhitsVar_stability.plotdxSliceHistos();
  TCanvas* c10 = BBgem_nhitsVar_stability.plotdxSliceResidual();
  TCanvas* c11 = BBgem_nhitsVar_stability.plot2DCutOverlay();
  TCanvas* c12 = BBgem_nhitsVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c13 = BBgem_nhitsVar_stability.plot1DCutOverlay();
  TCanvas* c14 = BBgem_nhitsVar_stability.plotNEntries();

  //draw the histos for vertex plots
  TCanvas* c15 = BBgem_chi2ndfVar_stability.plotRsfTGraphError();
  TCanvas* c16 = BBgem_chi2ndfVar_stability.plotdxSliceHistos();
  TCanvas* c17 = BBgem_chi2ndfVar_stability.plotdxSliceResidual();
  TCanvas* c18 = BBgem_chi2ndfVar_stability.plot2DCutOverlay();
  TCanvas* c19 = BBgem_chi2ndfVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c20 = BBgem_chi2ndfVar_stability.plot1DCutOverlay();
  TCanvas* c21 = BBgem_chi2ndfVar_stability.plotNEntries();

  //draw the histos for vertex plots
  TCanvas* c22 = vertVar_stability.plotRsfTGraphError();
  TCanvas* c23 = vertVar_stability.plotdxSliceHistos();
  TCanvas* c24 = vertVar_stability.plotdxSliceResidual();
  TCanvas* c25 = vertVar_stability.plot2DCutOverlay();
  TCanvas* c26 = vertVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c27 = vertVar_stability.plot1DCutOverlay();
  TCanvas* c28 = vertVar_stability.plotNEntries();
  
  //draw the histos for W2 plots
  TCanvas* c29 = W2Var_stability.plotRsfTGraphError();
  TCanvas* c30 = W2Var_stability.plotdxSliceHistos();
  TCanvas* c31 = W2Var_stability.plotdxSliceResidual();
  TCanvas* c32 = W2Var_stability.plot2DCutOverlay();
  TCanvas* c33 = W2Var_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c34 = W2Var_stability.plot1DCutOverlay();
  TCanvas* c35 = W2Var_stability.plotNEntries();

  //draw the histos for E over p plots
  TCanvas* c36 = BB_E_over_pVar_stability.plotRsfTGraphError();
  TCanvas* c37 = BB_E_over_pVar_stability.plotdxSliceHistos();
  TCanvas* c38 = BB_E_over_pVar_stability.plotdxSliceResidual();
  TCanvas* c39 = BB_E_over_pVar_stability.plot2DCutOverlay();
  TCanvas* c40 = BB_E_over_pVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c41 = BB_E_over_pVar_stability.plot1DCutOverlay();
  TCanvas* c42 = BB_E_over_pVar_stability.plotNEntries();

  //draw the histos for HCal E
  TCanvas* c43 = ehcalVar_stability.plotRsfTGraphError();
  TCanvas* c44 = ehcalVar_stability.plotdxSliceHistos();
  TCanvas* c45 = ehcalVar_stability.plotdxSliceResidual();
  TCanvas* c46 = ehcalVar_stability.plot2DCutOverlay();
  TCanvas* c47 = ehcalVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c48 = ehcalVar_stability.plot1DCutOverlay();
  TCanvas* c49 = ehcalVar_stability.plotNEntries();

  //draw the histos for HCal sh analog time diff
  TCanvas* c50 = hcal_sh_atime_diffVar_stability.plotRsfTGraphError();
  TCanvas* c51 = hcal_sh_atime_diffVar_stability.plotdxSliceHistos();
  TCanvas* c52 = hcal_sh_atime_diffVar_stability.plotdxSliceResidual();
  TCanvas* c53 = hcal_sh_atime_diffVar_stability.plot2DCutOverlay();
  TCanvas* c54 = hcal_sh_atime_diffVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c55 = hcal_sh_atime_diffVar_stability.plot1DCutOverlay();
  TCanvas* c56 = hcal_sh_atime_diffVar_stability.plotNEntries();

  //draw the histos for dy
  TCanvas* c57 = dyVar_stability.plotRsfTGraphError();
  TCanvas* c58 = dyVar_stability.plotdxSliceHistos();
  TCanvas* c59 = dyVar_stability.plotdxSliceResidual();
  TCanvas* c60 = dyVar_stability.plot2DCutOverlay();
  TCanvas* c61 = dyVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c62 = dyVar_stability.plot1DCutOverlay();
  TCanvas* c63 = dyVar_stability.plotNEntries();

   //draw the histos for Fid X, x expect
  TCanvas* c64 = xexpVar_stability.plotRsfTGraphError();
  TCanvas* c65 = xexpVar_stability.plotdxSliceHistos();
  TCanvas* c66 = xexpVar_stability.plotdxSliceResidual();
  TCanvas* c67 = xexpVar_stability.plot2DCutOverlay();
  TCanvas* c68 = xexpVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c69 = xexpVar_stability.plot1DCutOverlay();
  TCanvas* c70 = xexpVar_stability.plotNEntries();

  //draw the histos for Fid Y, Y expect
  TCanvas* c71 = yexpVar_stability.plotRsfTGraphError();
  TCanvas* c72 = yexpVar_stability.plotdxSliceHistos();
  TCanvas* c73 = yexpVar_stability.plotdxSliceResidual();
  TCanvas* c74 = yexpVar_stability.plot2DCutOverlay();
  TCanvas* c75 = yexpVar_stability.plot2DCutOverlay_CutRegion();
  TCanvas* c76 = yexpVar_stability.plot1DCutOverlay();
  TCanvas* c77 = yexpVar_stability.plotNEntries();

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
  c77->Print(end.Data(),"pdf");
 
  fout->Write();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
