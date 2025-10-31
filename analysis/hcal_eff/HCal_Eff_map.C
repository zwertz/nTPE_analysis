//Author: Ezekiel Wertz
//10/29/24
//Purpose: This script currently relies on root files from scripts like HCal_p_Eff_uniformity_data.C or HCal_p_Eff_uniformity_MC.C. The goal is to aggregate multiple types of histograms from multiple magnetic field settings within a single kinematic to determine a pseudo-universal HCal efficiency map. This efficiency map will then be used to determine a relative efficiency correction factor to mimic the non-uniformity in HCal for the MC files like we see in the real data. A caveat to this analysis is that the binning and ranges of all the histograms need to be consisent. Similarly all the cuts applied to the histograms need to be consistent. If these conditions are not met, results may not be reasonable. This script is essentially the second step in the HCal non-uniformity analysis. Reference my thesis for more details.
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
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/fits.C"
#include "../../src/data_object.C"
#include "../../src/cuts.C"
#include "../../src/physics.C"
#include "../../src/parse_config.C"
#include "../../src/plots.C"
#include "../../src/calc_FFs_RCS_obj.C"
#include "../../src/fit_histogram.C"

void HCal_Eff_map(const char *setup_file_name){
TH1::SetDefaultSumw2(kTRUE);
//Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//parse object to get in the information that The One Config file has and is manipulated
parse_config mainConfig(setup_file_name);

//Repurpose the run num functionality for multiple magnetic field settings
vector<int> magSets = mainConfig.getRunNums();
TString exp = mainConfig.getExp();
TString kin = mainConfig.getKin();
TString pass = mainConfig.getPass();
TString target = mainConfig.getTarg();
double fitx_low = mainConfig.get_fitxlow();
double fitx_high = mainConfig.get_fitxhigh();
double fity_low = mainConfig.get_fitylow();
double fity_high = mainConfig.get_fityhigh();
double dxsig_n = mainConfig.get_dxsign();
double dysig_n = mainConfig.get_dysign();
double dxsig_p = mainConfig.get_dxsigp();
double dysig_p = mainConfig.get_dysigp();
double dxsig_n_fac = mainConfig.get_dxSignFac();
double dxsig_p_fac = mainConfig.get_dxSigpFac();
double dysig_n_fac = mainConfig.get_dySignFac();
double dysig_p_fac = mainConfig.get_dySigpFac();
double fidx_min = mainConfig.get_fidxmin();
double fidx_max = mainConfig.get_fidxmax();
double fidy_min = mainConfig.get_fidymin();
double fidy_max = mainConfig.get_fidymax();

TString spot_choice = mainConfig.get_spot_choice();

vector<double> hcalpos;
vector<double> hcalaa;

if(pass == "mc"){
 //setup hcal physical bounds that match database
  hcalpos = cuts::hcal_Position_MC();

  //setup hcal active area with bounds that match database
  hcalaa = cuts::hcal_ActiveArea_MC(1,1);
}else{

//setup hcal physical bounds that match database for each pass
hcalpos = cuts::hcal_Position_data(pass);

//setup hcal active area with bounds that match database depending on pass
hcalaa = cuts::hcal_ActiveArea_data(1,1,pass);

}


//setup fiducial region based on dx and dy spot information
vector<double> hcalfid = cuts::hcalfid(dxsig_p,dxsig_n,dysig_p,hcalaa,dxsig_p_fac,dysig_p_fac);

//vector to hold input file names
vector<TString> data_file_names;
vector<TFile*> data_file_vec;

//loop over the entries in magnetic field settings to create the input file names
	for(int k=0; k<magSets.size();k++){
	//The name of the input file is the same as the output file but for the HCal uniformity analysis
	TString input_file_temp = utility::makeOutputFileNameHCalEffUniformity(exp,pass,kin,magSets[k],target);
	//store them for later just in case
	data_file_names.push_back(input_file_temp);
	
	//Create a new TFile and store it
	TFile *data_file_temp = new TFile(input_file_temp.Data());
	data_file_vec.push_back(data_file_temp);
	}

//outfile setup. Also needs automated
TString outfile = utility::makeOutputFileName_HCalEffMap(exp,pass,kin,target);
TFile *fout = new TFile(outfile,"RECREATE");

//vectors for the diff sets of histograms we need
vector<TH1D*> x_expect_num;
vector<TH1D*> x_expect_denom;

vector<TH1D*> y_expect_num;
vector<TH1D*> y_expect_denom;

vector<TH1D*> W2_num;
vector<TH1D*> W2_denom;

vector<TH1D*> x_hcal_num;
vector<TH1D*> x_hcal_denom;

vector<TH1D*> y_hcal_num;
vector<TH1D*> y_hcal_denom;

vector<TH2D*> xy_expect_num;
vector<TH2D*> xy_expect_denom;
vector<TH2D*> xy_expect_fail;

vector<TH2D*> xy_hcal_num;
vector<TH2D*> xy_hcal_denom;

vector<TH2D*> rowcol_hcal_num;
vector<TH2D*> rowcol_hcal_denom;

vector<TH1D*> dx_num;
vector<TH1D*> dy_num;
vector<TH2D*> dxdy_num;

vector<TH1D*> dx_denom;
vector<TH1D*> dy_denom;
vector<TH2D*> dxdy_denom;

//loop over the data files. Clone the histogram, put it in the vector for later
	for(int i=0; i<data_file_vec.size();i++){
	//Make copies of histograms from the corresponding files
	TH1D *h_x_expect_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hx_expect_hcalcut"));
	TH1D *h_x_expect_denom = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hx_expect_all"));

	TH1D *h_y_expect_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hy_expect_hcalcut"));
        TH1D *h_y_expect_denom = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hy_expect_all"));

	TH1D *h_x_hcal_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hx_HCal_cut"));
        TH1D *h_x_hcal_denom = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hx_HCal_all"));

        TH1D *h_y_hcal_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hy_HCal_cut"));
        TH1D *h_y_hcal_denom = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hy_HCal_all"));
	
	TH1D *h_W2_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hW2_hcalcut"));
	TH1D *h_W2_denom = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hW2_globcut"));

	TH2D *h_xy_expect_num = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hxy_expect_hcalcut"));
        TH2D *h_xy_expect_denom = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hxy_expect_all"));
	TH2D *h_xy_expect_fail = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hxy_expect_anticut"));

	TH2D *h_xy_hcal_num = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hxy_HCal_cut"));
        TH2D *h_xy_hcal_denom = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hxy_HCal_all"));

	TH2D *h_rowcol_hcal_num = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hrowcolHCal_cut"));
        TH2D *h_rowcol_hcal_denom = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hrowcolHCal_all"));

	TH1D *h_dx_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hdx_cut"));
	TH1D *h_dy_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hdy_cut"));
	TH2D *h_dxdy_num = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hdxdy_cut"));

	TH1D *h_dx_denom = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hdx_all"));
        TH1D *h_dy_denom = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hdy_all"));
        TH2D *h_dxdy_denom = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hdxdy_all"));

	//store them vecs
	x_expect_num.push_back(h_x_expect_num);
	x_expect_denom.push_back(h_x_expect_denom);

	y_expect_num.push_back(h_y_expect_num);
        y_expect_denom.push_back(h_y_expect_denom);

	x_hcal_num.push_back(h_x_hcal_num);
        x_hcal_denom.push_back(h_x_hcal_denom);

        y_hcal_num.push_back(h_y_hcal_num);
        y_hcal_denom.push_back(h_y_hcal_denom);

	W2_num.push_back(h_W2_num);
	W2_denom.push_back(h_W2_denom);

	xy_expect_num.push_back(h_xy_expect_num);
        xy_expect_denom.push_back(h_xy_expect_denom);
	xy_expect_fail.push_back(h_xy_expect_fail);

	xy_hcal_num.push_back(h_xy_hcal_num);
        xy_hcal_denom.push_back(h_xy_hcal_denom);

	rowcol_hcal_num.push_back(h_rowcol_hcal_num);
        rowcol_hcal_denom.push_back(h_rowcol_hcal_denom);

	dx_num.push_back(h_dx_num);
	dy_num.push_back(h_dy_num);
	dxdy_num.push_back(h_dxdy_num);

	dx_denom.push_back(h_dx_denom);
        dy_denom.push_back(h_dy_denom);
        dxdy_denom.push_back(h_dxdy_denom);

	}

//Start with the zeroth entry. Make a histogram copy for each histogram and then add the histograms together
TH1D *h_x_expect_num_total = new TH1D(*x_expect_num[0]);
h_x_expect_num_total->SetName("h_x_expect_num_total");
h_x_expect_num_total->SetTitle("h_x_expect_num_total");
TH1D *h_x_expect_denom_total = new TH1D(*x_expect_denom[0]);
h_x_expect_denom_total->SetName("h_x_expect_denom_total");
h_x_expect_denom_total->SetTitle("h_x_expect_denom_total");

TH1D *h_y_expect_num_total = new TH1D(*y_expect_num[0]);
h_y_expect_num_total->SetName("h_y_expect_num_total");
h_y_expect_num_total->SetTitle("h_y_expect_num_total");
TH1D *h_y_expect_denom_total = new TH1D(*y_expect_denom[0]);
h_y_expect_denom_total->SetName("h_y_expect_denom_total");
h_y_expect_denom_total->SetTitle("h_y_expect_denom_total");

TH1D *h_x_hcal_num_total = new TH1D(*x_hcal_num[0]);
h_x_hcal_num_total->SetName("h_x_hcal_num_total");
h_x_hcal_num_total->SetTitle("h_x_hcal_num_total");
TH1D *h_x_hcal_denom_total = new TH1D(*x_hcal_denom[0]);
h_x_hcal_denom_total->SetName("h_x_hcal_denom_total");
h_x_hcal_denom_total->SetTitle("h_x_hcal_denom_total");

TH1D *h_y_hcal_num_total = new TH1D(*y_hcal_num[0]);
h_y_hcal_num_total->SetName("h_y_hcal_num_total");
h_y_hcal_num_total->SetTitle("h_y_hcal_num_total");
TH1D *h_y_hcal_denom_total = new TH1D(*y_hcal_denom[0]);
h_y_hcal_denom_total->SetName("h_y_hcal_denom_total");
h_y_hcal_denom_total->SetTitle("h_y_hcal_denom_total");


TH1D *h_W2_num_total = new TH1D(*W2_num[0]);
h_W2_num_total->SetName("h_W2_num_total");
h_W2_num_total->SetTitle("h_W2_num_total");
TH1D *h_W2_denom_total = new TH1D(*W2_denom[0]);
h_W2_denom_total->SetName("h_W2_denom_total");
h_W2_denom_total->SetTitle("h_W2_denom_total");

TH2D *h_xy_expect_num_total = new TH2D(*xy_expect_num[0]);
h_xy_expect_num_total->SetName("h_xy_expect_num_total");
h_xy_expect_num_total->SetTitle("h_xy_expect_num_total");
TH2D *h_xy_expect_denom_total = new TH2D(*xy_expect_denom[0]);
h_xy_expect_denom_total->SetName("h_xy_expect_denom_total");
h_xy_expect_denom_total->SetTitle("h_xy_expect_denom_total");
TH2D *h_xy_expect_fail_total = new TH2D(*xy_expect_fail[0]);
h_xy_expect_fail_total->SetName("h_xy_expect_fail_total");
h_xy_expect_fail_total->SetTitle("h_xy_expect_fail_total");

TH2D *h_xy_hcal_num_total = new TH2D(*xy_hcal_num[0]);
h_xy_hcal_num_total->SetName("h_xy_hcal_num_total");
h_xy_hcal_num_total->SetTitle("h_xy_hcal_num_total");
TH2D *h_xy_hcal_denom_total = new TH2D(*xy_hcal_denom[0]);
h_xy_hcal_denom_total->SetName("h_xy_hcal_denom_total");
h_xy_hcal_denom_total->SetTitle("h_xy_hcal_denom_total");

TH2D *h_rowcol_hcal_num_total = new TH2D(*rowcol_hcal_num[0]);
h_rowcol_hcal_num_total->SetName("h_rowcol_hcal_num_total");
h_rowcol_hcal_num_total->SetTitle("h_rowcol_hcal_num_total");
TH2D *h_rowcol_hcal_denom_total = new TH2D(*rowcol_hcal_denom[0]);
h_rowcol_hcal_denom_total->SetName("h_rowcol_hcal_denom_total");
h_rowcol_hcal_denom_total->SetTitle("h_rowcol_hcal_denom_total");

TH1D *h_dx_num_total = new TH1D(*dx_num[0]);
h_dx_num_total->SetName("h_dx_num_total");
h_dx_num_total->SetTitle("h_dx_num_total");
TH1D *h_dy_num_total = new TH1D(*dy_num[0]);
h_dy_num_total->SetName("h_dy_num_total");
h_dy_num_total->SetTitle("h_dy_num_total");
TH2D *h_dxdy_num_total = new TH2D(*dxdy_num[0]);
h_dxdy_num_total->SetName("h_dxdy_num_total");
h_dxdy_num_total->SetTitle("h_dxdy_num_total");

TH1D *h_dx_denom_total = new TH1D(*dx_denom[0]);
h_dx_denom_total->SetName("h_dx_denom_total");
h_dx_denom_total->SetTitle("h_dx_denom_total");
TH1D *h_dy_denom_total = new TH1D(*dy_denom[0]);
h_dy_denom_total->SetName("h_dy_denom_total");
h_dy_denom_total->SetTitle("h_dy_denom_total");
TH2D *h_dxdy_denom_total = new TH2D(*dxdy_denom[0]);
h_dxdy_denom_total->SetName("h_dxdy_denom_total");
h_dxdy_denom_total->SetTitle("h_dxdy_denom_total");


	//exclude the zeroth entry in the add since we started with that
	for(int j=1; j< data_file_vec.size(); j++){
	h_x_expect_num_total->Add(x_expect_num[j]);
	h_x_expect_denom_total->Add(x_expect_denom[j]);

	h_y_expect_num_total->Add(y_expect_num[j]);
        h_y_expect_denom_total->Add(y_expect_denom[j]);

	h_x_hcal_num_total->Add(x_hcal_num[j]);
        h_x_hcal_denom_total->Add(x_hcal_denom[j]);

        h_y_hcal_num_total->Add(y_hcal_num[j]);
        h_y_hcal_denom_total->Add(y_hcal_denom[j]);

	h_W2_num_total->Add(W2_num[j]);
	h_W2_denom_total->Add(W2_denom[j]);

	h_xy_expect_num_total->Add(xy_expect_num[j]);
        h_xy_expect_denom_total->Add(xy_expect_denom[j]);
	h_xy_expect_fail_total->Add(xy_expect_fail[j]);

	h_xy_hcal_num_total->Add(xy_hcal_num[j]);
        h_xy_hcal_denom_total->Add(xy_hcal_denom[j]);

	h_rowcol_hcal_num_total->Add(rowcol_hcal_num[j]);
        h_rowcol_hcal_denom_total->Add(rowcol_hcal_denom[j]);

	h_dx_num_total->Add(dx_num[j]);
	h_dy_num_total->Add(dy_num[j]);
	h_dxdy_num_total->Add(dxdy_num[j]);

	h_dx_denom_total->Add(dx_denom[j]);
        h_dy_denom_total->Add(dy_denom[j]);
        h_dxdy_denom_total->Add(dxdy_denom[j]);

	}
//Make histogram for position dependent efficiency in the x-direction
TH1D *heff_vs_xexpect = new TH1D( *h_x_expect_num_total );
heff_vs_xexpect->SetName( "heff_vs_xexpect_total" );
heff_vs_xexpect->SetTitle( "heff_vs_xexpect_total" );
heff_vs_xexpect->Divide( h_x_expect_num_total, h_x_expect_denom_total );
//calculate the corresponding bin by bin error for the x-direction efficiency
for( int i=1; i<=heff_vs_xexpect->GetNbinsX(); i++ ){
    double eff = heff_vs_xexpect->GetBinContent(i);
    double N = std::max(1.0,h_x_expect_denom_total->GetBinContent(i));
    double err = sqrt(eff*(1.0-eff)/N); //binomial error
    heff_vs_xexpect->SetBinError(i,err);
  }

//Make histogram for position dependent efficiency in the y-direction
  TH1D *heff_vs_yexpect = new TH1D( *h_y_expect_num_total);
  heff_vs_yexpect->SetName( "heff_vs_yexpect_total" );
  heff_vs_yexpect->SetTitle( "heff_vs_yexpect_total" );
  heff_vs_yexpect->Divide( h_y_expect_num_total, h_y_expect_denom_total );
  //calculate the corresponding bin by bin error for the y-direction efficiency
  for( int j=1; j<=heff_vs_yexpect->GetNbinsX(); j++ ){
    double eff = heff_vs_yexpect->GetBinContent(j);
    double N = std::max(1.0,h_y_expect_denom_total->GetBinContent(j));
    double err = sqrt(eff*(1.0-eff)/N); //binomial error
    heff_vs_yexpect->SetBinError(j,err);
  }

  TH2D *heff_vs_xyexpect = new TH2D( *h_xy_expect_num_total );
  heff_vs_xyexpect->SetName("heff_vs_xyexpect_total" );
  heff_vs_xyexpect->SetTitle("heff_vs_xyexpect_total" );
  heff_vs_xyexpect->Divide( h_xy_expect_num_total, h_xy_expect_denom_total );

  for( int i=1; i<=heff_vs_xyexpect->GetNbinsX(); i++ ){
    for( int j=1; j<=heff_vs_xyexpect->GetNbinsY(); j++ ){

      int bin = heff_vs_xyexpect->GetBin(i,j);

      double eff = heff_vs_xyexpect->GetBinContent(bin);
      double N = std::max(1.0,h_xy_expect_denom_total->GetBinContent(bin));
      double err = sqrt(eff*(1.0-eff)/N);
      heff_vs_xyexpect->SetBinError(bin,err);
    }
  }

  TH2D *heff_vs_rowcol = new TH2D( *h_rowcol_hcal_num_total );
  heff_vs_rowcol->SetName( "heff_vs_rowcol_total" );
  heff_vs_rowcol->Divide( h_rowcol_hcal_num_total, h_rowcol_hcal_denom_total );


  double xexpect_max = heff_vs_xexpect->GetMaximum();
  double yexpect_max = heff_vs_yexpect->GetMaximum();

  double fidx_val_max = 1.05*xexpect_max;
  double fidy_val_max = 1.05*yexpect_max;

  //diff lines for the fiduical region on 1D histos
  TLine *LineL_FidX = plots::setupLine_Vert(0.0,fidx_val_max,fidx_min,2,kMagenta,2);
  TLine *LineR_FidX = plots::setupLine_Vert(0.0,fidx_val_max,fidx_max,2,kMagenta,2);
  TLine *LineL_FidY = plots::setupLine_Vert(0.0,fidy_val_max,fidy_min,2,kMagenta,2);
  TLine *LineR_FidY = plots::setupLine_Vert(0.0,fidy_val_max,fidy_max,2,kMagenta,2);

  vector<TLine*> Lines_Fid_diff;
  Lines_Fid_diff.push_back(LineL_FidX);
  Lines_Fid_diff.push_back(LineR_FidX);
  Lines_Fid_diff.push_back(LineL_FidY);
  Lines_Fid_diff.push_back(LineR_FidY);

  //Make the canvases which hold the fitted histogram
  TCanvas* c0 = plots::plotHCalEff(heff_vs_xexpect,"c0","heff_vs_xexpect_total","heff_vs_xexpect_total_wfit",fitx_low,fitx_high,Lines_Fid_diff);

  TCanvas* c2 = plots::plotHCalEff(heff_vs_yexpect,"c2","heff_vs_yexpect_total","heff_vs_yexpect_total_wfit",fity_low,fity_high,Lines_Fid_diff);

  TCanvas* c3 = plots::plot_Comp(h_x_expect_denom_total,h_x_expect_num_total,"c3","h_x_expect_denom_total","h_x_expect_num_total");

  TCanvas* c4 = plots::plot_Comp(h_y_expect_denom_total,h_y_expect_num_total,"c4","h_y_expect_denom_total","h_y_expect_num_total");

  TCanvas* c5 = plots::plot_Comp(h_W2_denom_total,h_W2_num_total,"c5","h_W2_denom_total","h_W2_num_total");
  
  TCanvas* c6 = plots::plot_HCalEffMap(heff_vs_xyexpect,"c6","heff_vs_xyexpect_total",target,spot_choice);

  //make lines for physical HCal position
  vector<TLine*> Lines_pos = plots::setupLines(hcalpos,4,kBlack);

  //make lines for fiducial region
  vector<TLine*> Lines_Fid = plots::setupLines(hcalfid,4,kMagenta);

  TCanvas* c7 = plots::plot_HCalEffMap_overlay(heff_vs_xyexpect,"c7","heff_vs_xyexpect_total_overlay",Lines_pos,Lines_Fid);

  TCanvas* c8 = plots::plot_HCalEffMap(h_xy_expect_num_total,"c8","h_xy_expect_num_total",target,spot_choice);

  TCanvas* c9 = plots::plot_HCalEffMap(h_xy_expect_denom_total,"c9","h_xy_expect_denom_total",target,spot_choice);

  TCanvas* c10 = plots::plot_HCalEffMap(h_xy_expect_fail_total,"c10","h_xy_expect_fail_total",target,spot_choice);

  TCanvas* c11 = plots::plot_HCalEffMap(heff_vs_rowcol,"c11","heff_vs_rowcol_total",target,spot_choice);

  TCanvas* c12 = plots::plot_HCalEffMap(h_dxdy_num_total,"c12","h_dxdy_num_total",target,spot_choice);

  TCanvas* c13 = plots::plot_HCalEffMap(h_dxdy_denom_total,"c13","h_dxdy_denom_total",target,spot_choice);

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
  c13->Print(end.Data(),"pdf");
  


//Write the info to file
fout->Write();

// Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
