//Author: Ezekiel Wertz
//10/29/24
//Purpose: This script currently relies on root files from scripts like HCal_p_Eff_uniformity_data.C or HCal_p_Eff_uniformity_MC.C. The goal is to aggregate multiple types of histograms from multiple magnetic field settings within a single kinematic to determine a pseudo-universal HCal efficiency map. This efficiency map will then be used to determine a relative efficiency correction factor to mimic the non-uniformity in HCal for the MC files like we see in the real data. A caveat to this analysis is that the binning and ranges of all the histograms need to be consisent. Similarly all the cuts applied to the histograms need to be consistent. If these conditions are not met, results may not be reasonable.
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

vector<TH2D*> xy_expect_num;
vector<TH2D*> xy_expect_denom;

//loop over the data files. Clone the histogram, put it in the vector for later
	for(int i=0; i<data_file_vec.size();i++){
	//Make copies of histograms from the corresponding files
	TH1D *h_x_expect_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hx_expect_hcalcut"));
	TH1D *h_x_expect_denon = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hx_expect_all"));

	TH1D *h_y_expect_num = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hy_expect_hcalcut"));
        TH1D *h_y_expect_denon = dynamic_cast<TH1D*> (data_file_vec[i]->Get("hy_expect_all"));

	TH2D *h_xy_expect_num = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hxy_expect_hcalcut"));
        TH2D *h_xy_expect_denon = dynamic_cast<TH2D*> (data_file_vec[i]->Get("hxy_expect_all"));


	//store them vecs
	x_expect_num.push_back(h_x_expect_num);
	x_expect_denom.push_back(h_x_expect_denon);

	y_expect_num.push_back(h_y_expect_num);
        y_expect_denom.push_back(h_y_expect_denon);

	xy_expect_num.push_back(h_xy_expect_num);
        xy_expect_denom.push_back(h_xy_expect_denon);

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

TH2D *h_xy_expect_num_total = new TH2D(*xy_expect_num[0]);
h_xy_expect_num_total->SetName("h_xy_expect_num_total");
h_xy_expect_num_total->SetTitle("h_xy_expect_num_total");
TH2D *h_xy_expect_denom_total = new TH2D(*xy_expect_denom[0]);
h_xy_expect_denom_total->SetName("h_xy_expect_denom_total");
h_xy_expect_denom_total->SetTitle("h_xy_expect_denom_total");


	//exclude the zeroth entry in the add since we started with that
	for(int j=1; j< data_file_vec.size(); j++){
	h_x_expect_num_total->Add(x_expect_num[j]);
	h_x_expect_denom_total->Add(x_expect_denom[j]);

	h_y_expect_num_total->Add(y_expect_num[j]);
        h_y_expect_denom_total->Add(y_expect_denom[j]);

	h_xy_expect_num_total->Add(xy_expect_num[j]);
        h_xy_expect_denom_total->Add(xy_expect_denom[j]);

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

  //Make the canvases which hold the fitted histogram
  TCanvas* c0 = plots::plotHCalEff(heff_vs_xexpect,"c0","heff_vs_xexpect_total","heff_vs_xexpect_total_wfit",fitx_low,fitx_high);

  TCanvas* c2 = plots::plotHCalEff(heff_vs_yexpect,"c2","heff_vs_yexpect_total","heff_vs_yexpect_total_wfit",fity_low,fity_high);

  TCanvas* c3 = plots::plot_Comp(h_x_expect_denom_total,h_x_expect_num_total,"c3","h_x_expect_denom_total","h_x_expect_num_total");

  TCanvas* c4 = plots::plot_Comp(h_y_expect_denom_total,h_y_expect_num_total,"c4","h_y_expect_denom_total","h_y_expect_num_total");

  TCanvas* c5 = plots::plot_HCalEffMap(heff_vs_xyexpect,"c5","heff_vs_xyexpect_total");

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
  c5->Print(end.Data(),"pdf");
  


//Write the info to file
fout->Write();

// Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
