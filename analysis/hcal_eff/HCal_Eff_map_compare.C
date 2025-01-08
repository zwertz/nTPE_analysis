//Author: Ezekiel Wertz
//12/03/2024
//Purpose: A script to compare HCal efficiency maps between two different files. Should be compatible with output from either HCal_Eff_map.C or HCal_p_Eff_uniformity_data/MC.C.Functionality should include comparing data to MC, but also across kinematics. Basic diagramtic comparision of both common 1-dimensional and 2-dimensional plots. 


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

void HCal_Eff_map_compare(const char *setup_file_name){
//Set this default to true so that way fits to histogram should be more correct. This effects statistical error
TH1::SetDefaultSumw2(kTRUE);

//Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//parse object to get in the information that The One Config file has and is manipulated
parse_config mainConfig(setup_file_name);

TString exp = mainConfig.getExp();
TString kin_1 = mainConfig.get_kin_1();
TString kin_2 = mainConfig.get_kin_2();
TString pass_1 = mainConfig.get_pass_1();
TString pass_2 = mainConfig.get_pass_2();
int sbs_field_1 = mainConfig.getSBS_field_1();
int sbs_field_2 = mainConfig.getSBS_field_2();
TString target_1 = mainConfig.get_targ_1();
TString target_2 = mainConfig.get_targ_2();
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

//setup input file info
TString Comp_file_1_name = mainConfig.get_Comp_file_1();
TString Comp_file_2_name = mainConfig.get_Comp_file_2();

vector<double> hcalpos;
vector<double> hcalaa;

if(pass_1 == "mc" || pass_2 == "mc"){
 //setup hcal physical bounds that match database
  hcalpos = cuts::hcal_Position_MC();

  //setup hcal active area with bounds that match database
  hcalaa = cuts::hcal_ActiveArea_MC(1,1);
}else{

//setup hcal physical bounds that match database for each pass
hcalpos = cuts::hcal_Position_data(pass_1);

//setup hcal active area with bounds that match database depending on pass
hcalaa = cuts::hcal_ActiveArea_data(1,1,pass_1);

}


//setup fiducial region based on dx and dy spot information
vector<double> hcalfid = cuts::hcalfid(dxsig_p,dxsig_n,dysig_p,hcalaa,dxsig_p_fac,dysig_p_fac);

TFile *Comp_file_1 = new TFile(Comp_file_1_name);
TFile *Comp_file_2 = new TFile(Comp_file_2_name);

//Setup output file
TString outfile = utility::makeOutputFileNameHCalEffMapComp(exp,pass_1,pass_2,kin_1,kin_2,sbs_field_1,sbs_field_2,target_1,target_2);
TFile *fout = new TFile(outfile,"RECREATE");

//conditional to determine histogram names because naming convention is not uniform
TString eff_xexpect_name_1;
TString eff_yexpect_name_1;
TString eff_xyexpect_name_1;

TString eff_xexpect_name_2;
TString eff_yexpect_name_2;
TString eff_xyexpect_name_2;
//Use the sbs_field parameters which correllate to the naming scheme
	if(sbs_field_1 < 0 && sbs_field_2 < 0 ){
	//Both are some form of collection of multiple magnetic field settings
	eff_xexpect_name_1 = "heff_vs_xexpect_total";
	eff_yexpect_name_1 = "heff_vs_yexpect_total";
	eff_xyexpect_name_1 = "heff_vs_xyexpect_total";

	eff_xexpect_name_2 = "heff_vs_xexpect_total";
        eff_yexpect_name_2 = "heff_vs_yexpect_total";
        eff_xyexpect_name_2 = "heff_vs_xyexpect_total";

	}else if(sbs_field_1 < 0 && sbs_field_2 >= 0){
	//sbs_field_2 is a single magnetic field setting
	eff_xexpect_name_1 = "heff_vs_xexpect_total";
        eff_yexpect_name_1 = "heff_vs_yexpect_total";
        eff_xyexpect_name_1 = "heff_vs_xyexpect_total";

	eff_xexpect_name_2 = "heff_vs_xexpect";
        eff_yexpect_name_2 = "heff_vs_yexpect";
        eff_xyexpect_name_2 = "heff_vs_xyexpect";

	}else if(sbs_field_1 >= 0 && sbs_field_2 < 0){
	//sbs_field_1 is a single magnetic field setting
	eff_xexpect_name_1 = "heff_vs_xexpect";
        eff_yexpect_name_1 = "heff_vs_yexpect";
        eff_xyexpect_name_1 = "heff_vs_xyexpect";

	eff_xexpect_name_2 = "heff_vs_xexpect_total";
        eff_yexpect_name_2 = "heff_vs_yexpect_total";
        eff_xyexpect_name_2 = "heff_vs_xyexpect_total";

	}else{
	//Both are some form of single magnetic field setting
	eff_xexpect_name_1 = "heff_vs_xexpect";
        eff_yexpect_name_1 = "heff_vs_yexpect";
        eff_xyexpect_name_1 = "heff_vs_xyexpect";

	eff_xexpect_name_2 = "heff_vs_xexpect";
        eff_yexpect_name_2 = "heff_vs_yexpect";
        eff_xyexpect_name_2 = "heff_vs_xyexpect";

	}




//Get histograms we will need for comparison
TH1D *heff_vs_xexpect_1 = dynamic_cast<TH1D*>(Comp_file_1->Get(eff_xexpect_name_1));
TH1D *heff_vs_yexpect_1 = dynamic_cast<TH1D*>(Comp_file_1->Get(eff_yexpect_name_1));
TH2D *heff_vs_xyexpect_1 = dynamic_cast<TH2D*>(Comp_file_1->Get(eff_xyexpect_name_1));

TH1D *heff_vs_xexpect_2 = dynamic_cast<TH1D*>(Comp_file_2->Get(eff_xexpect_name_2));
TH1D *heff_vs_yexpect_2 = dynamic_cast<TH1D*>(Comp_file_2->Get(eff_yexpect_name_2));
TH2D *heff_vs_xyexpect_2 = dynamic_cast<TH2D*>(Comp_file_2->Get(eff_xyexpect_name_2));

TString eff_xexpect_name_1_new = Form("%s_1",eff_xexpect_name_1.Data());
TString eff_yexpect_name_1_new = Form("%s_1",eff_yexpect_name_1.Data());
TString eff_xyexpect_name_1_new = Form("%s_1",eff_xyexpect_name_1.Data());

TString eff_xexpect_name_2_new = Form("%s_2",eff_xexpect_name_2.Data());
TString eff_yexpect_name_2_new = Form("%s_2",eff_yexpect_name_2.Data());
TString eff_xyexpect_name_2_new = Form("%s_2",eff_xyexpect_name_2.Data());



//Make clones of the histograms
TH1D *heff_vs_xexpect_1_clone = (TH1D*)(heff_vs_xexpect_1->Clone(eff_xexpect_name_1_new));
TH1D *heff_vs_yexpect_1_clone = (TH1D*)(heff_vs_yexpect_1->Clone(eff_yexpect_name_1_new));
TH2D *heff_vs_xyexpect_1_clone = (TH2D*)(heff_vs_xyexpect_1->Clone(eff_xyexpect_name_1_new));

TH1D *heff_vs_xexpect_2_clone = (TH1D*)(heff_vs_xexpect_2->Clone(eff_xexpect_name_2_new));
TH1D *heff_vs_yexpect_2_clone = (TH1D*)(heff_vs_yexpect_2->Clone(eff_yexpect_name_2_new));
TH2D *heff_vs_xyexpect_2_clone = (TH2D*)(heff_vs_xyexpect_2->Clone(eff_xyexpect_name_2_new));

//Make comparison histograms
TH1D *heff_vs_xexpect_comp = new TH1D( *heff_vs_xexpect_1_clone);
heff_vs_xexpect_comp->SetName(Form("%s_%s_compare",eff_xexpect_name_1_new.Data(),eff_xexpect_name_2_new.Data()));
heff_vs_xexpect_comp->SetTitle(Form("%s_%s_compare",eff_xexpect_name_1_new.Data(),eff_xexpect_name_2_new.Data()));
heff_vs_xexpect_comp->Divide(heff_vs_xexpect_1_clone,heff_vs_xexpect_2_clone);

TH1D *heff_vs_yexpect_comp = new TH1D(*heff_vs_yexpect_1_clone);
heff_vs_yexpect_comp->SetName(Form("%s_%s_compare",eff_yexpect_name_1_new.Data(),eff_yexpect_name_2_new.Data()));
heff_vs_yexpect_comp->SetTitle(Form("%s_%s_compare",eff_yexpect_name_1_new.Data(),eff_yexpect_name_2_new.Data()));
heff_vs_yexpect_comp->Divide(heff_vs_yexpect_1_clone,heff_vs_yexpect_2_clone);

TH2D *heff_vs_xyexpect_comp = new TH2D(*heff_vs_xyexpect_1_clone);
heff_vs_xyexpect_comp->SetName(Form("%s_%s_compare",eff_xyexpect_name_1_new.Data(),eff_xyexpect_name_2_new.Data()));
heff_vs_xyexpect_comp->SetTitle(Form("%s_%s_compare",eff_xyexpect_name_1_new.Data(),eff_xyexpect_name_2_new.Data()));
heff_vs_xyexpect_comp->Divide(heff_vs_xyexpect_1_clone,heff_vs_xyexpect_2_clone);

//make lines for physical HCal position
  vector<TLine*> Lines_pos = plots::setupLines(hcalpos,4,kBlack);

  //make lines for fiducial region
  vector<TLine*> Lines_Fid = plots::setupLines(hcalfid,4,kMagenta);


  //diff lines for the fiduical region on 1D histos
  TLine *LineL_FidX = plots::setupLine_Vert(0.0,1.9,fidx_min,2,kMagenta,2);
  TLine *LineR_FidX = plots::setupLine_Vert(0.0,1.9,fidx_max,2,kMagenta,2);
  TLine *LineL_FidY = plots::setupLine_Vert(0.0,1.9,fidy_min,2,kMagenta,2);
  TLine *LineR_FidY = plots::setupLine_Vert(0.0,1.9,fidy_max,2,kMagenta,2);

  vector<TLine*> Lines_Fid_diff;
  Lines_Fid_diff.push_back(LineL_FidX);
  Lines_Fid_diff.push_back(LineR_FidX);
  Lines_Fid_diff.push_back(LineL_FidY);
  Lines_Fid_diff.push_back(LineR_FidY);
  
  //diff lines for the fiduical region on 1D histos
  TLine *LineL_FidX_2 = plots::setupLine_Vert(0.0,1.0,fidx_min,2,kMagenta,2);
  TLine *LineR_FidX_2 = plots::setupLine_Vert(0.0,1.0,fidx_max,2,kMagenta,2);
  TLine *LineL_FidY_2 = plots::setupLine_Vert(0.0,1.0,fidy_min,2,kMagenta,2);
  TLine *LineR_FidY_2 = plots::setupLine_Vert(0.0,1.0,fidy_max,2,kMagenta,2);

  vector<TLine*> Lines_Fid_diff_2;
  Lines_Fid_diff_2.push_back(LineL_FidX_2);
  Lines_Fid_diff_2.push_back(LineR_FidX_2);
  Lines_Fid_diff_2.push_back(LineL_FidY_2);
  Lines_Fid_diff_2.push_back(LineR_FidY_2);



TString label_1 = Form("%s_%s_%s_%i",pass_1.Data(),kin_1.Data(),target_1.Data(),sbs_field_1);
TString label_2 = Form("%s_%s_%s_%i",pass_2.Data(),kin_2.Data(),target_2.Data(),sbs_field_2);


TCanvas* c0 = plots::plot_HCalEffMap_1D(heff_vs_xexpect_comp,"c0",Form("%s_%s_compare",eff_xexpect_name_1_new.Data(),eff_xexpect_name_2_new.Data()),Form("xexpect_%s_divide_%s",label_1.Data(),label_2.Data()),Lines_Fid_diff);
TCanvas* c2 = plots::plot_HCalEffMap_1D(heff_vs_yexpect_comp,"c2",Form("%s_%s_compare",eff_yexpect_name_1_new.Data(),eff_yexpect_name_2_new.Data()),Form("yexpect_%s_divide_%s",label_1.Data(),label_2.Data()),Lines_Fid_diff);

TCanvas* c3 = plots::plot_HCalEffMap_Comp(heff_vs_xyexpect_comp,"c3",Form("%s_%s_xycompare",label_1.Data(),label_2.Data()));
TCanvas* c4 = plots::plot_HCalEffMap_overlay_Comp(heff_vs_xyexpect_comp,"c4",Form("%s_%s_xycompare",label_1.Data(),label_2.Data()),Lines_pos,Lines_Fid);

TCanvas* c5 = plots::plot_Comp_1DEff(heff_vs_xexpect_1_clone,heff_vs_xexpect_2_clone,"c5",label_1,label_2,eff_xexpect_name_1_new,Lines_Fid_diff_2);
TCanvas* c6 = plots::plot_Comp_1DEff(heff_vs_yexpect_1_clone,heff_vs_yexpect_2_clone,"c6",label_1,label_2,eff_yexpect_name_1_new,Lines_Fid_diff_2);



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
 c6->Print(end.Data(),"pdf");


//Write the info to file
fout->Write();

// Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;


} //end main
