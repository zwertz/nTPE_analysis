
//Author Ezekiel Wertz
//04/09/2024
//Purpose of the script is to extract n/p yield ratio by comparing real with MC. This will be controlled by a particular kinematic and field setting. So it will not be implemented over multiple kinematics or field settings. The script will take in histograms and information from separate already parsed data and MC files. These files should already implement elastic cuts, an intime-clustering for HCal,and a fiducial cut. The goal is to then respectively fit the dx distributions to get a best fit. Then compose a sum of these fits allowing for scaling parameters which map to the total MC yields. Need to determine how best to handle the background, that is grounded in physics. Ultimately will apply the three floating parameter fit to the data and check residuals and chi-square. Should also figure out how to vary the fit function, ideally to yield better results. The ratio will then be extracted.

#include "TFile.h"
#include <vector>
#include <utility>
#include "../../src/utility.C"
#include "../../src/exp_constants.C"
#include "../../src/fits.C"
#include "../../src/parse_config.C"
#include "../../src/plots.C"


//Histograms needed for interpolation of elastic signal
TH1D* hdx_p;
TH1D* hdx_n;
TH1D* hdx_p_nofid;
TH1D* hdx_n_nofid;


//This fit function  must sit outside of main and has to be in this file so that way it has some knowledge of the histogram

//Total function for MC + background.
double fitFull(double *x, double *param){
double dx = x[0];
double proton_scale = param[0];
double neutron_scale = param[1];

double proton = proton_scale * hdx_p->Interpolate(dx);
double neutron = neutron_scale * hdx_n->Interpolate(dx);

double fitFull = proton + neutron + fits::poly4_fit(x, &param[2]);
return fitFull;
}

//Full fit of proton and neutron peaks. No background
double fitFullNoBG(double *x, double *param){
double dx = x[0];
double proton_scale = param[0];
double neutron_scale = param[1];

double proton = proton_scale * hdx_p->Interpolate(dx);
double neutron = neutron_scale * hdx_n->Interpolate(dx);

double fitFull = proton + neutron;
return fitFull;
}

//Total function for MC + background. Independently shifts the neutron and proton peaks
double fitFullShift(double *x, double *param){
//MC float parameters
double proton_scale = param[0];
double neutron_scale = param[1];

double dx_shift_p = param[2];
double dx_shift_n = param[3];

//Apply the shifts before any interpolation
double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);

//The total function is the proton and neutron peaks + a 4th order polynomial for background
double fitFullshift = proton + neutron + fits::poly4_fit(x, &param[4]);
return fitFullshift;
}

//Total function for MC + background. Independently shifts the neutron and proton peaks
double fitFullShiftNoFid(double *x, double *param){
//MC float parameters
double proton_scale = param[0];
double neutron_scale = param[1];

double dx_shift_p = param[2];
double dx_shift_n = param[3];

//Apply the shifts before any interpolation
double proton = proton_scale * hdx_p_nofid->Interpolate(x[0] - dx_shift_p);
double neutron = neutron_scale * hdx_n_nofid->Interpolate(x[0] - dx_shift_n);

//The total function is the proton and neutron peaks + a 4th order polynomial for background
double fitFullshift = proton + neutron + fits::poly4_fit(x, &param[4]);
return fitFullshift;
}

//Full fit of proton and neutron peaks. No background. Independently shifts the neutron and proton peaks
double fitFullShiftNoBG(double *x, double *param){
//MC float parameters
double proton_scale = param[0];
double neutron_scale = param[1];

double dx_shift_p = param[2];
double dx_shift_n = param[3];

//Apply the shifts before any interpolation
double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);

double fitFullshift = proton + neutron;
return fitFullshift;
}


//Main
void data_mc_compare(const char *setup_file_name){

//Set this default to true so that way fits to histogram should be more correct. This effects statistical error
TH1::SetDefaultSumw2(kTRUE);

//Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//parse object to get in the information that The One Config file has and is manipulated
parse_config mainConfig(setup_file_name);
//mainConfig.printDataYields();

//store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
TString exp = mainConfig.getExp();
TString kin = mainConfig.getKin();
TString pass = mainConfig.getPass();
int sbs_field = mainConfig.getSBSField();
TString target = mainConfig.getTarg();
TString fitopt = mainConfig.getFitOpt();

double hcalfit_low = exp_constants::hcalposXi_mc; //lower fit/bin limit for hcal dx plots. 
double hcalfit_high = exp_constants::hcalposXf_mc; //higher fit/bin limit for hcal dx plots.
double hcal_fitrange = exp_constants::hcal_vrange; //Full range of hcal dx plots


//setup input file info
TString Data_input_file_name = mainConfig.getDataFile();
TString MC_input_file_name = mainConfig.getMCFileName();

TFile *data_file = new TFile(Data_input_file_name.Data());
TFile *mc_file = new TFile(MC_input_file_name.Data());

//Setup output file
TString outfile = utility::makeOutputFileName_DataMCComp(exp,pass,kin,sbs_field,target);
TString reportfile = utility::makeReportFileName_DataMCComp(exp,pass,kin,sbs_field,target);
TFile *fout = new TFile(outfile,"RECREATE");

//Get Histograms we will need from the respective data or MC files
TH1D *hdx_data = dynamic_cast<TH1D*>(data_file->Get("dx_cut")); //The data dx histogram with all cuts

hdx_p = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_p")); //the MC dx histrogram for protons with all cuts

hdx_n = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_n")); //the MC dx histrogram for neutrons with all cuts

TH1D *hdx_data_nofid = dynamic_cast<TH1D*>(data_file->Get("dx_cut_nofid")); //The data dx histogram with all cuts but fid

hdx_p_nofid = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_nofid_p")); //the MC dx histrogram for protons with all cuts but fid

hdx_n_nofid = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_nofid_n")); //the MC dx histrogram for neutrons with all cuts but fid


//Histogram clones for this analysis
TH1D *hdx_data_clone = (TH1D*)(hdx_data->Clone("dx_data_clone"));

TH1D *hdx_data_clone_bg = (TH1D*)(hdx_data->Clone("dx_data_clone_bg"));

TH1D *hdx_data_bg = (TH1D*)(hdx_data->Clone("dx_data_bg"));

TH1D *hdx_p_clone = (TH1D*)(hdx_p->Clone("dx_mc_p_clone"));

TH1D *hdx_n_clone = (TH1D*)(hdx_n->Clone("dx_mc_n_clone"));

TH1D *hdx_p_nofid_clone = (TH1D*)(hdx_p_nofid->Clone("dx_mc_p_nofid_clone"));

TH1D *hdx_n_nofid_clone = (TH1D*)(hdx_n_nofid->Clone("dx_mc_n_nofid_clone"));

TH1D *hdx_data_fit_error = (TH1D*)(hdx_data->Clone("dx_data_fit_error"));

TH1D *hdx_data_nofid_clone = (TH1D*)(hdx_data_nofid->Clone("dx_data_nofid_clone"));

//do some fitting and get the fit parameters, Error, ChiSquare, and NDF from the fit. This function is not perfect as it returns the parameter and Error pair vector and updates the pair to hold the ChiSquare and NDF. Better way to code this would be a structure that holds all of that information. 
pair<double,double> shiftQual;
const char* fitType = "fitFullShift";
auto shiftpar_vec = fits::fitAndFineFit(hdx_data_clone,"shiftFit",fitType,9,hcalfit_low,hcalfit_high,shiftQual,fitopt.Data());

//fit parameters for dx plot without fid cut. A check essentially
pair<double,double> shiftQual_nofid;
const char* fitType_nofid = "fitFullShiftNoFid";
auto shiftpar_vec_nofid = fits::fitAndFineFit(hdx_data_nofid_clone,"shiftFit_nofid",fitType_nofid,9,hcalfit_low,hcalfit_high,shiftQual_nofid,fitopt.Data());

//Make background functions
TF1 *bg_shiftFit = new TF1("bg_shiftFit",fits::poly4_fit,hcalfit_low,hcalfit_high,5);
//set the background parameters
	for(int j=0; j<5; ++j ){
	bg_shiftFit->SetParameter(j,shiftpar_vec[j+4].first);
	bg_shiftFit->SetParError(j,shiftpar_vec[j+4].second);
	}
//Make background function clone
TF1 *bg_shiftFit_clone = new TF1("bg_shiftFit_clone",fits::poly4_fit,hcalfit_low,hcalfit_high,5);
//set the background parameters
 	for(int j=0; j<5; ++j ){
        bg_shiftFit_clone->SetParameter(j,shiftpar_vec[j+4].first);
        bg_shiftFit_clone->SetParError(j,shiftpar_vec[j+4].second);
        }

TF1 *total_fit_bg_error = new TF1("total_fit_bg_error",fitType,hcalfit_low,hcalfit_high,9);
TFitResultPtr totfit_ptr = hdx_data_fit_error->Fit(total_fit_bg_error,"QS");
//set the background parameters
	for(int k=0; k<9; ++k){
	total_fit_bg_error->SetParameter(k,shiftpar_vec[k].first);
	total_fit_bg_error->SetParError(k,shiftpar_vec[k].second);
	}

//background for the nofid plot
TF1 *bg_shiftFit_nofid = new TF1("bg_shiftFit_nofid",fits::poly4_fit,hcalfit_low,hcalfit_high,5);
	for(int j=0; j<5; ++j ){
        bg_shiftFit_nofid->SetParameter(j,shiftpar_vec_nofid[j+4].first);
        bg_shiftFit_nofid->SetParError(j,shiftpar_vec_nofid[j+4].second);
        }

//handle the information sent to the report file
//Declare output file and open it
ofstream report;
report.open(reportfile);

//make canvas to show data and MC compare plot
TCanvas* c0 = plots::plotDataMCFitsResiduals(hdx_data_clone,hdx_p_clone,hdx_n_clone,bg_shiftFit,"c0","shiftfit poly4 BG",fitType,shiftpar_vec,shiftQual,hcalfit_low,hcalfit_high,true,report);

//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg = plots::subtractBG(hdx_data_clone_bg,bg_shiftFit,totfit_ptr);
pair<double,double> shiftQual_nobg;
const char* fitType_nobg = "fitFullShiftNoBG";
auto shiftpar_nobg_vec = fits::fitAndFineFit(hdx_data_nobg,"shiftFitNoBG",fitType_nobg,4,hcalfit_low,hcalfit_high,shiftQual_nobg,fitopt.Data());

//make canvas to show data and MC compare plot. For background subtracted version
TCanvas* c2 = plots::plotDataMCFitsResiduals_NoBG(hdx_data_nobg,hdx_p_clone,hdx_n_clone,"c2","shiftfit BG subtracted",fitType_nobg,shiftpar_nobg_vec,shiftQual_nobg,hcalfit_low,hcalfit_high,true);

TCanvas* c3 = plots::plotBGResiduals(hdx_data_bg,hdx_p_clone,hdx_n_clone,bg_shiftFit_clone,"c3","poly4 BG",fitType,shiftpar_vec,shiftQual,hcalfit_low,hcalfit_high,true);

//Canvas for no fid
TCanvas* c4 = plots::plotDataMCFitsResiduals(hdx_data_nofid_clone,hdx_p_nofid_clone,hdx_n_nofid_clone,bg_shiftFit_nofid,"c4","shiftfit poly4 BG nofid",fitType_nofid,shiftpar_vec_nofid,shiftQual_nofid,hcalfit_low,hcalfit_high,true,report);

report.close();

//Write stuff to a pdf
TString plotname = outfile;
plotname.ReplaceAll(".root",".pdf");
TString start = Form("%s%s",plotname.Data(),"(");
//middle is the same as the name
TString end = Form("%s%s",plotname.Data(),")");  

c0->Print(start.Data(),"pdf");
c2->Print(plotname.Data(),"pdf");
c3->Print(plotname.Data(),"pdf");
c4->Print(end.Data(),"pdf");

fout->Write();

// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
