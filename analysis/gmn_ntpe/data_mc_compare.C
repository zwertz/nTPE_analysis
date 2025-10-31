
//Author Ezekiel Wertz
//04/09/2024
//The purose of the script is to extract the neutron to proton scale factor ratio, called Rsf, by comparing real data to MC simulated events for delta x distributions. The overall fit of the data delta x distribution takes into account the shape of the individual interpolated neutron and proton MC delta x distributions and a choice of background fit that is either a specific polynomial or a scaled interpolated background from a histogram. The neutron to proton cross-section ratio as modeled in the SIMC MC simulation is equivalently calculated using the calc_FFs_RCS_obj class. If one multiplies the extracted Rsf value and the MC neutron to proton cross-section ratio, then the experimental neutron to proton cross-section ratio is determined. This or these extracted values can then be used to determine Form Factors or ratios of FFs for the SBS GMn and nTPE experiments.This will be controlled by a particular kinematic and field setting. So it will not be implemented over multiple kinematics or field settings. The script will take in histograms and information from separate already parsed data and MC files. These files should already implement elastic cuts, an intime-clustering for HCal,and a fiducial cut. Currently this script considers multiple different shapes of the background and redetermines the value of Rsf for each different background shape. The Rsf values from the various different background shapes are then used to address the systematic effect due to choice of background shape.

//The exact ordering of this matters. ROOT for some reason cannot handle calls for files that have already been included.
#include "TFile.h"
#include <vector>
#include <utility>
#include "../../src/utility.C"
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/fits.C"
#include "../../src/parse_config.C"
#include "../../src/plots.C"
#include "../../src/calc_FFs_RCS_obj.C"
#include "../../src/fit_histogram.C"

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
TString kinematic_file = mainConfig.getKinFileName();

double hcalfit_low = exp_constants::hcalposXi_mc; //lower fit/bin limit for hcal dx plots. 
double hcalfit_high = exp_constants::hcalposXf_mc; //higher fit/bin limit for hcal dx plots.
double hcal_fitrange = exp_constants::hcal_vrange; //Full range of hcal dx plots

//Introduce an object that has information about the cross-section from the MC
calc_FFs_RCS_obj myFFs("SIMC",kinematic_file,kin);


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

TH1D *hdx_p = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_p")); //the MC dx histrogram for protons with all cuts

TH1D *hdx_n = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_n")); //the MC dx histrogram for neutrons with all cuts

TH1D *hdx_data_nofid = dynamic_cast<TH1D*>(data_file->Get("dx_cut_nofid")); //The data dx histogram with all cuts but fid

TH1D *hdx_p_nofid = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_nofid_p")); //the MC dx histrogram for protons with all cuts but fid

TH1D *hdx_n_nofid = dynamic_cast<TH1D*>(mc_file->Get("dx_cut_nofid_n")); //the MC dx histrogram for neutrons with all cuts but fid

TH1D *hdx_dt_bkgd = dynamic_cast<TH1D*>(data_file->Get("dx_coin_anticut")); //the data dx histogram for anti coin cut, a type of bkgd

TH1D *hdx_W2_anticut = dynamic_cast<TH1D*>(data_file->Get("dx_W2_anticut_morehigh")); //the data dx histogram to include some of the W2 inlastic background

//Histogram clones for this analysis
/*TH1D *hdx_data_clone = (TH1D*)(hdx_data->Clone("dx_data_clone"));

TH1D *hdx_data_clone_bg = (TH1D*)(hdx_data->Clone("dx_data_clone_bg"));
*/
TH1D *hdx_data_bg_gaussBG = (TH1D*)(hdx_data->Clone("dx_data_bg_gaussBG"));

TH1D *hdx_p_clone_gaussBG = (TH1D*)(hdx_p->Clone("dx_mc_p_clone_gaussBG"));

TH1D *hdx_n_clone_gaussBG = (TH1D*)(hdx_n->Clone("dx_mc_n_clone_gaussBG"));

TH1D *hdx_data_bg_poly2 = (TH1D*)(hdx_data->Clone("dx_data_bg_poly2"));

TH1D *hdx_p_clone_poly2 = (TH1D*)(hdx_p->Clone("dx_mc_p_clone_poly2"));

TH1D *hdx_n_clone_poly2 = (TH1D*)(hdx_n->Clone("dx_mc_n_clone_poly2"));

TH1D *hdx_data_bg_poly3 = (TH1D*)(hdx_data->Clone("dx_data_bg_poly3"));

TH1D *hdx_p_clone_poly3 = (TH1D*)(hdx_p->Clone("dx_mc_p_clone_poly3"));

TH1D *hdx_n_clone_poly3 = (TH1D*)(hdx_n->Clone("dx_mc_n_clone_poly3"));

TH1D *hdx_data_bg_dt_bkgd = (TH1D*)(hdx_data->Clone("dx_data_bg_dt_bkgd"));

TH1D *hdx_p_clone_dt_bkgd = (TH1D*)(hdx_p->Clone("dx_mc_p_clone_dt_bkgd"));

TH1D *hdx_n_clone_dt_bkgd = (TH1D*)(hdx_n->Clone("dx_mc_n_clone_dt_bkgd"));

TH1D *hdx_data_bg_W2_bkgd = (TH1D*)(hdx_data->Clone("dx_data_bg_W2_bkgd"));

TH1D *hdx_p_clone_W2_bkgd = (TH1D*)(hdx_p->Clone("dx_mc_p_clone_W2_bkgd"));

TH1D *hdx_n_clone_W2_bkgd = (TH1D*)(hdx_n->Clone("dx_mc_n_clone_W2_bkgd"));

/*
TH1D *hdx_p_nofid_clone = (TH1D*)(hdx_p_nofid->Clone("dx_mc_p_nofid_clone"));

TH1D *hdx_n_nofid_clone = (TH1D*)(hdx_n_nofid->Clone("dx_mc_n_nofid_clone"));

TH1D *hdx_data_fit_error = (TH1D*)(hdx_data->Clone("dx_data_fit_error"));

TH1D *hdx_data_nofid_clone = (TH1D*)(hdx_data_nofid->Clone("dx_data_nofid_clone"));
*/
TH1D *hdx_data_gaussBG_clone = (TH1D*)(hdx_data->Clone("dx_data_gaussBG_clone"));

TH1D *hdx_p_gaussBG_clone = (TH1D*)(hdx_p->Clone("dx_mc_p_gaussBG_clone"));

TH1D *hdx_n_gaussBG_clone = (TH1D*)(hdx_n->Clone("dx_mc_n_gaussBG_clone"));

TH1D *hdx_data_gaussBG_clone_bg = (TH1D*)(hdx_data->Clone("dx_data_gaussBG_clone_bg"));

TH1D *hdx_data_totfit_gaussBG = (TH1D*)(hdx_data->Clone("dx_data_totfit_gaussBG"));

TH1D *hdx_data_poly2_clone = (TH1D*)(hdx_data->Clone("dx_data_poly2_clone"));

TH1D *hdx_p_poly2_clone = (TH1D*)(hdx_p->Clone("dx_mc_p_poly2_clone"));

TH1D *hdx_n_poly2_clone = (TH1D*)(hdx_n->Clone("dx_mc_n_poly2_clone"));

TH1D *hdx_data_poly2_clone_bg = (TH1D*)(hdx_data->Clone("dx_data_poly2_clone_bg"));

TH1D *hdx_data_totfit_poly2 = (TH1D*)(hdx_data->Clone("dx_data_totfit_poly2"));

TH1D *hdx_data_poly3_clone = (TH1D*)(hdx_data->Clone("dx_data_poly3_clone"));

TH1D *hdx_p_poly3_clone = (TH1D*)(hdx_p->Clone("dx_mc_p_poly3_clone"));

TH1D *hdx_n_poly3_clone = (TH1D*)(hdx_n->Clone("dx_mc_n_poly3_clone"));

TH1D *hdx_data_poly3_clone_bg = (TH1D*)(hdx_data->Clone("dx_data_poly3_clone_bg"));

TH1D *hdx_data_totfit_poly3 = (TH1D*)(hdx_data->Clone("dx_data_totfit_poly3"));

TH1D *hdx_dt_bkgd_clone = (TH1D*)(hdx_dt_bkgd->Clone("hdx_dt_bkgd_clone"));

hdx_dt_bkgd_clone->Smooth();

TH1D *hdx_data_dt_bkgd_clone = (TH1D*)(hdx_data->Clone("dx_data_dt_bkgd_clone"));

TH1D *hdx_p_dt_bkgd_clone = (TH1D*)(hdx_p->Clone("dx_mc_p_dt_bkgd_clone"));

TH1D *hdx_n_dt_bkgd_clone = (TH1D*)(hdx_n->Clone("dx_mc_n_dt_bkgd_clone"));

TH1D *hdx_data_dt_bkgd_clone_bg = (TH1D*)(hdx_data->Clone("dx_data_dt_bkgd_clone_bg"));

TH1D *hdx_data_totfit_dt_bkgd = (TH1D*)(hdx_data->Clone("dx_data_totfit_dt_bkgd"));

TH1D *hdx_W2_bkgd_clone = (TH1D*)(hdx_W2_anticut->Clone("hdx_W2_bkgd_clone"));
hdx_W2_bkgd_clone->Smooth();

TH1D *hdx_data_W2_bkgd_clone = (TH1D*)(hdx_data->Clone("dx_data_W2_bkgd_clone"));

TH1D *hdx_p_W2_bkgd_clone = (TH1D*)(hdx_p->Clone("dx_mc_p_W2_bkgd_clone"));

TH1D *hdx_n_W2_bkgd_clone = (TH1D*)(hdx_n->Clone("dx_mc_n_W2_bkgd_clone"));

TH1D *hdx_data_W2_bkgd_clone_bg = (TH1D*)(hdx_data->Clone("dx_data_W2_bkgd_clone_bg"));

TH1D *hdx_data_totfit_W2_bkgd = (TH1D*)(hdx_data->Clone("dx_data_totfit_W2_bkgd"));

/*
//do some fitting and get the fit parameters, Error, ChiSquare, and NDF from the fit. This function is not perfect as it returns the parameter and Error pair vector and updates the pair to hold the ChiSquare and NDF. Better way to code this would be a structure that holds all of that information. 
const char* fitType = "fitFullShift_polyBG";
//Create a fit_histogram object which is a class that handles fitting dx histograms depending on input
fit_histogram *dxhisto_allcuts = new fit_histogram(hdx_data_clone,hdx_p,hdx_n,"shiftFit",fitType,4,9,hcalfit_low,hcalfit_high,fitopt.Data() ); 
vector<double> shiftpar_bkgd = dxhisto_allcuts->get_bkgd_params();
vector<double> shiftpar_bkgd_err = dxhisto_allcuts->get_bkgd_param_errs();
string fitTypeString = fitType;
vector<pair<double,double>> shiftpar_vec = dxhisto_allcuts->get_fitParamsErrs();

//fit parameters for dx plot without fid cut. A check essentially
const char* fitType_nofid = "fitFullShift_polyBG";
//Creat a fit_histogram object for the dx no fid
fit_histogram *dxhisto_nofid = new fit_histogram(hdx_data_nofid_clone,hdx_p_nofid,hdx_n_nofid,"shiftFit_nofid",fitType,4,9,hcalfit_low,hcalfit_high,fitopt.Data() );
vector<double> shiftpar_bkgd_nofid = dxhisto_nofid->get_bkgd_params();
vector<double> shiftpar_bkgd_nofid_err = dxhisto_nofid->get_bkgd_param_errs();

//Make background functions
TF1 *bg_shiftFit = new TF1("bg_shiftFit",dxhisto_allcuts,&fit_histogram::polyN_fit,hcalfit_low,hcalfit_high,5,"fit_histogram","polyN_fit");
//set the background parameters
	for(int j=0; j<5; ++j ){
	bg_shiftFit->SetParameter(j,shiftpar_bkgd[j]);
	bg_shiftFit->SetParError(j,shiftpar_bkgd_err[j]);
	}

//Make background function clone
TF1 *bg_shiftFit_clone = new TF1("bg_shiftFit_clone",dxhisto_allcuts,&fit_histogram::polyN_fit,hcalfit_low,hcalfit_high,5,"fit_histogram","polyN_fit");
//set the background parameters
 	for(int j=0; j<5; ++j ){
        bg_shiftFit_clone->SetParameter(j,shiftpar_bkgd[j]);
        bg_shiftFit_clone->SetParError(j,shiftpar_bkgd_err[j]);
        }

	
TF1 *total_fit_bg_error;
	//This conditional works but we will need to modify it if we need any new types of functions
        if(fitTypeString == "fitFull_polyBG"){
        total_fit_bg_error = new TF1("total_fit_bg_error",dxhisto_allcuts,&fit_histogram::fitFull_polyBG,hcalfit_low,hcalfit_high,9,"fit_histogram","fitFull_polyBG");
        }else if(fitTypeString =="fitFullShift_polyBG"){
        total_fit_bg_error = new TF1("total_fit_bg_error",dxhisto_allcuts,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,9,"fit_histogram","fitFullShift_polyBG");
        }else{
        cout << "The plot function you are trying to implement " << fitType << " is no good! Figure it out now!" << endl;
        total_fit_bg_error = new TF1("total_fit_bg_error",dxhisto_allcuts,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,9,"fit_histogram","fitFullShift_polyBG");
        }

TFitResultPtr totfit_ptr = hdx_data_fit_error->Fit(total_fit_bg_error,"QS");

//set the background parameters
	for(int k=0; k<shiftpar_vec.size(); ++k){
	total_fit_bg_error->SetParameter(k,shiftpar_vec[k].first);
	total_fit_bg_error->SetParError(k,shiftpar_vec[k].second);
	}

//background for the nofid plot
TF1 *bg_shiftFit_nofid = new TF1("bg_shiftFit_nofid",dxhisto_nofid,&fit_histogram::polyN_fit,hcalfit_low,hcalfit_high,5,"fit_histogram","polyN_fit");
	for(int j=0; j<5; ++j ){
        bg_shiftFit_nofid->SetParameter(j,shiftpar_bkgd_nofid[j]);
        bg_shiftFit_nofid->SetParError(j,shiftpar_bkgd_nofid_err[j]);
        }
*/
//*********************************************************
//Information for Systematic Study and/or final result
//********Gauss Background***
const char* fitType_gaussBG = "fitFullShift_gaussBG";
//Create a fit_histogram object which is a class that handles fitting dx histograms depending on input
fit_histogram *dxhisto_gaussBG = new fit_histogram(hdx_data_gaussBG_clone,hdx_p_gaussBG_clone,hdx_n_gaussBG_clone,"shiftFit_gaussBG",fitType_gaussBG,7,hcalfit_low,hcalfit_high,fitopt.Data() );
vector<double> gaussBG_param = dxhisto_gaussBG->get_bkgd_params();
vector<double> gaussBG_param_err = dxhisto_gaussBG->get_bkgd_param_errs();
string fitTypeString_gaussBG = fitType_gaussBG;
vector<pair<double,double>> gaussBG_vec = dxhisto_gaussBG->get_fitParamsErrs();

//Make background functions
TF1 *gaussBG_shiftFit = new TF1("gaussBG_shiftFit",dxhisto_gaussBG,&fit_histogram::Gauss,hcalfit_low,hcalfit_high,3,"fit_histogram","Gauss");

//set the background parameters
        for(int j=0; j<3; ++j ){
	gaussBG_shiftFit->SetParameter(j,gaussBG_param[j]);
        gaussBG_shiftFit->SetParError(j,gaussBG_param_err[j]);
        }

TF1 *gaussBG_shiftFit_clone = new TF1("gaussBG_shiftFit_clone",dxhisto_gaussBG,&fit_histogram::Gauss,hcalfit_low,hcalfit_high,3,"fit_histogram","Gauss");
//set the background parameters
        for(int j=0; j<3; ++j ){
        gaussBG_shiftFit_clone->SetParameter(j,gaussBG_param[j]);
        gaussBG_shiftFit_clone->SetParError(j,gaussBG_param_err[j]);
        }

TF1 *total_fit_gaussBG;
        //This conditional works but we will need to modify it if we need any new types of functions
        if(fitTypeString_gaussBG =="fitFullShift_gaussBG"){
	total_fit_gaussBG = new TF1("total_fit_gaussBG",dxhisto_gaussBG,&fit_histogram::fitFullShift_gaussBG,hcalfit_low,hcalfit_high,7,"fit_histogram","fitFullShift_gaussBG");
	}else{
        cout << "The plot function you are trying to implement " << fitType_gaussBG << " is no good! Figure it out now!" << endl;
        total_fit_gaussBG = new TF1("total_fit_gaussBG",dxhisto_gaussBG,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,7,"fit_histogram","fitFullShift_polyBG");
	}

TFitResultPtr totfit_ptr_gaussBG = hdx_data_totfit_gaussBG->Fit(total_fit_gaussBG,"QS");

//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg_gaussBG = plots::subtractBG(hdx_data_gaussBG_clone_bg,gaussBG_shiftFit,totfit_ptr_gaussBG);
const char* fitType_nobg_gaussBG = "fitFullShiftNoBG";
fit_histogram *dxhisto_gaussBG_nobg = new fit_histogram(hdx_data_nobg_gaussBG,hdx_p_gaussBG_clone,hdx_n_gaussBG_clone,"shiftFit_gaussBG no bg",fitType_nobg_gaussBG,7,hcalfit_low,hcalfit_high,fitopt.Data() );

/**********************************************************
//Information for Systematic Study and/or final result
*******Poly2 Background*********/
const char* fitType_poly2 = "fitFullShift_polyBG";
//Create a fit_histogram object which is a class that handles fitting dx histograms depending on input
fit_histogram *dxhisto_poly2 = new fit_histogram(hdx_data_poly2_clone,hdx_p_poly2_clone,hdx_n_poly2_clone,"fitFullShift_poly2",fitType_poly2,2,7,hcalfit_low,hcalfit_high,fitopt.Data() );
vector<double> poly2_param = dxhisto_poly2->get_bkgd_params();
vector<double> poly2_param_err = dxhisto_poly2->get_bkgd_param_errs();
string fitTypeString_poly2 = fitType_poly2;
vector<pair<double,double>> poly2_vec = dxhisto_poly2->get_fitParamsErrs();

//Make background functions
TF1 *poly2_shiftFit = new TF1("poly2_shiftFit",dxhisto_poly2,&fit_histogram::polyN_fit,hcalfit_low,hcalfit_high,3,"fit_histogram","polyN_fit");

//set the background parameters
        for(int j=0; j<3; ++j ){
        poly2_shiftFit->SetParameter(j,poly2_param[j]);
        poly2_shiftFit->SetParError(j,poly2_param_err[j]);
        }

TF1 *poly2_shiftFit_clone = new TF1("poly2_shiftFit_clone",dxhisto_poly2,&fit_histogram::polyN_fit,hcalfit_low,hcalfit_high,3,"fit_histogram","polyN_fit");
//set the background parameters
        for(int j=0; j<3; ++j ){
        poly2_shiftFit_clone->SetParameter(j,poly2_param[j]);
        poly2_shiftFit_clone->SetParError(j,poly2_param_err[j]);
        }

TF1 *total_fit_poly2;
        //This conditional works but we will need to modify it if we need any new types of functions
        if(fitTypeString_poly2 =="fitFullShift_polyBG"){
        total_fit_poly2 = new TF1("total_fit_poly2",dxhisto_poly2,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,7,"fit_histogram","fitFullShift_polyBG");
        }else{
        cout << "The plot function you are trying to implement " << fitType_poly2 << " is no good! Figure it out now!" << endl;
        total_fit_poly2 = new TF1("total_fit_poly2",dxhisto_poly2,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,7,"fit_histogram","fitFullShift_polyBG");
        }

TFitResultPtr totfit_ptr_poly2 = hdx_data_totfit_poly2->Fit(total_fit_poly2,"QS");

//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg_poly2 = plots::subtractBG(hdx_data_poly2_clone_bg,poly2_shiftFit,totfit_ptr_poly2);
const char* fitType_nobg_poly2 = "fitFullShiftNoBG";
fit_histogram *dxhisto_poly2_nobg = new fit_histogram(hdx_data_nobg_poly2,hdx_p_poly2_clone,hdx_n_poly2_clone,"shiftFit_poly2 no bg",fitType_nobg_poly2,2,7,hcalfit_low,hcalfit_high,fitopt.Data() );


/**********************************************************
//Information for Systematic Study and/or final result
*******Poly3 Background*********/
const char* fitType_poly3 = "fitFullShift_polyBG";
//Create a fit_histogram object which is a class that handles fitting dx histograms depending on input
fit_histogram *dxhisto_poly3 = new fit_histogram(hdx_data_poly3_clone,hdx_p_poly3_clone,hdx_n_poly3_clone,"fitFullShift_poly3",fitType_poly3,3,8,hcalfit_low,hcalfit_high,fitopt.Data() );
vector<double> poly3_param = dxhisto_poly3->get_bkgd_params();
vector<double> poly3_param_err = dxhisto_poly3->get_bkgd_param_errs();
string fitTypeString_poly3 = fitType_poly3;
vector<pair<double,double>> poly3_vec = dxhisto_poly3->get_fitParamsErrs();

//Make background functions
TF1 *poly3_shiftFit = new TF1("poly3_shiftFit",dxhisto_poly3,&fit_histogram::polyN_fit,hcalfit_low,hcalfit_high,4,"fit_histogram","polyN_fit");

//set the background parameters
        for(int j=0; j<4; ++j ){
        poly3_shiftFit->SetParameter(j,poly3_param[j]);
        poly3_shiftFit->SetParError(j,poly3_param_err[j]);
        }

TF1 *poly3_shiftFit_clone = new TF1("poly3_shiftFit_clone",dxhisto_poly3,&fit_histogram::polyN_fit,hcalfit_low,hcalfit_high,4,"fit_histogram","polyN_fit");
//set the background parameters
        for(int j=0; j<4; ++j ){
        poly3_shiftFit_clone->SetParameter(j,poly3_param[j]);
        poly3_shiftFit_clone->SetParError(j,poly3_param_err[j]);
        }

TF1 *total_fit_poly3;
        //This conditional works but we will need to modify it if we need any new types of functions
        if(fitTypeString_poly3 =="fitFullShift_polyBG"){
        total_fit_poly3 = new TF1("total_fit_poly3",dxhisto_poly3,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,8,"fit_histogram","fitFullShift_polyBG");
        }else{
        cout << "The plot function you are trying to implement " << fitType_poly3 << " is no good! Figure it out now!" << endl;
        total_fit_poly3 = new TF1("total_fit_poly3",dxhisto_poly3,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,8,"fit_histogram","fitFullShift_polyBG");
        }

TFitResultPtr totfit_ptr_poly3 = hdx_data_totfit_poly3->Fit(total_fit_poly3,"QS");

//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg_poly3 = plots::subtractBG(hdx_data_poly3_clone_bg,poly3_shiftFit,totfit_ptr_poly3);
const char* fitType_nobg_poly3 = "fitFullShiftNoBG";
fit_histogram *dxhisto_poly3_nobg = new fit_histogram(hdx_data_nobg_poly3,hdx_p_poly3_clone,hdx_n_poly3_clone,"shiftFit_poly3 no bg",fitType_nobg_poly3,3,8,hcalfit_low,hcalfit_high,fitopt.Data() );

/**********************************************************
//Information for Systematic Study and/or final result
*******Anti Coin Background*********/

const char* fitType_dt_bkgd = "fitFullShift_InterpolateBG";
//Create a fit_histogram object which is a class that handles fitting dx histograms depending on input
fit_histogram *dxhisto_dt_bkgd = new fit_histogram(hdx_data_dt_bkgd_clone,hdx_p_dt_bkgd_clone,hdx_n_dt_bkgd_clone,hdx_dt_bkgd_clone,"fitFullShift_dt_bkgd",fitType_dt_bkgd,5,hcalfit_low,hcalfit_high,fitopt.Data() );
vector<double> dt_bkgd_param = dxhisto_dt_bkgd->get_bkgd_params();
vector<double> dt_bkgd_param_err = dxhisto_dt_bkgd->get_bkgd_param_errs();
string fitTypeString_dt_bkgd = fitType_dt_bkgd;
vector<pair<double,double>> dt_bkgd_vec = dxhisto_dt_bkgd->get_fitParamsErrs();

//Make background functions
TF1 *dt_bkgd_shiftFit = new TF1("dt_bkgd_shiftFit",dxhisto_dt_bkgd,&fit_histogram::InterpolateBG,hcalfit_low,hcalfit_high,1,"fit_histogram","InterpolateBG");
//set the background parameters
        for(int j=0; j<1; ++j ){
        dt_bkgd_shiftFit->SetParameter(j,dt_bkgd_param[j]);
        dt_bkgd_shiftFit->SetParError(j,dt_bkgd_param_err[j]);
        }
TF1 *dt_bkgd_shiftFit_clone = new TF1("dt_bkgd_shiftFit_clone",dxhisto_dt_bkgd,&fit_histogram::InterpolateBG,hcalfit_low,hcalfit_high,1,"fit_histogram","InterpolateBG");
//set the background parameters
        for(int j=0; j<1; ++j ){
        dt_bkgd_shiftFit_clone->SetParameter(j,dt_bkgd_param[j]);
        dt_bkgd_shiftFit_clone->SetParError(j,dt_bkgd_param_err[j]);
        }

TF1 *total_fit_dt_bkgd;
        //This conditional works but we will need to modify it if we need any new types of functions
        if(fitTypeString_dt_bkgd == "fitFullShift_InterpolateBG"){
        total_fit_dt_bkgd = new TF1("total_fit_dt_bkgd",dxhisto_dt_bkgd,&fit_histogram::InterpolateBG,hcalfit_low,hcalfit_high,5,"fit_histogram","InterpolateBG");
        }else{
        cout << "The plot function you are trying to implement " << fitType_dt_bkgd << " is no good! Figure it out now!" << endl;
        total_fit_dt_bkgd = new TF1("total_fit_dt_bkgd",dxhisto_dt_bkgd,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,5,"fit_histogram","fitFullShift_polyBG");
        }

TFitResultPtr totfit_ptr_dt_bkgd = hdx_data_totfit_dt_bkgd->Fit(total_fit_dt_bkgd,"QS");

//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg_dt_bkgd = plots::subtractBG(hdx_data_dt_bkgd_clone_bg,dt_bkgd_shiftFit,totfit_ptr_dt_bkgd);
const char* fitType_nobg_dt_bkgd = "fitFullShiftNoBG";
fit_histogram *dxhisto_dt_bkgd_nobg = new fit_histogram(hdx_data_nobg_dt_bkgd,hdx_p_dt_bkgd_clone,hdx_n_dt_bkgd_clone,"shiftFit_dt_bkgd no bg",fitType_nobg_dt_bkgd,5,hcalfit_low,hcalfit_high,fitopt.Data() );

/**********************************************************
//Information for Systematic Study and/or final result
*******Anti W2 cut Background*********/
const char* fitType_W2_bkgd = "fitFullShift_InterpolateBG";
//Create a fit_histogram object which is a class that handles fitting dx histograms depending on input
fit_histogram *dxhisto_W2_bkgd = new fit_histogram(hdx_data_W2_bkgd_clone,hdx_p_W2_bkgd_clone,hdx_n_W2_bkgd_clone,hdx_W2_bkgd_clone,"fitFullShift_W2_bkgd",fitType_W2_bkgd,5,hcalfit_low,hcalfit_high,fitopt.Data() );
vector<double> W2_bkgd_param = dxhisto_W2_bkgd->get_bkgd_params();
vector<double> W2_bkgd_param_err = dxhisto_W2_bkgd->get_bkgd_param_errs();
string fitTypeString_W2_bkgd = fitType_W2_bkgd;
vector<pair<double,double>> W2_bkgd_vec = dxhisto_W2_bkgd->get_fitParamsErrs();

//Make background functions
TF1 *W2_bkgd_shiftFit = new TF1("W2_bkgd_shiftFit",dxhisto_W2_bkgd,&fit_histogram::InterpolateBG,hcalfit_low,hcalfit_high,1,"fit_histogram","InterpolateBG");
//set the background parameters
        for(int j=0; j<1; ++j ){
        W2_bkgd_shiftFit->SetParameter(j,W2_bkgd_param[j]);
        W2_bkgd_shiftFit->SetParError(j,W2_bkgd_param_err[j]);
        }
TF1 *W2_bkgd_shiftFit_clone = new TF1("W2_bkgd_shiftFit_clone",dxhisto_W2_bkgd,&fit_histogram::InterpolateBG,hcalfit_low,hcalfit_high,1,"fit_histogram","InterpolateBG");
//set the background parameters
        for(int j=0; j<1; ++j ){
        W2_bkgd_shiftFit_clone->SetParameter(j,W2_bkgd_param[j]);
        W2_bkgd_shiftFit_clone->SetParError(j,W2_bkgd_param_err[j]);
        }

TF1 *total_fit_W2_bkgd;
        //This conditional works but we will need to modify it if we need any new types of functions
        if(fitTypeString_W2_bkgd == "fitFullShift_InterpolateBG"){
        total_fit_W2_bkgd = new TF1("total_fit_W2_bkgd",dxhisto_W2_bkgd,&fit_histogram::InterpolateBG,hcalfit_low,hcalfit_high,5,"fit_histogram","InterpolateBG");
        }else{
        cout << "The plot function you are trying to implement " << fitType_W2_bkgd << " is no good! Figure it out now!" << endl;
        total_fit_W2_bkgd = new TF1("total_fit_W2_bkgd",dxhisto_W2_bkgd,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,5,"fit_histogram","fitFullShift_polyBG");
        }

TFitResultPtr totfit_ptr_W2_bkgd = hdx_data_totfit_W2_bkgd->Fit(total_fit_W2_bkgd,"QS");

//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg_W2_bkgd = plots::subtractBG(hdx_data_W2_bkgd_clone_bg,W2_bkgd_shiftFit,totfit_ptr_W2_bkgd);
const char* fitType_nobg_W2_bkgd = "fitFullShiftNoBG";
fit_histogram *dxhisto_W2_bkgd_nobg = new fit_histogram(hdx_data_nobg_W2_bkgd,hdx_p_W2_bkgd_clone,hdx_n_W2_bkgd_clone,"shiftFit_W2_bkgd no bg",fitType_nobg_W2_bkgd,5,hcalfit_low,hcalfit_high,fitopt.Data() );

//handle the information sent to the report file
//Declare output file and open it
ofstream report;
report.open(reportfile);

//Print out the FF object stuff to the report file
 myFFs.print(report);


//Old but worth not deleting
//make canvas to show data and MC compare plot
/*TCanvas* c0 = plots::plotDataMCFitsResiduals(dxhisto_allcuts,bg_shiftFit,"c0",true,report,myFFs);
//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg = plots::subtractBG(hdx_data_clone_bg,bg_shiftFit,totfit_ptr);
const char* fitType_nobg = "fitFullShiftNoBG";
fit_histogram *dxhisto_nobg = new fit_histogram(hdx_data_nobg,hdx_p,hdx_n,"shiftFit_nobg",fitType_nobg,4,9,hcalfit_low,hcalfit_high,fitopt.Data() );


//make canvas to show data and MC compare plot. For background subtracted version
TCanvas* c2 = plots::plotDataMCFitsResiduals_NoBG(dxhisto_nobg,"c2",true,report,myFFs);

TCanvas* c3 = plots::plotBGResiduals(hdx_data_bg,hdx_p_clone,hdx_n_clone,bg_shiftFit_clone,"c3","poly4 BG",fitType,shiftpar_vec,hcalfit_low,hcalfit_high,true);

//Canvas for no fid
TCanvas* c4 = plots::plotDataMCFitsResiduals(dxhisto_nofid,bg_shiftFit_nofid,"c4",true,report,myFFs);
*/

//Canvas for Gauss BG
TCanvas* c0 = plots::plotDataMCFitsResiduals(dxhisto_gaussBG,gaussBG_shiftFit,"c0",true,report,myFFs);
TCanvas* c2 = plots::plotDataMCFitsResiduals_NoBG(dxhisto_gaussBG_nobg,"c2",true,report,myFFs);
TCanvas* c3 = plots::plotBGResiduals(hdx_data_bg_gaussBG,hdx_p_clone_gaussBG,hdx_n_clone_gaussBG,gaussBG_shiftFit_clone,"c3","gauss BG",fitType_gaussBG,gaussBG_vec,hcalfit_low,hcalfit_high,true);

//Canvas for poly2 BG
TCanvas* c4 = plots::plotDataMCFitsResiduals(dxhisto_poly2,poly2_shiftFit,"c4",true,report,myFFs);
TCanvas* c5 = plots::plotDataMCFitsResiduals_NoBG(dxhisto_poly2_nobg,"c5",true,report,myFFs);
TCanvas* c6 = plots::plotBGResiduals(hdx_data_bg_poly2,hdx_p_clone_poly2,hdx_n_clone_poly2,poly2_shiftFit_clone,"c6","poly2 BG",fitType_poly2,poly2_vec,hcalfit_low,hcalfit_high,true);

//Canvas for poly3 BG
TCanvas* c7 = plots::plotDataMCFitsResiduals(dxhisto_poly3,poly3_shiftFit,"c7",true,report,myFFs);
TCanvas* c8 = plots::plotDataMCFitsResiduals_NoBG(dxhisto_poly3_nobg,"c8",true,report,myFFs);
TCanvas* c9 = plots::plotBGResiduals(hdx_data_bg_poly3,hdx_p_clone_poly3,hdx_n_clone_poly3,poly3_shiftFit_clone,"c9","poly3 BG",fitType_poly3,poly3_vec,hcalfit_low,hcalfit_high,true);

//Canvas for anti-coin BG
TCanvas* c10 = plots::plotDataMCFitsResiduals(dxhisto_dt_bkgd,dt_bkgd_shiftFit,"c10",true,report,myFFs);
TCanvas* c11 = plots::plotDataMCFitsResiduals_NoBG(dxhisto_dt_bkgd_nobg,"c11",true,report,myFFs);
TCanvas* c12 = plots::plotBGResiduals(hdx_data_bg_dt_bkgd,hdx_p_clone_dt_bkgd,hdx_n_clone_dt_bkgd,dt_bkgd_shiftFit_clone,"c12","dt BG",fitType_dt_bkgd,dt_bkgd_vec,hcalfit_low,hcalfit_high,true);

//Canvas for anti-W2 BG
TCanvas* c13 = plots::plotDataMCFitsResiduals(dxhisto_W2_bkgd,W2_bkgd_shiftFit,"c13",true,report,myFFs);
TCanvas* c14 = plots::plotDataMCFitsResiduals_NoBG(dxhisto_W2_bkgd_nobg,"c14",true,report,myFFs);
TCanvas* c15 = plots::plotBGResiduals(hdx_data_bg_W2_bkgd,hdx_p_clone_W2_bkgd,hdx_n_clone_W2_bkgd,W2_bkgd_shiftFit_clone,"c15","W2 BG",fitType_W2_bkgd,W2_bkgd_vec,hcalfit_low,hcalfit_high,true);

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
c15->Print(end.Data(),"pdf");


fout->Write();

// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
