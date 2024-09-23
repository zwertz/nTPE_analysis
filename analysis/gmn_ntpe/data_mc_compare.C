
//Author Ezekiel Wertz
//04/09/2024
//Purpose of the script is to extract n/p cross-section ratio by comparing real with MC. This will be controlled by a particular kinematic and field setting. So it will not be implemented over multiple kinematics or field settings. The script will take in histograms and information from separate already parsed data and MC files. These files should already implement elastic cuts, an intime-clustering for HCal,and a fiducial cut. The goal is to then respectively fit the dx distributions to get a best fit. Then compose a sum of these fits allowing for scaling parameters which map to the total MC yields. Need to determine how best to handle the background, that is grounded in physics. Ultimately will apply the three floating parameter fit to the data and check residuals and chi-square. Should also figure out how to vary the fit function, ideally to yield better results. The ratio will then be extracted.

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

//handle the information sent to the report file
//Declare output file and open it
ofstream report;
report.open(reportfile);

//Print out the FF object stuff to the report file
 myFFs.print(report);

//make canvas to show data and MC compare plot
TCanvas* c0 = plots::plotDataMCFitsResiduals(dxhisto_allcuts,bg_shiftFit,"c0",true,report,myFFs);

//Do some fitting and get fit parameters. Do this for the background subtracted histogram
TH1D* hdx_data_nobg = plots::subtractBG(hdx_data_clone_bg,bg_shiftFit,totfit_ptr);
const char* fitType_nobg = "fitFullShiftNoBG";
fit_histogram *dxhisto_nobg = new fit_histogram(hdx_data_nobg,hdx_p,hdx_n,"shiftFit_nobg",fitType_nobg,4,9,hcalfit_low,hcalfit_high,fitopt.Data() );


//make canvas to show data and MC compare plot. For background subtracted version
TCanvas* c2 = plots::plotDataMCFitsResiduals_NoBG(dxhisto_nobg,"c2",true,report,myFFs);

TCanvas* c3 = plots::plotBGResiduals(hdx_data_bg,hdx_p_clone,hdx_n_clone,bg_shiftFit_clone,"c3","poly4 BG",fitType,shiftpar_vec,hcalfit_low,hcalfit_high,true);

//Canvas for no fid
TCanvas* c4 = plots::plotDataMCFitsResiduals(dxhisto_nofid,bg_shiftFit_nofid,"c4",true,report,myFFs);

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
