//author Ezekiel Wertz
//The script requires multiple output root files from sf_dxdy.C, for different MC scale field values. The goal of this script is to get the proton peak position from the respective delta x graphs and plot them on a single TGraph as a function of the scale field value. Then fit the TGraph and get the fit parameters to ultimately see the correlation between the delta x proton peak position and the scale field value. To calibrate the MC scale field variable to realistically generate future MC files which are most identical to the real data delta x distributions. It should also be noted that this script was created prior to the common framework of this analysis. If the script would be used in the future it would probably need to be updated to be compatible with this framework. Likely would not work as is.


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
#include <iostream>
#include <fstream>
#include <sstream>
#include "analysis_utility_functions.h"
#include "BlindFactor.h"
#include <vector>
#include "TGraphErrors.h"
#include "TStopwatch.h"

const int SBS_Field = 100;
const string Kin = "SBS8";
const string targ = "LD2";
//number of field settings we will be evaluating
const int course_num_files = 12;
const int fine_num_files = 11; 
//generic sigma to be use in fitting later of the Gaussian
const double gen_sig = 0.2;
//scale fields for course study
const double fs_course[course_num_files]={10,20,30,40,50,60,70,80,90,100,110,120};
const string fs_course_string[course_num_files]={"10","20","30","40","50","60","70","80","90","100","110","120"};
//scale fields for fine study
const double fs_fine_50p[fine_num_files]={37.5,42.5,45.0,46.0,46.2,46.3,46.4,46.5,47.5,50.0,55.0};
const double fs_fine_70p[fine_num_files]={55.7,60.7,63.2,64.2,64.3,64.4,64.5,64.6,65.6,68.1,73.1};
const double fs_fine_100p[fine_num_files]={82.3,85.8,89.3,90.3,90.8,91.3,91.8,92.3,93.3,96.8,100.3};
double fs_fine[fine_num_files];

const string fs_fine_string_50p[fine_num_files]={"37.5","42.5","45.0","46.0","46.2","46.3","46.4","46.5","47.5","50.0","55.0"};
const string fs_fine_string_70p[fine_num_files]={"55.7","60.7","63.2","64.2","64.3","64.4","64.5","64.6","65.6","68.1","73.1"};
const string fs_fine_string_100p[fine_num_files]={"82.3","85.8","89.3","90.3","90.8","91.3","91.8","92.3","93.3","96.8","100.3"};
string fs_fine_string[fine_num_files];

//helper function to fit a Gaussian twice to the distribution
vector<double> fitGaussian_GetParams(TH1D* hist, double sig, double low, double high){

vector<double> daParams (4); //stores amplitude, mean, sigma, and error on the mean
//cout << low << " " << high <<  endl;
//Find the bin numbers corresponding to the specified range
int bin_low = hist->FindBin(low);
int bin_high = hist->FindBin(high);
//cout << low << " " << high <<  " " << bin_low << " " << bin_high << endl;
//Find the max content value and bin;
double maxContent = 0;
int maxBin = -1;

//loop to search and determine these values

for(Int_t l=bin_low; l<=bin_high;++l){
double bin_con = hist->GetBinContent(l);
	if(bin_con>maxContent){
	maxContent=bin_con;
	maxBin=l;
	}
//cout << l << " " << maxBin << " " << maxContent << endl;
}// end loop

double xMax = hist->GetXaxis()->GetBinCenter(maxBin);

//define initial fit range
double fitMin = xMax-sig;
double fitMax = xMax+sig;

//Now fit the histogram in the range
TF1 *gausFit = new TF1("gausFit","gaus",fitMin,fitMax);
hist->Fit(gausFit,"QR");

//get the parameters
double amplitude = gausFit->GetParameter(0);
double mean = gausFit->GetParameter(1);
double sigma = gausFit->GetParameter(2);
double mean_error = gausFit->GetParError(1);

//clean up
delete gausFit;

//Fit the histogram with the new fit info
double fitMin_2 = mean - sigma; 
double fitMax_2 = mean + sigma;

TF1 *gausFit_2 = new TF1("gausFit_2","gaus",fitMin_2,fitMax_2);
gausFit_2->SetParameters(amplitude,mean,sigma);
hist->Fit(gausFit_2,"R");

//Store the new parameters in the vector 
daParams[0] = gausFit_2->GetParameter(0); //Amplitude
daParams[1] = gausFit_2->GetParameter(1); //Mean
daParams[2] = gausFit_2->GetParameter(2); //Sigma
daParams[3] = gausFit_2->GetParError(1); //Error on Mean

return daParams;
}//End helper function

//main function we will do everything here
void sf_graph(bool course = false){
// Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//logistical book keeping for file i/o
double x_min_p = -3;
double x_max_p = -0.2;

double x_min_n = -0.2;
double x_max_n = 0.2;

string SBS_Field_strg;
string course_string;
int num_files;

//Conditional for SBS field info
if(SBS_Field == 50){
SBS_Field_strg = "50";
}else if(SBS_Field == 70){
SBS_Field_strg = "70";
}else if(SBS_Field == 100){
SBS_Field_strg = "100";
}else{
//We should never get here under normal stuff
cout << "SBS Field " << SBS_Field << " Something went wrong, fix it!" << endl;
return;
}//end of conditional

//conditional to determine if course or fine
if(course){
num_files = course_num_files;
course_string="course";
}else{//fine, that's the only other option
num_files = fine_num_files;
course_string="fine";
	//populate the fs_fine array with whatever array we need to
	for(int i=0; i<num_files;i++){
		if(SBS_Field == 50){
		fs_fine[i]=fs_fine_50p[i];
		fs_fine_string[i]=fs_fine_string_50p[i];
		}else if(SBS_Field == 70){
		fs_fine[i]=fs_fine_70p[i];
                fs_fine_string[i]=fs_fine_string_70p[i];
		}else if(SBS_Field == 100){
		fs_fine[i]=fs_fine_100p[i];
                fs_fine_string[i]=fs_fine_string_100p[i];
		}else{
		//We should never get here under normal stuff
		cout << "SBS Field " << SBS_Field << " Something went wrong in the loop, fix it!" << endl;
		return;
		}		
	}//end for loop
}//end of conditional

//logistical book keeping for file i/o
string field_scale_string[num_files];
double field_scale[num_files];

for(int j=0; j<num_files; j++){
	//So if course is true
	if(course){
	field_scale_string[j]=fs_course_string[j];
	field_scale[j] = fs_course[j];
	}
	//In this case it should be binary so this is if course is false or is fine
	else{
	field_scale_string[j]=fs_fine_string[j];
	field_scale[j]=fs_fine[j];
	}
}//end for loop

//define input directory path
string inputfiledir = "/work/halla/sbs/ewertz/nTPE_analysis/yields_output/";

//define the output file name
string outputfilename = inputfiledir + Form("MC_Zeke_Gmn_nTPE_%s_%s_%s_sf%sgraph.root",Kin.c_str(),targ.c_str(),SBS_Field_strg.c_str(),course_string.c_str());


//create some vectors to hold the various string, files, histograms we care about
vector<string> input_file_names (num_files);
vector<TFile*> input_files (num_files);
vector<TH1D*> hdx (num_files);
vector<TH1D*> hdx_fid (num_files);
vector<double> graph_x (num_files);
vector<double> graph_dx (num_files);
vector<double> graph_y_p (num_files);
vector<double> graph_dy_p (num_files);
vector<double> graph_x_fid (num_files);
vector<double> graph_dx_fid (num_files);
vector<double> graph_y_fid_p (num_files);
vector<double> graph_dy_fid_p (num_files);

vector<double> graph_y_n (num_files);
vector<double> graph_dy_n (num_files);
vector<double> graph_y_fid_n (num_files);
vector<double> graph_dy_fid_n (num_files);

//loop to populate the vectors and get the information we care about
for(int k=0; k<num_files; k++){
//cout << k << endl; 
string inputfile_name = inputfiledir + Form("MC_Zeke_Gmn_nTPE_%s_%s_%s_sf%s.root",Kin.c_str(),targ.c_str(),SBS_Field_strg.c_str(),field_scale_string[k].c_str());
//cout << inputfile_name << endl;
TFile *sf_file = new TFile(inputfile_name.c_str(),"READ");
TH1D *hist_dx_p = (TH1D*)sf_file->Get("dx");
TH1D *hist_dx_fid_p = (TH1D*)sf_file->Get("dx_fidcut");

TH1D *hist_dx_n = (TH1D*)sf_file->Get("dx");
TH1D *hist_dx_fid_n = (TH1D*)sf_file->Get("dx_fidcut");

input_file_names[k]=inputfile_name;
input_files[k]=sf_file;
hdx[k]=hist_dx_p;
hdx_fid[k]=hist_dx_fid_p;

vector<double> sf_params_p = fitGaussian_GetParams(hist_dx_p,gen_sig,x_min_p,x_max_p);
vector<double> sf_params_fid_p = fitGaussian_GetParams(hist_dx_fid_p,gen_sig,x_min_p,x_max_p);

vector<double> sf_params_n = fitGaussian_GetParams(hist_dx_n,gen_sig,x_min_n,x_max_n);
vector<double> sf_params_fid_n = fitGaussian_GetParams(hist_dx_fid_n,gen_sig,x_min_n,x_max_n);

graph_x[k] = field_scale[k]; // Value fixed by MC
graph_dx[k] = 0; //No error, cause value is set
graph_y_p[k] = sf_params_p[1]; //mean
graph_y_n[k] = sf_params_n[1]; //mean

graph_x_fid[k] = field_scale[k]; // Value fixed by MC
graph_dx_fid[k] = 0; //No error, cause value is set
graph_y_fid_p[k] = sf_params_fid_p[1]; //mean
graph_y_fid_n[k] = sf_params_fid_n[1]; //mean

if(course){
graph_dy_p[k] = sf_params_p[2]; //sigma
graph_dy_fid_p[k] = sf_params_fid_p[2]; //sigma
graph_dy_n[k] = sf_params_n[2]; //sigma
graph_dy_fid_n[k] = sf_params_fid_n[2]; //sigma

}else{ //only fine
graph_dy_p[k] = sf_params_p[3]; //error on mean
graph_dy_fid_p[k] = sf_params_fid_p[3]; //error on mean
graph_dy_n[k] = sf_params_n[3]; //error on mean
graph_dy_fid_n[k] = sf_params_fid_n[3]; //error on mean

}

}//end for loop

//setup output file

TFile *fout = new TFile(outputfilename.c_str(),"RECREATE");

//convert the vectors to array because for some reason TGraph does not like them
double * graph_x_arr = graph_x.data();
double * graph_dx_arr = graph_dx.data();
double * graph_y_p_arr = graph_y_p.data();
double * graph_dy_p_arr = graph_dy_p.data();
double * graph_y_n_arr = graph_y_n.data();
double * graph_dy_n_arr = graph_dy_n.data();


double * graph_x_fid_arr = graph_x_fid.data();
double * graph_dx_fid_arr = graph_dx_fid.data();
double * graph_y_fid_p_arr = graph_y_fid_p.data();
double * graph_dy_fid_p_arr = graph_dy_fid_p.data();
double * graph_y_fid_n_arr = graph_y_fid_n.data();
double * graph_dy_fid_n_arr = graph_dy_fid_n.data();


//Create graph
gStyle->SetEndErrorSize(0);
TGraph *graph_p = new TGraphErrors(num_files,graph_x_arr,graph_y_p_arr,graph_dx_arr,graph_dy_p_arr);
graph_p->SetMarkerStyle(20); // e.g., 20 is a full circle
graph_p->SetMarkerSize(1.0); // Size of the marker
graph_p->SetTitle(Form("%s %s Proton Peak Means vs Scale Field",Kin.c_str(),targ.c_str()));

// Fit with a straight line
TF1 *fitFunc_p = new TF1("fitFunc", "pol1"); // pol1 for linear fit
graph_p->Fit(fitFunc_p, "Q"); // "Q" for quiet mode (no print)

// Print slope and y-intercept
double slope_p = fitFunc_p->GetParameter(1);
double yIntercept_p = fitFunc_p->GetParameter(0);

//set up the canvas
TCanvas *c0 = new TCanvas("c0","scale field comparisons",1600,1200);

graph_p->Draw("AP");

// Create and populate the legend
TLegend *legend_p = new TLegend(0.5, 0.7, 0.8, 0.9); // Adjust these coordinates as needed
TString slopeStr_p;
slopeStr_p.Form("Slope: %0.5f", slope_p);
TString interceptStr_p;
interceptStr_p.Form("Y-Intercept: %0.5f", yIntercept_p);
legend_p->AddEntry(graph_p, slopeStr_p, "l");
legend_p->AddEntry(graph_p, interceptStr_p, "l");
legend_p->Draw();

c0->Update();
c0->Write();

TGraph *graph_fid_p = new TGraphErrors(num_files,graph_x_fid_arr,graph_y_fid_p_arr,graph_dx_fid_arr,graph_dy_fid_p_arr);
graph_fid_p->SetMarkerStyle(20); // e.g., 20 is a full circle
graph_fid_p->SetMarkerSize(1.0); // Size of the marker
graph_fid_p->SetTitle(Form("%s %s Proton Peak Means vs Scale Field, fidcut",Kin.c_str(),targ.c_str()));

// Fit with a straight line
TF1 *fitFunc_fid_p = new TF1("fitFunc_fid_p", "pol1"); // pol1 for linear fit
graph_fid_p->Fit(fitFunc_fid_p, "Q"); // "Q" for quiet mode (no print)

// Print slope and y-intercept
double slope_fid_p = fitFunc_fid_p->GetParameter(1);
double yIntercept_fid_p = fitFunc_fid_p->GetParameter(0);

//set up the canvas
TCanvas *c2 = new TCanvas("c2","scale field comparisons fid",1600,1200);

graph_fid_p->Draw("AP");

// Create and populate the legend
TLegend *legend_fid_p = new TLegend(0.5, 0.7, 0.8, 0.9); // Adjust these coordinates as needed
TString slopeStr_fid_p;
slopeStr_fid_p.Form("Slope: %0.5f", slope_fid_p);
TString interceptStr_fid_p;
interceptStr_fid_p.Form("Y-Intercept: %0.5f", yIntercept_fid_p);
legend_fid_p->AddEntry(graph_fid_p, slopeStr_fid_p, "l");
legend_fid_p->AddEntry(graph_fid_p, interceptStr_fid_p, "l");
legend_fid_p->Draw();

c2->Update();
c2->Write();

//Create graph
gStyle->SetEndErrorSize(0);
TGraph *graph_n = new TGraphErrors(num_files,graph_x_arr,graph_y_n_arr,graph_dx_arr,graph_dy_n_arr);
graph_n->SetMarkerStyle(20); // e.g., 20 is a full circle
graph_n->SetMarkerSize(1.0); // Size of the marker
graph_n->SetTitle(Form("%s %s Proton Neutron Means vs Scale Field",Kin.c_str(),targ.c_str()));

// Fit with a straight line
TF1 *fitFunc_n = new TF1("fitFunc_n", "pol1"); // pol1 for linear fit
graph_n->Fit(fitFunc_n, "Q"); // "Q" for quiet mode (no print)

// Print slope and y-intercept
double slope_n = fitFunc_n->GetParameter(1);
double yIntercept_n = fitFunc_n->GetParameter(0);

//set up the canvas
TCanvas *c3 = new TCanvas("c3","scale field comparisons",1600,1200);

graph_n->Draw("AP");

// Create and populate the legend
TLegend *legend_n = new TLegend(0.5, 0.7, 0.8, 0.9); // Adjust these coordinates as needed
TString slopeStr_n;
slopeStr_n.Form("Slope: %0.5f", slope_n);
TString interceptStr_n;
interceptStr_n.Form("Y-Intercept: %0.5f", yIntercept_n);
legend_n->AddEntry(graph_n, slopeStr_n, "l");
legend_n->AddEntry(graph_n, interceptStr_n, "l");
legend_n->Draw();

c3->Update();
c3->Write();

TGraph *graph_fid_n = new TGraphErrors(num_files,graph_x_fid_arr,graph_y_fid_n_arr,graph_dx_fid_arr,graph_dy_fid_n_arr);
graph_fid_n->SetMarkerStyle(20); // e.g., 20 is a full circle
graph_fid_n->SetMarkerSize(1.0); // Size of the marker
graph_fid_n->SetTitle(Form("%s %s Neutron Peak Means vs Scale Field, fidcut",Kin.c_str(),targ.c_str()));

// Fit with a straight line
TF1 *fitFunc_fid_n = new TF1("fitFunc_fid_n", "pol1"); // pol1 for linear fit
graph_fid_n->Fit(fitFunc_fid_n, "Q"); // "Q" for quiet mode (no print)

// Print slope and y-intercept
double slope_fid_n = fitFunc_fid_n->GetParameter(1);
double yIntercept_fid_n = fitFunc_fid_n->GetParameter(0);

//set up the canvas
TCanvas *c4 = new TCanvas("c4","scale field comparisons fid",1600,1200);

graph_fid_n->Draw("AP");

// Create and populate the legend
TLegend *legend_fid_n = new TLegend(0.5, 0.7, 0.8, 0.9); // Adjust these coordinates as needed
TString slopeStr_fid_n;
slopeStr_fid_n.Form("Slope: %0.5f", slope_fid_n);
TString interceptStr_fid_n;
interceptStr_fid_n.Form("Y-Intercept: %0.5f", yIntercept_fid_n);
legend_fid_n->AddEntry(graph_fid_n, slopeStr_fid_n, "l");
legend_fid_n->AddEntry(graph_fid_n, interceptStr_fid_n, "l");
legend_fid_n->Draw();

c4->Update();
c4->Write();



//output the info to some files
TString plotname = outputfilename;
plotname.ReplaceAll(".root",".pdf");
c0->Print(plotname + "(");
c2->Print(plotname + "");
c3->Print(plotname + "");
c4->Print(plotname + ")");
fout->Write();

watch->Stop();
// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end of sf_graph function
