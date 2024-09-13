#include "../include/fit_histogram.h"

//Author: Ezekiel Wertz
//A companion implenetation for the fit_histogram class.A class to handle fitting histograms that are used when comparing Data to MC for dx histograms. This is class adapted from Maria, but is maybe a little more basic. It should just take in input information about the histograms and return information about the fit. To be used in later analysis. This is useful because of implementing Interpolate in most of these fitting procedures. This should take as input a type specifier to allow for modularity.
//// It allows you to fit the dx distribution from data to scaled version of montecarlo plus background,
//// This should remove the necessity of the simulation variables being globals

//constructor to handle polynomial backgrounds. This implementation inherently assumes you are shifting the histograms as well
fit_histogram::fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, string fit_name, const char *fitType, int poly_ord, int param_count, double xmin, double xmax, string& fit_options){
	
	
	//First check that the histograms are not null
	if(!h_data || !h_p || !h_n){
	cerr << "At least one of the histograms are null in the constructor!" << endl;
	}

	//convert the const char to a string
	string fitTypeString = fitType;

	//check that the fitType is implemented or not. Will need to update that if more fit types/functions are implemented
	if(fitTypeString != "fitFull_polyBG" || fitTypeString != "fitFullNoBG" || fitTypeString != "fitFullShift_polyBG" || fitTypeString !="fitFullShiftNoBG" || fitTypeString != "polyN_fit"){
	cerr << "The fitType: " << fitType << " is not properly implemented! Find the error or implement the fit!" << endl;
	}

	//Store the class variables we can right a way.
	hist_data = (TH1D*) h_data->Clone(Form("h_dx_data_%s",fitType));
	hist_p = (TH1D*) h_p->Clone(Form("h_dx_mc_p_%s",fitType));
	hist_n = (TH1D*) h_n->Clone(Form("h_dx_mc_n_%s",fitType));
	fitFunc = fitType;
	polyorder = poly_ord;
	paramCount = param_count;
	xMin = xmin;
	xMax = xmax;
	fitName = fit_name;
	fitOptions = fit_options;

	//Now fit the data histo with our function with a polynomial background
	vector<pair<double,double>> paramsErrsChi = fitAndFineFit(hist_data,fitName,fitFunc,paramCount,xMin,xMax,fitOptions);

	//Initialize the fit params errs vector. So we can store the original fit information
	for(int i=0; i<paramCount; i++){
	fitParamsErrs.first[i] = paramsErrsChi[i].first;
	fitParamsErrs.second[i] = paramsErrsChi[i].second;
	}
	
	//Now that we have the fit parameters. Make it easy to access the parts we will need
	scale_p = paramsErrsChi[0].first;
	scale_p_err = paramsErrsChi[0].second;

	Rsf = paramsErrsChi[1].first;
        Rsf_err = paramsErrsChi[1].second;

	scale_n = Rsf * scale_p;
	scale_n_err = scale_n*sqrt(pow(Rsf_error/Rsf,2)+pow(scale_p_err/scale_p,2)); //Treating these as indpendent parameters and adding in quadrature is not quite right. One should really do a dependent error analsis. However, this is not going to be use in an extraction. Could be worse.

	shift_p = paramsErrsChi[2].first;
	shift_p_err = paramsErrsChi[2].second;

	shift_n = paramsErrsChi[3].first;
        shift_n_err = paramsErrsChi[3].second;
	
	//Handle the polynomial parameters. This should be compatible with N order polynomial
	for(int j = 4; j < paramCount; j++){
	bkgd_params[j-4] = paramsErrsChi[j].first;
	bkgd_param_errs[j-4] = paramsErrsChi[j].second;
	}
	Chisq = paramsErrsChi[paramCount+1].first;
	NDF = paramsErrsChi[paramCount+1].second;

	//This should initizalie everything but bgscale and bgscale error
} //end of a constructor

//constructor to handle user-input (essentially anticut or MC inelastic) backgrounds
fit_histogram::fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, string fit_name, const char *fitType, vector<double> bkgd_coeffs, int param_count, double xmin, double xmax, string& fit_options){
        //First check that the histograms are not null
        if(!hist_data || !hist_p || !hist_n){
        cerr << "At least one of the histograms are null in the constructor!" << endl;
        }

        //convert the const char to a string
        string fitTypeString = fitType;

        //check that the fitType is implemented or not. Will need to update that if more fit types/functions are implemented
        if(fitTypeString != "fitFull_polyBG" || fitTypeString != "fitFullNoBG" || fitTypeString != "fitFullShift_polyBG" || fitTypeString !="fitFullShiftNoBG" || fitTypeString != "polyN_fit"){
        cerr << "The fitType: " << fitType << " is not properly implemented! Find the error or implement the fit!" << endl;
        }

	//Store the class variables we can right a way.
        hist_data = (TH1D*) h_data->Clone(Form("h_dx_data_%s",fitType));
        hist_p = (TH1D*) h_p->Clone(Form("h_dx_mc_p_%s",fitType));
        hist_n = (TH1D*) h_n->Clone(Form("h_dx_mc_n_%s",fitType));
        fitFunc = fitType;
        bkgd_params = bkgd_coeffs;
        paramCount = param_count;
        xMin = xmin;
        xMax = xmax;
        fitName = fit_name;
        fitOptions = fit_options;

	//Now fit the data histo with our function with a polynomial background
        vector<pair<double,double>> paramsErrsChi = fitAndFineFit(hist_data,fitName,fitFunc,paramCount,xMin,xMax,fitOptions);

	//Now that we have the fit parameters. Make it easy to access the parts we will need
        scale_p = paramsErrsChi[0].first;
        scale_p_err = paramsErrsChi[0].second;

        Rsf = paramsErrsChi[1].first;
        Rsf_err = paramsErrsChi[1].second;

        scale_n = Rsf * scale_p;
        scale_n_err = scale_n*sqrt(pow(Rsf_error/Rsf,2)+pow(scale_p_err/scale_p,2)); //Treating these as indpendent parameters and adding in quadrature is not quite right. One should really do a dependent error analsis. However, this is not going to be use in an extraction. Could be worse.

        shift_p = paramsErrsChi[2].first;
        shift_p_err = paramsErrsChi[2].second;

        shift_n = paramsErrsChi[3].first;
        shift_n_err = paramsErrsChi[3].second;

	BGscale = paramsErrsChi[4].first;
	BGscale_err = paramsErrsChi[4].second;


        //Handle the polynomial parameters. This should be compatible with N order polynomial
        for(int j = 5; j < paramCount; j++){
        bkgd_params[j-5] = paramsErrsChi[j].first;
        bkgd_param_errs[j-5] = paramsErrsChi[j].second;
        }
        Chisq = paramsErrsChi[paramCount+1].first;
        NDF = paramsErrsChi[paramCount+1].second;

        //This should initizalie everything but poly order

}// end of a constructor

//destructor
//We will need to delete the dynamically allocated memory
fit_histogram::~fit_histogram(){
delete hist_data;
delete hist_p;
delete hist_n;
}

//Implemented no longer really used
//Implemented so Rsf is a fit parameter, which should better handle correllated error analysis
//Total function for MC + background.
double fit_histogram::fitFull_polyBG(double *x, double *param){
double dx = x[0];
double proton_scale = param[0];
double R_sf = param[1];

double proton = hist_p->Interpolate(dx);
double neutron = hist_n->Interpolate(dx);

//This is equivalent to neutron_scale*neutron + proton_scale*proton + background
double fitFull = proton_scale*( R_sf*neutron + proton) + polyN_fit(x, &param[2]);
return fitFull;
}

//Implemented no longer really used
//Full fit of proton and neutron peaks. No background
double fit_histogram::fitFullNoBG(double *x, double *param){
double dx = x[0];
double proton_scale = param[0];
double R_sf = param[1];

double proton = hist_p->Interpolate(dx);
double neutron = hist_n->Interpolate(dx);

//This is equivalent to neutron_scale*neutron + proton_scale*proton
double fitFull = proton_scale*( R_sf*neutron + proton);
return fitFull;
}

//Implemented so Rsf is a fit parameter, which should better handle correllated error analysis
//Total function for MC + background. Independently shifts the neutron and proton peaks
double fit_histogram::fitFullShift_polyBG(double *x, double *param){
//MC float parameters
double proton_scale = param[0];
double R_sf = param[1];

double dx_shift_p = param[2];
double dx_shift_n = param[3];

//Apply the shifts before any interpolation
double proton = hist_p->Interpolate(x[0] - dx_shift_p);
double neutron = hist_n->Interpolate(x[0] - dx_shift_n);

//The total function is the proton and neutron peaks + a 4th order polynomial for background
//This is equivalent to neutron_scale*neutron + proton_scale*proton + background
double fitFullshift = proton_scale*( R_sf*neutron + proton) + polyN_fit(x, &param[4]);
return fitFullshift;
}

//Implemented so Rsf is a fit parameter, which should better handle correllated error analysis
//Total function for MC + background. Independently shifts the neutron and proton peaks
double fit_histogram::fitFullShift_scale_polyBG(double *x, double *param){
//MC float parameters
double proton_scale = param[0];
double R_sf = param[1];

double dx_shift_p = param[2];
double dx_shift_n = param[3];

double bg_scale = param[5];

//Apply the shifts before any interpolation
double proton = hist_p->Interpolate(x[0] - dx_shift_p);
double neutron = hist_n->Interpolate(x[0] - dx_shift_n);

//The total function is the proton and neutron peaks + a Nth order polynomial which has some user-input shape
//This is equivalent to neutron_scale*neutron + proton_scale*proton + background
double val = proton_scale*( R_sf*neutron + proton);

	//loop over the known bkgd coefficients and scale the background accordingly
	for(int i =0, i < bkgd_params.size(); i++ ){
	val += bg_scale * (bkgd_params[i] * pow(x[0],i));
	}

return val;
}


//Implemented so Rsf is a fit parameter, which should better handle correllated error analysis
//Full fit of proton and neutron peaks. No background. Independently shifts the neutron and proton peaks
double fit_histogram::fitFullShiftNoBG(double *x, double *param){
//MC float parameters
double proton_scale = param[0];
double R_sf = param[1];

double dx_shift_p = param[2];
double dx_shift_n = param[3];

//Apply the shifts before any interpolation
double proton =  hist_p->Interpolate(x[0] - dx_shift_p);
double neutron = hist_n->Interpolate(x[0] - dx_shift_n);

//This is equivalent to neutron_scale*neutron + proton_scale*proton
double fitFullshift = proton_scale*( R_sf*neutron + proton);
return fitFullshift;
}

//A function that should be able to generate any order polynomail. We are most interested in even order polynomials. A prerequiste is that polyorder must be not null. Or this will most likely crash
double fit_histogram::polyN_fit(double *x, double *param){

double poly_func = 0.0;

	for(int i=0; i<=polyorder; i++){
	poly_func += (param[i]*pow(x[0],i));
	}
return poly_func;
}




//function to fit the histogram then fine fit and return all desired parameters.
vector<pair<double,double>> fit_histogram::fitAndFineFit(TH1D* histogram, const string& fitName, const string& fitFormula, int paramCount, double hcalfit_low, double hcalfit_high, const string& fitOptions = "RBMQ0"){

//make the vector we will eventually return. It will hold the parameters and the last pair will be chi2 and ndf
vector<pair<double,double>> paramsAndErrs(paramCount+1);

//Make a fit based on the information provided
TF1* datFit = new TF1(fitName.c_str(),fitFormula.c_str(),hcalfit_low,hcalfit_high,paramCount);

//check if the fitFormula has shifted parameters
	if(fitFormula.contains("Shift")){
	//verify that hte shifts are going to be within the fit parameter
	datFit->SetParLimits(2,hcalfit_low,hcalfit_high);
	datFit->SetParLimits(3,hcalfit_low,hcalfit_high);
	|
//Fit the provided histogram with the fit we have
histogram->Fit(datFit,fitOptions.c_str());

//vector to hold fine fit initial parameters
vector<double> fineFitParams(paramCount);
        //loop over the information from the fit and store it in our parameter and error vector
        for(int j=0; j<paramCount; ++j){
        paramsAndErrs[j].first = datFit->GetParameter(j);
        paramsAndErrs[j].second = datFit->GetParError(j);
        fineFitParams[j] = paramsAndErrs[j].first;
        //cout << datFit->GetParameter(j) << endl;
        }

//Fine fit to get a better set of parameters
TF1* datFineFit = new TF1((fitName+"_fine").c_str(),fitFormula.c_str(),hcalfit_low,hcalfit_high,paramCount);

//set the fine fit intiall parameters to be the ones we just found
datFineFit->SetParameters(fineFitParams.data());

//Fit the histogram again but with these better intial parameters
histogram->Fit(datFineFit,fitOptions.c_str());

        //update the parameters and errors with the ones we just found from fine fit
        for(int k=0; k<paramCount; ++k){
        paramsAndErrs[k].first = datFineFit->GetParameter(k);
        paramsAndErrs[k].second = datFineFit->GetParError(k);
        //Diagnostic
        cout << datFineFit->GetParameter(k) << endl;
        }

//don't like these lines because it kind of circumvents the idea of a return
//get the Chisquare and NDF from the fit which is useful for analysis
paramsAndErrs[paramCount+1].first = datFineFit->GetChisquare();
paramsAndErrs[paramCount+1].second = datFineFit->GetNDF();

//cleanup
delete datFit;
delete datFineFit;
return paramsAndErrs;
}//end fit and fine fit func

//Getter functions for class variables
string& get_fitType();

string& get_fitName();

string& get_fitOptions();

double get_xMin();

double get_xMax();

double get_scale_p();

double get_scale_p_err();

double get_scale_n();

double get_scale_n_err();

double get_shift_p();

double get_shift_p_err();

double get_shift_n();

double get_shift_n_err();

double get_Rsf();

double get_Rsf_err();

double get_BGscale();

double get_BGscale_err();

int get_polyorder();

double get_ChiSq;()

double get_NDF();

int get_paramCount();

vector<pair<double,double>> get_fitParamsErrs();

vector<double> get_bkgd_params();

vector<double> get_bkgd_param_errs();
