#include "../include/fit_histogram.h"

//Author: Ezekiel Wertz
//A companion implenetation for the fit_histogram class.A class to handle fitting histograms that are used when comparing Data to MC for dx histograms. This is class adapted from Maria, but is maybe a little more basic. It should just take in input information about the histograms and return information about the fit. To be used in later analysis. This is useful because of implementing Interpolate in most of these fitting procedures. This should take as input a type specifier to allow for modularity.
//// It allows you to fit the dx distribution from data to scaled version of montecarlo plus background,
//// This should remove the necessity of the simulation variables being globals

//constructor to handle polynomial backgrounds. This implementation inherently assumes you are shifting the histograms as well
fit_histogram::fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, const char* fit_name, const char* fit_type, int poly_ord, int param_count, double xmin, double xmax,const string& fit_options){
	
	
	//First check that the histograms are not null
	if(!h_data || !h_p || !h_n){
	cerr << "At least one of the histograms are null in the constructor!" << endl;
	}

	fitType = fit_type;
	fitName = fit_name;

	//check that the fitType is implemented or not. Will need to update that if more fit types/functions are implemented
	if(fitType != "fitFull_polyBG" && fitType != "fitFullNoBG" && fitType != "fitFullShift_polyBG" && fitType !="fitFullShiftNoBG" && fitType != "polyN_fit"){
	cerr << "The fitType: " << fitType << " is not properly implemented! Find the error or implement the fit!" << endl;
	}

	//Store the class variables we can right a way.
	hist_data = (TH1D*) h_data->Clone(Form("h_dx_data_%s",fitName.c_str()));
	hist_p = (TH1D*) h_p->Clone(Form("h_dx_mc_p_%s",fitName.c_str()));
	hist_n = (TH1D*) h_n->Clone(Form("h_dx_mc_n_%s",fitName.c_str()));
	polyorder = poly_ord;
	paramCount = param_count;
	xMin = xmin;
	xMax = xmax;
	fitOptions = fit_options;

	
	//Now fit the data histo with our function with a polynomial background
	vector<pair<double,double>> paramsErrsChi = fitAndFineFit(hist_data,fitName,fitType,paramCount,xMin,xMax,fitOptions);

	//Initialize the fit params errs vector. So we can store the original fit information
	for(int l=0; l<paramCount; ++l){
	pair<double,double> mydoub(paramsErrsChi[l].first,paramsErrsChi[l].second);
	fitParamsErrs.push_back(mydoub);
	
	}
	
	//Now that we have the fit parameters. Make it easy to access the parts we will need
	scale_p = paramsErrsChi[0].first;
	scale_p_err = paramsErrsChi[0].second;
	
	Rsf = paramsErrsChi[1].first;
        Rsf_err = paramsErrsChi[1].second;

	scale_n = Rsf * scale_p;
	scale_n_err = scale_n*sqrt(pow(Rsf_err/Rsf,2)+pow(scale_p_err/scale_p,2)); //Treating these as indpendent parameters and adding in quadrature is not quite right. One should really do a dependent error analsis. However, this is not going to be use in an extraction. Could be worse.

	shift_p = paramsErrsChi[2].first;
	shift_p_err = paramsErrsChi[2].second;

	shift_n = paramsErrsChi[3].first;
        shift_n_err = paramsErrsChi[3].second;
	
	//Handle the polynomial parameters. This should be compatible with N order polynomial
	for(int j = 4; j < paramCount; j++){
	bkgd_params.push_back(paramsErrsChi[j].first);
	bkgd_param_errs.push_back(paramsErrsChi[j].second);
	}
	
	ChiSq = paramsErrsChi[paramCount].first;
	NDF = paramsErrsChi[paramCount].second;
	
	//This should initizalie everything but bgscale and bgscale error
} //end of a constructor

//constructor to handle user-input (essentially anticut or MC inelastic) backgrounds
fit_histogram::fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, const char* fit_name, const char* fit_type, vector<double> bkgd_coeffs, int param_count, double xmin, double xmax,const  string& fit_options){
        //First check that the histograms are not null
        if(!hist_data || !hist_p || !hist_n){
        cerr << "At least one of the histograms are null in the constructor!" << endl;
        }
	
	fitType = fit_type;
	fitName = fit_name;

        //check that the fitType is implemented or not. Will need to update that if more fit types/functions are implemented
        if(fitType != "fitFull_polyBG" && fitType != "fitFullNoBG" && fitType != "fitFullShift_polyBG" && fitType !="fitFullShiftNoBG" && fitType != "polyN_fit"){
        cerr << "The fitType: " << fitType << " is not properly implemented! Find the error or implement the fit!" << endl;
        }

	//Store the class variables we can right a way.
        hist_data = (TH1D*) h_data->Clone(Form("h_dx_data_%s",fitName.c_str()));
        hist_p = (TH1D*) h_p->Clone(Form("h_dx_mc_p_%s",fitName.c_str()));
        hist_n = (TH1D*) h_n->Clone(Form("h_dx_mc_n_%s",fitName.c_str()));
        bkgd_params = bkgd_coeffs;
        paramCount = param_count;
        xMin = xmin;
        xMax = xmax;
        fitOptions = fit_options;

	//Now fit the data histo with our function with a polynomial background
        vector<pair<double,double>> paramsErrsChi = fitAndFineFit(hist_data,fitName.c_str(),fitType.c_str(),paramCount,xMin,xMax,fitOptions);

	 //Initialize the fit params errs vector. So we can store the original fit information
        for(int l=0; l<paramCount; ++l){
        pair<double,double> mydoub(paramsErrsChi[l].first,paramsErrsChi[l].second);
        fitParamsErrs.push_back(mydoub);

        }

	//Now that we have the fit parameters. Make it easy to access the parts we will need
        scale_p = paramsErrsChi[0].first;
        scale_p_err = paramsErrsChi[0].second;

        Rsf = paramsErrsChi[1].first;
        Rsf_err = paramsErrsChi[1].second;

        scale_n = Rsf * scale_p;
        scale_n_err = scale_n*sqrt(pow(Rsf_err/Rsf,2)+pow(scale_p_err/scale_p,2)); //Treating these as indpendent parameters and adding in quadrature is not quite right. One should really do a dependent error analsis. However, this is not going to be use in an extraction. Could be worse.

        shift_p = paramsErrsChi[2].first;
        shift_p_err = paramsErrsChi[2].second;

        shift_n = paramsErrsChi[3].first;
        shift_n_err = paramsErrsChi[3].second;

	BGscale = paramsErrsChi[4].first;
	BGscale_err = paramsErrsChi[4].second;


        //Handle the polynomial parameters. This should be compatible with N order polynomial
        for(int j = 5; j < paramCount; j++){
        bkgd_params.push_back(paramsErrsChi[j].first);
        bkgd_param_errs.push_back(paramsErrsChi[j].second);
        }
        ChiSq = paramsErrsChi[paramCount].first;
        NDF = paramsErrsChi[paramCount].second;

        //This should initizalie everything but poly order

}// end of a constructor

//destructor
//We will need to delete the dynamically allocated memory
fit_histogram::~fit_histogram(){
delete hist_data;
delete hist_p;
delete hist_n;
delete fitFunc;
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
	for(int i =0; i < bkgd_params.size(); i++ ){
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
vector<pair<double,double>> fit_histogram::fitAndFineFit(TH1D* histogram,  string fitName,  string fitFormula, int paramCount, double hcalfit_low, double hcalfit_high, const string& fitOptions = "RBMQ0"){

//make the vector we will eventually return. It will hold the parameters and the last pair will be chi2 and ndf
vector<pair<double,double>> paramsAndErrs(paramCount+1);

//Will need io make sure new fitTypes get implemented if necessary
//Explicilty implement the lambda function, need to see if this works, these are tricky
//Make a fit based on the information provided
TF1 *datFit = new TF1(fitName.c_str(),[this,fitFormula](double *x, double *par) -> double {
	if(fitFormula == "fitFull_polyBG"){
	return this->fitFull_polyBG(x, par);
	}else if(fitFormula == "fitFullNoBG"){
	return this->fitFullNoBG(x, par);
	}else if(fitFormula =="fitFullShift_polyBG"){
	return this->fitFullShift_polyBG(x, par);
	}else if(fitFormula =="fitFullShiftNoBG"){
	return this->fitFullShiftNoBG(x, par);
	}else{
	cout << "The Lambda function you are trying to implement " << fitFormula << " is no good! Figure it out now!" << endl;
	return this->fitFullShift_polyBG(x, par);
	}
  },hcalfit_low,hcalfit_high,paramCount);


string mystring("Shift");

//check if the fitFormula has shifted parameters
	if(fitFormula.find(mystring) != string::npos){
	
	//verify that hte shifts are going to be within the fit parameter
	datFit->SetParLimits(2,hcalfit_low,hcalfit_high);
	datFit->SetParLimits(3,hcalfit_low,hcalfit_high);
	}

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

//Fine fit to get a better set of parameters. Make it equal to the fit function in the class. This function should only get called by the constructor
fitFunc = new TF1((fitName+"_fine").c_str(),[this,fitFormula](double *x, double *par) -> double {
        if(fitFormula == "fitFull_polyBG"){
        return this->fitFull_polyBG(x, par);
        }else if(fitFormula == "fitFullNoBG"){
        return this->fitFullNoBG(x, par);
        }else if(fitFormula =="fitFullShift_polyBG"){
        return this->fitFullShift_polyBG(x, par);
        }else if(fitFormula =="fitFullShiftNoBG"){
        return this->fitFullShiftNoBG(x, par);
        }else{
        cout << "The Lambda function you are trying to implement " << fitFormula << " is no good! Figure it out now!" << endl;
        return this->fitFullShift_polyBG(x, par);
        }
	},hcalfit_low,hcalfit_high,paramCount);

//set the fine fit intiall parameters to be the ones we just found
fitFunc->SetParameters(fineFitParams.data());

//Fit the histogram again but with these better intial parameters
histogram->Fit(fitFunc,fitOptions.c_str());

        //update the parameters and errors with the ones we just found from fine fit
        for(int k=0; k<paramCount; ++k){
        paramsAndErrs[k].first = fitFunc->GetParameter(k);
        paramsAndErrs[k].second = fitFunc->GetParError(k);
        //Diagnostic
        cout << fitFunc->GetParameter(k) << endl;
        }

//don't like these lines because it kind of circumvents the idea of a return
//get the Chisquare and NDF from the fit which is useful for analysis
paramsAndErrs[paramCount].first = fitFunc->GetChisquare();
paramsAndErrs[paramCount].second = fitFunc->GetNDF();


return paramsAndErrs;
}//end fit and fine fit func

//Getter functions for class variables
string fit_histogram::get_fitType(){return fitType;}

string fit_histogram::get_fitName(){return fitName;}

string& fit_histogram::get_fitOptions(){return fitOptions;}

double fit_histogram::get_xMin(){return xMin;}

double fit_histogram::get_xMax(){return xMax;}

double fit_histogram::get_scale_p(){return scale_p;}

double fit_histogram::get_scale_p_err(){return scale_p_err;}

double fit_histogram::get_scale_n(){return scale_n;}

double fit_histogram::get_scale_n_err(){return scale_n_err;}

double fit_histogram::get_shift_p(){return shift_p;}

double fit_histogram::get_shift_p_err(){return shift_p_err;}

double fit_histogram::get_shift_n(){return shift_n;}

double fit_histogram::get_shift_n_err(){return shift_n_err;}

double fit_histogram::get_Rsf(){return Rsf;}

double fit_histogram::get_Rsf_err(){return Rsf_err;}

double fit_histogram::get_BGscale(){return BGscale;}

double fit_histogram::get_BGscale_err(){return BGscale_err;}

int fit_histogram::get_polyorder(){return polyorder;}

double fit_histogram::get_ChiSq(){return ChiSq;}

double fit_histogram::get_NDF(){return NDF;}

int fit_histogram::get_paramCount(){return paramCount;}

vector<pair<double,double>> fit_histogram::get_fitParamsErrs(){return fitParamsErrs;}

vector<double> fit_histogram::get_bkgd_params(){return bkgd_params;}

vector<double> fit_histogram::get_bkgd_param_errs(){return bkgd_param_errs;}

TH1D* fit_histogram::get_hist_data(){return hist_data;}

TH1D* fit_histogram::get_hist_p(){return hist_p;}

TH1D* fit_histogram::get_hist_n(){return hist_n;}
