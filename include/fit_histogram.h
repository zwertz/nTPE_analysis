#ifndef FIT_HISTOGRAM_H
#define FIT_HISTOGRAM_H

//Author: Ezekiel Wertz
//A class to handle fitting histograms that are used when comparing Data to MC for dx histograms.  It should just take in input information about the histograms and return information about the fit. To be used in later analysis. This is useful because of implementing Interpolate in most of these fitting procedures. This should take as input a type specifier to allow for modularity.
//// It allows you to fit the dx distribution from data to scaled version of montecarlo plus background,
//// This should remove the necessity of the simulation variables being globals

#include <string>
#include <vector>
class fit_histogram{

private:

string fitType;
string fitName;
string fitOptions;

TH1D *hist_data; //data histogram for dx containing both n and p info
TH1D *hist_p;  // proton simulation histogram
TH1D *hist_n; // neutron simulation histogram
TH1D *hist_bkgd; //candidate background histogram
TF1 *fitFunc; // Custom fit function

double xMin; // min x on fit range
double xMax; // max x on fit range
double scale_p; // proton scale factor from fit
double scale_p_err; //
double scale_n; // neutron scale factor from fit
double scale_n_err;
double shift_p;// how much the neutron histogram should be shifted from fit
double shift_p_err;
double shift_n;// how much the neutron histogram should be shifted from fit
double shift_n_err;
double Rsf; // ratio n/p scale from fit
double Rsf_err; //
double BGscale; // scale for the provided bg polynomial used in mc_p_n_BG_R
double BGscale_err;//
int polyorder; // order for the polynomial background used in mc_p_n_poly_R
double ChiSq; // chi^2 from fit
double NDF; // number of degrees of freedom from fit.
int paramCount; //number of parameters in the total fit

vector<pair<double,double>> fitParamsErrs;

vector<double> bkgd_params;
vector<double> bkgd_param_errs;

vector<pair<double,double>> fitAndFineFit(TH1D* histogram, string fitName, string fitFormula, int paramCount, double hcalfit_low, double hcalfit_high, const string& fitOptions = "RBMQ0");

public:
//constructor to handle polynomial backgrounds. This implementation inherently assumes you are shifting the histograms as well
fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, const char* fit_name, const char* fitType, int poly_ord, int param_count, double xmin, double xmax, const string& fit_options);

//constructor to handle user-input (essentially anticut or MC inelastic) backgrounds
fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, const char* fit_name, const char* fitType, vector<double> bkgd_coeffs, int param_count, double xmin, double xmax, const string& fit_options);

//constructor to handle gauss backgrounds. This implementation inherently assumes you are shifting the histograms as well
fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, const char* fit_name, const char* fit_type, int param_count, double xmin, double xmax,const string& fit_options);

//constructor to handle user-input (essentially anticut or MC inelastic) backgrounds
fit_histogram(TH1D *h_data,TH1D *h_p,TH1D *h_n, TH1D *h_bkgd, const char* fit_name, const char* fit_type, int param_count, double xmin, double xmax,const  string& fit_options);

//A function that should be able to generate any order polynomail. We are most interested in even order polynomials. A prerequiste is that polyorder must be not null. Or this will most likely crash
double polyN_fit(double *x, double *param);

//Function to generate Gaussian same as the total fit function
double Gauss(double *x, double *param);

//A function to mimic the user-input background
double InterpolateBG(double *x, double *param);

//destructor
//We will need to delete the dynamically allocated memory
~fit_histogram();

//No proton and neutron central value shifts. MC proton and neutron shape via the Interpolate function. Assumes a polynomial background and should be able to handle any order. 2nd and 4th order are most common.
double fitFull_polyBG(double *x, double *param);

//No proton and neutron central value shifts. MC proton and neutron shape via the Interpolate function. Not background at all
double fitFullNoBG(double *x, double *param);

//Shift proton and neutron central values. MC proton and neutron shape via the Interpolate function. Assumes a polynomial background and should be able to handle any order. 2nd and 4th order are most common.
double fitFullShiftNoBG(double *x, double *param);

double fitFullShift_scale_polyBG(double *x, double *param);

//Function that will take the histogram shape for the bkgd and interpolate
double fitFullShift_InterpolateBG(double *x, double *param);

//Implemented so Rsf is a fit parameter, which should better handle correllated error analysis
//Total function for MC + background. Independently shifts the neutron and proton peaks
double fitFullShift_polyBG(double *x, double *param);

double fitFullShift_gaussBG(double *x, double *param);

TH1D* get_hist_data();

TH1D* get_hist_p();

TH1D* get_hist_n();

string get_fitType();

string get_fitName();

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

double get_ChiSq();

double get_NDF();

int get_paramCount();

vector<pair<double,double>> get_fitParamsErrs();

vector<double> get_bkgd_params();

vector<double> get_bkgd_param_errs();

};//end class
#endif
