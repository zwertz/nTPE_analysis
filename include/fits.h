#ifndef FITS_H
#define FITS_h

#include <vector>
#include <utility>

//Author: Ezekiel Wertz
//namespace to hold common fit definitions
namespace fits{

//4th-order polynomial
double poly4_fit(double *x, double *param);

//function to fit the histogram then fine fit and return all desired parameters.
vector<pair<double,double>> fitAndFineFit(TH1D* histogram, const string& fitName, const string& fitFormula, int paramCount, double hcalfit_low, double hcalfit_high, pair<double,double>& fitqual, const string& fitOptions = "RBMQ0");

//Get the variance on a fit from the quadrature sum of the parameter errors
double getFitError(TF1* fit);

//Get the variance on a fit from the quadrature sum of the diagonal entries. Which the diagonal entries should be the variance on each parameter.
double getFitError(TFitResultPtr fit_ptr);

}//end namespace
