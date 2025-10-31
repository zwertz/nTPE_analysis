#ifndef FITS_H
#define FITS_H

#include <vector>
#include <utility>

//Author: Ezekiel Wertz
//namespace to hold common fit definitions
//Largely outdated and depreciated due to the fit_histogram class. Ideally I should come back to this class and see if anything still uses it.
namespace fits{

//Get the variance on a fit from the quadrature sum of the parameter errors
double getFitError(TF1* fit);

//Get the variance on a fit from the quadrature sum of the diagonal entries. Which the diagonal entries should be the variance on each parameter.
double getFitError(TFitResultPtr fit_ptr);


}//end namespace
#endif
