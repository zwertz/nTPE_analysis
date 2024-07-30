#ifndef BRADFORD_PARAM_H
#define BRADFORD_PARAM_H

#include "/work/halla/sbs/jboyd/include/experimental_constants.h"

vector<double> bradford_gmn_a = {1.00, 1.81, 0.00};
vector<double> bradford_gmn_b = {-9.86, 305.00, -758.00, 802.00};

Double_t bradford_GMn_parameterization( Double_t *x, Double_t *par_G ){
	double GMn = calc_bradford_parameterization( 4, x[0], true );
	return GMn;
}

void general_plot_bradford(TCanvas *c_to_plot_to){

}

#endif
