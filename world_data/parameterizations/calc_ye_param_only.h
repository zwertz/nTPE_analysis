#ifndef CALC_KELLY_PARAMETERIZATION_H
#define CALC_KELLY_PARAMETERIZATION_H

#include "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/world_data/parameterizations/GetFFfunctions.h"
//-------------------------------
//FF 0: All, 1: GEp, 2: GMp, 3: GEn, 4: GMn

Double_t calc_ye_parameterization( int FF, double q2, bool normed = true ){

	double delta_norm = 0.71;
	double G_D_norm = 1.0/( pow( (1 + (q2/delta_norm)), 2 ) );
	double mu_neutron = 1.9130427;
	double mu_proton = 2.792847351;

	double FF_lookup = 0.0;

	FF_lookup = GetFFvalue( FF, q2 );

	if( normed == false ){
		if( FF == 1 || FF == 3 ){
			FF_lookup = FF_lookup*G_D_norm;
		}
		else if( FF == 2 ){
			FF_lookup = FF_lookup*G_D_norm*mu_proton;
		}
		else{
			FF_lookup = FF_lookup*G_D_norm*mu_neutron;
		}
	}

	return FF_lookup;

}

Double_t calc_ye_error( int FF, double q2, bool normed = true){

	double delta_norm = 0.71;
	double G_D_norm = 1/( pow( (1 + (q2/delta_norm)), 2 ) );
	double mu_neutron = 1.91;

	double FF_error_lookup = 0.0;

	FF_error_lookup = GetFFerror( FF, q2 );

	if( normed == false ){
		if( FF == 1 || FF == 2 ){
			FF_error_lookup = FF_error_lookup*G_D_norm;
		}
		if( FF == 3 || FF == 4 ){
			FF_error_lookup = FF_error_lookup*G_D_norm*mu_neutron;
		}
	}


	return FF_error_lookup;

}

#endif