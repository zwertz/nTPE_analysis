#ifndef KELLY_PARAM_ONLY_H
#define KELLY_PARAM_ONLY_H

#include "/work/halla/sbs/jboyd/include/calc_errors.h"
#include "/work/halla/sbs/jboyd/include/experimental_constants.h"

double par_G[4], par_kelly_G_error[4];

double par_GEp[5] = {1.0, -0.24, 10.98, 12.82, 21.97};
double par_GEp_error[5] = {0.0, 0.12, 0.19, 1.1, 6.8};

double par_GMp[5] = {1.0, 0.12, 10.97, 18.86, 6.55};
double par_GMp_error[5] = {0.0, 0.04, 0.11, 0.28, 1.2};

double par_GMn[5] = {1.0, 2.33, 14.72, 24.20, 84.1};
double par_GMn_error[5] = {0.0, 1.4, 1.7, 9.8, 41};

double par_GEn[2] = {1.70, 3.30};
double par_GEn_error[2] = {0.04, 0.32};

Double_t calc_kelly_parameterization( int FF, double q2, bool normed = true ){
	int n = 1;
	
	double Mp = 0.938272;
	double Mn = 0.939565;
	double MN = 0.0;
	if( FF == 1 || FF == 2 ){
		MN = Mp;
	}
	if( FF == 3 || FF == 4 ){
		MN = Mn;
	}

	double tau = q2/(4*MN*MN);
	double delta2 = 0.71;

	double G_D = pow( (1 + (q2/delta2)), -2);

	double G_calc = 0.0;

	double a_sum = 0.0, b_sum = 0.0;

	double par_G_select[5];

//----GEp
	if( FF == 1 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GEp[i];
		}
	}
//----GMp
	if( FF == 2 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GMp[i];
		}
	}
//----GEn
	if( FF == 3 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GEn[i];
		}
		for( int i = 2; i < 6; i++ ){
			par_G_select[i] = 0.0;
		}
	}
//----GMn
	if( FF == 4 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GMn[i];
		}
	}

//----Calculation
	for( int i = 0; i <= n; i++ ){
		a_sum += par_G_select[i]*pow( tau, i );
	}

	for( int i = 1; i <= n+2; i++ ){
		b_sum += par_G_select[i+1]*pow( tau, i);
	}
	b_sum += 1.0;

	//--- For everything but GEn
	if( FF != 3 ){
		G_calc = a_sum/b_sum;	
	}

	//--- For GEn
	if( FF == 3 ){
		G_calc = ( ( par_GEn[0]*tau )/( 1 + par_GEn[1]*tau ) )*G_D;
	}


//-----Return normalized value or not
	if( !normed ){
		if( FF == 2 ){ //GMp
			return G_calc*mu_p*mu_p;
		}
		if( FF == 4 ){ //GMn
			return G_calc*mu_n;
		}
		else{
			return G_calc;	
		}
			
	}
	if( normed ){
		if( FF == 3 ){
			return G_calc*mu_n/G_D;
		}
		else{
			return G_calc/G_D;			
		}
	
	}

}

Double_t calc_kelly_errors( int FF, double q2, bool normed = true ){
	double delta2 = 0.71;
	double G_D = pow( (1 + (q2/delta2)), -2);

	double Mp = 0.938272;
	double Mn = 0.939565;
	double MN = 0.0;
	if( FF == 1 || FF == 2 ){
		MN = Mp;
	}
	if( FF == 3 || FF == 4 ){
		MN = Mn;
	}

	double tau = q2/(4*MN*MN);

	double numerator = 0.0;
	double numerator_error = 0.0;

	double denominator = 0.0;
	double denominator_error = 0.0;

	double division = 0.0;
	double division_error = 0.0;

	double final_error = 0.0;

	double par_G_select[5], par_G_select_error[5];

//----GEp
	if( FF == 1 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GEp[i];
			par_G_select_error[i] = par_GEp_error[i];
		}
	}
//----GMp
	if( FF == 2 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GMp[i];
			par_G_select_error[i] = par_GMp_error[i];
		}
	}
//----GEn
	if( FF == 3 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GEn[i];
			par_G_select_error[i] = par_GEn_error[i];
		}
		for( int i = 2; i < 6; i++ ){
			par_G_select[i] = 0.0;
			par_G_select_error[i] = 0.0;
		}
	}
//----GMn
	if( FF == 4 ){
		for( int i = 0; i < 6; i++ ){
			par_G_select[i] = par_GMn[i];
			par_G_select_error[i] = par_GMn_error[i];
		}
	}

	if( FF != 3 ){
		numerator = par_G_select[1];
		numerator_error = par_G_select_error[1];

		denominator = 1 + par_G_select[2] + par_G_select[3] + par_G_select[4];
		denominator_error = sqrt( pow(par_G_select_error[2], 2) + pow(par_G_select_error[3], 2) + pow(par_G_select_error[4], 2) );

		division = numerator/denominator;
		division_error = CalculateErrorMultiplicationDivision( numerator, numerator_error, denominator, denominator_error, division);
	}

	if( FF == 3 ){
		numerator = par_GEn[0];
		numerator_error = par_GEn_error[0];

		denominator = 1+ par_GEn[1];
		denominator_error = par_GEn_error[1];

		division = G_D*numerator/denominator;
		division_error = CalculateErrorMultiplicationDivision( numerator, numerator_error, denominator, denominator_error, division);
		division_error = division_error;		
	}

	if( normed ){
		division_error = division_error/G_D;
	}

	// cout << "division error: " << division_error << endl;
	return division_error;
}

#endif