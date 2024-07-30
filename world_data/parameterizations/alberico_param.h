#ifndef ALBERICO_PARAM_H
#define ALBERICO_PARAM_H

double bMp1 = 12.31;
double bMp2 = 25.57;
double bMp3 = 30.61;
double aMp1 = 1.09;

double bEp1 = 11.12;
double bEp2 = 15.16;
double bEp3 = 21.25;
double aEp1 = -0.19;

double aEn = -0.1;
double bEn1 = 2.83;
double bEn2 = 0.43;

double bMn1 = 21.30;
double bMn2 = 77;
double bMn3 = 238;
double aMn1 = 8.28;

double calc_alberico_param( int FF = 4, double Q2 = 4.5, bool normalize = true ){
	// FF = 1 -> GEp
	// FF = 2 -> GMp
	// FF = 3 -> GEn
	// FF = 4 -> GMn

	double proton_mass = 0.938272;
	double neutron_mass = 0.939565;
	double tau = 0.0;

	double mu_proton = 1.91;
	double mu_neutron = 2.79;
	double mu = 0.0;

	if( ( FF == 1 ) || ( FF == 2 ) ){
		mu = mu_proton;
		tau = 1/(4*pow(mu_proton, 2));
	}
	else{
		mu = mu_neutron;
		tau = 1/(4*pow(mu_neutron, 2));
	}

	double delta_norm = 0.71;
	double G_D_norm = 1.0/pow( (1 + (Q2/delta_norm)), 2);


	double numerator = 0.0, denominator = 0.0, FF_calc = 0.0;

	if( FF == 1 ){
		numerator = 1.0 + aEp1*tau;
		denominator = 1 + bEp1*tau + bEp2*tau*tau + bEp3*tau*tau*tau;

		FF_calc = numerator/denominator;
	}
	if( FF == 2 ){
		numerator = 1.0 + aMp1*tau;
		denominator = 1 + bMp1*tau + bMp2*tau*tau + bMp3*tau*tau*tau;

		FF_calc = numerator/denominator;
	}
	if( FF == 3 ){
		numerator = aEn;
		double denominator_1 = pow( (1 + bEn1*Q2), 2);
		double denominator_2 = pow( (1 + bEn2*Q2), 2);


		FF_calc = (numerator/denominator_1) - (numerator/denominator_2);
	}
	if( FF == 4 ){
		numerator = 1.0 + aMn1*tau;
		denominator = 1 + bMn1*tau + bMn2*tau*tau + bMn3*tau*tau*tau;

		FF_calc = numerator/denominator;
	}

	if( normalize ){
		if( ( FF == 1 ) || ( FF == 2 ) ){
			FF_calc = FF_calc/G_D_norm;
		}
		else{
			FF_calc = FF_calc/(mu*G_D_norm);
		}
	}

	return FF_calc;
}


#endif