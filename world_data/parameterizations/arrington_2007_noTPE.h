#ifndef ARRINGTON_2007_NOTPE_H
#define ARRINGTON_2007_NOTPE_H

#include "/work/halla/sbs/jboyd/include/utility_functions.h"
#include "/work/halla/sbs/jboyd/include/plot_pimping.h"

double arr07_par_G[4], arr07_par_G_error[4];

double am1 = -1.465;
double am2 = 1.260;
double am3 = 0.262;
double bm1 = 9.627;
double bm2 = 0.000;
double bm3 = 0.000;
double bm4 = 11.179;
double bm5 = 13.245;

double ae1 = 3.439;
double ae2 = -1.602;
double ae3 = 0.068;
double be1 = 15.055;
double be2 = 48.061;
double be3 = 99.304;
double be4 = 0.012;
double be5 = 8.650;

double arr07_proton_mass = 0.938272, arr07_mu_p, arr07_delta, arr07_G_D, tau, numerator = 0.0, denominator = 0.0;

double FF_calc_val = 0.0;

vector<vector<double>> arr07_GEp_vals_and_errs = { 
	// Q2   GE/GD   Err
	{0.007, 1.000, 0.006}, {0.012, 0.996, 0.006}, {0.017, 0.995, 0.003}, {0.022, 0.993, 0.003}, {0.030, 0.988, 0.007},
	{0.038, 0.987, 0.004}, {0.041, 0.974, 0.007}, {0.048, 0.057, 0.005}, {0.057, 0.980, 0.005}, {0.061, 0.990, 0.005},
	{0.069, 0.977, 0.007}, {0.081, 0.977, 0.005}, {0.098, 0.973, 0.007}, {0.115, 0.974, 0.011}, {0.138, 0.964, 0.012},
	{0.171, 0.986, 0.012}, {0.199, 0.979, 0.011}, {0.234, 0.978, 0.011}, {0.273, 0.955, 0.010}, {0.304, 0.960, 0.010},
	{0.350, 0.939, 0.014}, {0.390, 0.965, 0.011}, {0.428, 0.961, 0.015}, {0.473, 0.970, 0.011}, {0.528, 0.984, 0.013},
	{0.584, 0.967, 0.013}, {0.622, 0.969, 0.014}, {0.689, 0.981, 0.023}, {0.779, 0.965, 0.011}, {0.853, 0.916, 0.022},
	{0.979, 0.933, 0.049}, {1.020, 0.920, 0.017}, {1.170, 0.922, 0.014}, {1.330, 0.936, 0.047}, {1.520, 0.889, 0.022},
	{1.740, 0.844, 0.022}, {1.830, 0.856, 0.038}, {2.070, 0.837, 0.038}, {2.500, 0.744, 0.038}, {2.690, 0.723, 0.075},
	{2.900, 0.676, 0.039}, {3.090, 0.673, 0.118}, {3.240, 0.673, 0.112}, {3.710, 0.652, 0.040}, {4.110, 0.546, 0.059},
	{5.010, 0.409, 0.057}, {5.850, 0.287, 0.095}
};

vector<vector<double>> arr07_GMp_vals_and_errs = {
	{0.007, 0.715, 0.393}, {0.012, 0.966, 0.152}, {0.017, 0.984, 0.031}, {0.022, 0.953, 0.019}, {0.030, 0.983, 0.052},
	{0.038, 0.985, 0.013}, {0.041, 1.016, 0.027}, {0.048, 0.982, 0.009}, {0.057, 0.992, 0.010}, {0.061, 0.978, 0.015},
	{0.069, 0.991, 0.016}, {0.081, 0.980, 0.008}, {0.098, 0.977, 0.007}, {0.115, 0.978, 0.007}, {0.138, 0.982, 0.008},
	{0.171, 0.977, 0.007}, {0.199, 0.981, 0.006}, {0.234, 0.978, 0.012}, {0.273, 0.972, 0.008}, {0.304, 0.981, 0.007},
	{0.350, 0.985, 0.013}, {0.390, 0.993, 0.008}, {0.428, 1.003, 0.014}, {0.473, 0.995, 0.008}, {0.528, 0.996, 0.009},
	{0.584, 1.007, 0.007}, {0.622, 1.007, 0.007}, {0.689, 1.017, 0.010}, {0.779, 1.021, 0.005}, {0.853, 1.045, 0.007},
	{0.979, 1.045, 0.017}, {1.020, 1.054, 0.006}, {1.170, 1.059, 0.055}, {1.330, 1.057, 0.010}, {1.520, 1.074, 0.006},
	{1.740, 1.077, 0.004}, {1.830, 1.076, 0.008}, {2.070, 1.075, 0.006}, {2.500, 1.077, 0.005}, {2.690, 1.076, 0.009},
	{2.900, 1.079, 0.007}, {3.090, 1.075, 0.010}, {3.240, 1.074, 0.010}, {3.710, 1.068, 0.005}, {4.110, 1.063, 0.007},
	{5.010, 1.046, 0.005}, {5.850, 1.027, 0.007} 
};

double interpolate_arrington07_error( int FF = 1, double Q2_val = 4.5, bool normalized = true, bool uncorrect_TPE = false ){
	double interpolated_error = 0.0;
	double tau = Q2_val/(4*pow(0.938272, 2));
	double arr07_G_D = 1/pow( (1 + (Q2_val/0.71)), 2);
	double arr07_mu_proton = 2.79;

	if( FF == 1 ){ //GEp
		for( size_t val = 0; val < arr07_GEp_vals_and_errs.size() - 1; val++ ){
			double x1 = arr07_GEp_vals_and_errs[val][0];
			double x2 = arr07_GEp_vals_and_errs[val+1][0];

			double err1 = arr07_GEp_vals_and_errs[val][2];
			double err2 = arr07_GEp_vals_and_errs[val+1][2];
			if( (Q2_val >= x1) && (Q2_val < x2) ){
				interpolated_error = linear_interpolation( x1, err1, x2, err2, Q2_val );
				break;
			}
		}
	} 
	else{ //GMp
		for( size_t val = 0; val < arr07_GMp_vals_and_errs.size() - 1; val++ ){
			double x1 = arr07_GMp_vals_and_errs[val][0];
			double x2 = arr07_GMp_vals_and_errs[val+1][0];

			double err1 = arr07_GMp_vals_and_errs[val][2];
			double err2 = arr07_GMp_vals_and_errs[val+1][2];
			if( (Q2_val >= x1) && (Q2_val < x2) ){
				cout << "x1: " << x1 << ", x2: " << x2 << ", err1: " << err1 << ", err2: " << err2 << endl;
				interpolated_error = linear_interpolation( x1, err1, x2, err2, Q2_val );
				break;
			}
		}
	} 

	if( !normalized ){
		interpolated_error = interpolated_error*arr07_mu_proton*arr07_G_D;
	}

	if( uncorrect_TPE ){
		interpolated_error = 0.995*interpolated_error;
	}

	return interpolated_error;

}

double calc_arrington_2007_proton( int FF = 1, double Q2_val = 4.5, bool return_norm = true ){
	/// FF = 1 is GEp,
	/// FF = 2 is GMp,

	arr07_mu_p = 2.792847351;
	arr07_delta = 0.71;
	arr07_G_D = 1/pow( (1 + (Q2_val/(arr07_delta))), 2);

	tau = (Q2_val)/(4*arr07_proton_mass*arr07_proton_mass);

	if( FF == 1 ){ //GEp
		//not normalized by default

		numerator = 1 + (ae1*tau) + (ae2*tau*tau) + (ae3*pow(tau, 3));
		denominator = 1 + (be1*tau) + (be2*tau*tau) + (be3*pow(tau, 3)) + (be4*pow(tau, 4)) + (be5*pow(tau, 5));
		FF_calc_val = numerator/denominator;
	}
	else{ //GMp 
		//Normalized by arr07_mu_p by default
		
		numerator = 1 + (am1*tau) + (am2*tau*tau) + (am3*pow(tau, 3));
		denominator = 1 + (bm1*tau) + (bm2*tau*tau) + (bm3*pow(tau, 3)) + (bm4*pow(tau, 4)) + (bm5*pow(tau, 5));
		FF_calc_val = numerator/denominator;
	}

	if( return_norm ){
		return FF_calc_val/arr07_G_D;
	}
	else{
		if( FF == 1 ){
			return FF_calc_val;
		}
		
		else{
			return (arr07_mu_p)*FF_calc_val;
		}
	}

}

double calc_arr07_reduced_CS( double Q2_val, double epsilon, bool return_norm = false ){
	double reduced_CS = 0.0;

	double arr07_tau = (Q2_val)/(4*arr07_proton_mass*arr07_proton_mass);

	double arr07_GEp = calc_arrington_2007_proton( 1, Q2_val, return_norm );
	double arr07_GMp = 2.792847351*calc_arrington_2007_proton( 2, Q2_val, return_norm );

	reduced_CS = epsilon*pow(arr07_GEp, 2) + arr07_tau*pow(arr07_GMp, 2);

	return reduced_CS;
}

//--------------------------------
Double_t arr07_GEp_param( Double_t *x, Double_t *arr07_par_G ){
	double arr07_GEp = calc_arrington_2007_proton( 1, x[0] );
	return arr07_GEp;
}
Double_t arr07_GEp_error_high_fit( Double_t *x, Double_t *arr07_par_G_error ){
	double arr07_GEp = calc_arrington_2007_proton( 1, x[0] );
	double error = interpolate_arrington07_error( 1, x[0] );

	return arr07_GEp + error;
}
Double_t arr07_GEp_error_low_fit( Double_t *x, Double_t *arr07_par_G_error ){
	double arr07_GEp = calc_arrington_2007_proton( 1, x[0] );
	double error = interpolate_arrington07_error( 1, x[0] );

	return arr07_GEp - error;
}
//--------------------------------
Double_t arr07_GMp_param( Double_t *x, Double_t *arr07_par_G ){
	double arr07_GMp = calc_arrington_2007_proton( 2, x[0] );
	return arr07_GMp;
}
Double_t arr07_GMp_error_high_fit( Double_t *x, Double_t *arr07_par_G_error ){
	double arr07_GMp = calc_arrington_2007_proton( 2, x[0] );
	double error = interpolate_arrington07_error( 2, x[0] );

	return arr07_GMp + error;
}
Double_t arr07_GMp_error_low_fit( Double_t *x, Double_t *arr07_par_G_error ){
	double arr07_GMp = calc_arrington_2007_proton( 2, x[0] );
	double error = interpolate_arrington07_error( 2, x[0] );

	return arr07_GMp - error;
}
//--------------------------------


TF1 *tf1_arr07_GEp, *tf1_arr07_GEp_error_low_fit, *tf1_arr07_GEp_error_high_fit;
TF1 *tf1_arr07_GMp, *tf1_arr07_GMp_error_low_fit, *tf1_arr07_GMp_error_high_fit;

void plot_arr07( int FF = 1 ){
	double min_x = 0.001;
	double max_x = 12.0;
	double min_y = 0.2;
	double max_y = 1.2;
	double fit_Npx = 5000;

		TH1D *h_arr07_GEp = new TH1D("h_arr07_GEp", "GEp - Arr07 Parameterization", fit_Npx, min_x, max_x);
		h_arr07_GEp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_arr07_GEp->GetXaxis()->SetTitleOffset(1.25f);

		h_arr07_GEp->GetYaxis()->SetTitle("G_{E}^{p}/G_{D}");
		h_arr07_GEp->GetYaxis()->SetTitleOffset(1.25f);
		h_arr07_GEp->SetLineColor(0);
		h_arr07_GEp->GetYaxis()->SetRangeUser(min_y, max_y);
		h_arr07_GEp->GetXaxis()->SetRangeUser(min_x, max_x);
		h_arr07_GEp->SetStats(0);
		h_arr07_GEp->Draw();
		
		TCanvas *c_arr07_GEp = new TCanvas("c_arr07_GEp", "c_arr07_GEp", 600, 500);
		c_arr07_GEp->SetLogx();
		c_arr07_GEp->Range(min_x, min_y, max_x, max_y);

		tf1_arr07_GEp = new TF1("tf1_arr07_GEp", arr07_GEp_param, min_x, max_x);
		tf1_arr07_GEp->SetNpx(fit_Npx);
		tf1_arr07_GEp->SetLineColor(1);
		tf1_arr07_GEp->SetLineWidth(2);
		tf1_arr07_GEp->Draw("same");

		tf1_arr07_GEp_error_low_fit = new TF1("tf1_arr07_GEp_error_low_fit", arr07_GEp_error_low_fit, min_x, max_x);
		tf1_arr07_GEp_error_low_fit->SetNpx(fit_Npx);
		tf1_arr07_GEp_error_low_fit->SetLineColor(1);
		tf1_arr07_GEp_error_low_fit->SetLineWidth(1);
		tf1_arr07_GEp_error_low_fit->SetLineStyle(7);
		tf1_arr07_GEp_error_low_fit->Draw("same");

		tf1_arr07_GEp_error_high_fit = new TF1("tf1_arr07_GEp_error_high_fit", arr07_GEp_error_high_fit, min_x, max_x);
		tf1_arr07_GEp_error_high_fit->SetNpx(fit_Npx);
		tf1_arr07_GEp_error_high_fit->SetLineColor(1);
		tf1_arr07_GEp_error_high_fit->SetLineWidth(1);
		tf1_arr07_GEp_error_high_fit->SetLineStyle(7);
		tf1_arr07_GEp_error_high_fit->Draw("same");

		fshade_btw_TF1s_w_3rd(c_arr07_GEp, h_arr07_GEp, tf1_arr07_GEp_error_high_fit, tf1_arr07_GEp_error_low_fit, tf1_arr07_GEp, min_y, max_y);



}



#endif
