#ifndef BRADFORD_PARAM_H
#define BRADFORD_PARAM_H

#include "/work/halla/sbs/jboyd/include/experimental_constants.h"

vector<double> bradford_gep_a = {1.00, -0.06, 0.00};
vector<double> bradford_gep_b = {11.10, 13.60, 33.0, 0.00};

vector<double> bradford_gmp_a = {1.00, 0.15, 0.00};
vector<double> bradford_gmp_b = {11.10, 19.60, 7.54, 0.00};

vector<double> bradford_gen_a = {0.00, 1.25, 1.30};
vector<double> bradford_gen_b = {-9.86, 305.00, -758.00, 802.00};

vector<double> bradford_gmn_a = {1.00, 1.81, 0.00};
vector<double> bradford_gmn_b = {14.10, 30.70, 68.70, 0.00};

TF1 *tf1_bradford_GMn;

Double_t calc_bradford_parameterization( int FF, double q2, bool normed = true ){
	int n = 1;
	double param_MN = 0.0;
	double param_mu = 0.0;

	if( FF == 1 || FF == 2 ){
		param_MN = Mp;
		param_mu = mu_p;
	}
	if( FF == 3 || FF == 4 ){
		param_MN = Mn;
		param_mu = mu_n;
	}

	double tau = q2/(4*MN*MN);
	double delta2 = 0.71;

	double G_D = pow( (1 + (q2/delta2)), -2);

	double G_calc = 0.0;

	double a_sum = 0.0, b_sum = 0.0;

//----GEp
	if( FF == 1 ){
		for( size_t k = 0; k < bradford_gep_a.size(); k++ ){
			a_sum += bradford_gep_a[k]*pow( tau, k );
		}
		for( size_t k = 0; k < bradford_gep_b.size(); k++ ){
			b_sum += bradford_gep_b[k]*pow( tau, k + 1 );
		}
		b_sum += 1.0;
	}

//----GMp
	if( FF == 2 ){
		for( size_t k = 0; k < bradford_gmp_a.size(); k++ ){
			a_sum += bradford_gmp_a[k]*pow( tau, k );
		}
		for( size_t k = 0; k < bradford_gmp_b.size(); k++ ){
			b_sum += bradford_gmp_b[k]*pow( tau, k + 1 );
		}
		b_sum += 1.0;
	}

//----GEn
	if( FF == 3 ){
		for( size_t k = 0; k < bradford_gen_a.size(); k++ ){
			a_sum += bradford_gen_a[k]*pow( tau, k );
		}
		for( size_t k = 0; k < bradford_gen_b.size(); k++ ){
			b_sum += bradford_gen_b[k]*pow( tau, k + 1 );
		}
		b_sum += 1.0;
	}

//----GMn
	if( FF == 4 ){
		for( size_t k = 0; k < bradford_gmn_a.size(); k++ ){
			a_sum += bradford_gmn_a[k]*pow( tau, k );
		}
		for( size_t k = 0; k < bradford_gmn_b.size(); k++ ){
			b_sum += bradford_gmn_b[k]*pow( tau, k + 1 );
		}
		b_sum += 1.0;
	}

	G_calc = a_sum/b_sum;

	if( normed ){
		if( FF == 1 || FF == 3 ){
			G_calc = G_calc/G_D;
		}
		if( FF == 2 || FF == 4 ){
			G_calc = G_calc/(G_D);
		}
	}

	return G_calc;

}

Double_t bradford_GMn_parameterization( Double_t *x, Double_t *par_G ){
	double GMn = calc_bradford_parameterization( 4, x[0], true );
	return GMn;
}

void general_plot_bradford(TCanvas *c_to_plot_on, TLegend *tleg_to_add = NULL, int setLineColor = 2, int setLineStyle = 9, int setLineWidth = 2, bool plot_with_error = false){
	double GMn_min_x = 0.01;
	double GMn_max_x = 60;
	double GMn_min_y = 0.4;
	double GMn_max_y = 1.15;
	double GMn_fit_Npx = 5000;

	c_to_plot_on->cd();
	c_to_plot_on->SetGrid();
	TH1D *h_GMn = new TH1D("h_GMn", "GMn - Bradford Parameterization", GMn_fit_Npx, GMn_min_x, GMn_max_x);
	h_GMn->Draw("same");

	tf1_bradford_GMn = new TF1("tf1_kelly_GMn", bradford_GMn_parameterization, GMn_min_x, GMn_max_x, 0);
	tf1_bradford_GMn->SetNpx(GMn_fit_Npx);
	tf1_bradford_GMn->SetLineColor(setLineColor);
	tf1_bradford_GMn->SetLineStyle(setLineStyle);
	tf1_bradford_GMn->SetLineWidth(setLineWidth);
	tf1_bradford_GMn->Draw("same");

	if( tleg_to_add != NULL ){
		tleg_to_add->AddEntry(tf1_kelly_GMn, "Kelly");
		c_to_plot_on->Update();
		c_to_plot_on->Modified();
	}
}

#endif
