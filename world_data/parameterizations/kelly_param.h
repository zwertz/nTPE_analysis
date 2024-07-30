#ifndef KELLY_PARAM_H
#define KELLY_PARAM_H

#include "/work/halla/sbs/jboyd/include/utility_functions.h"

double par_G[4], par_G_error[4];

double par_GEp[5] = {1.0, -0.24, 10.98, 12.82, 21.97};
double par_GEp_error[5] = {0.0, 0.12, 0.19, 1.1, 6.8};

double par_GMp[5] = {1.0, 0.12, 10.97, 18.86, 6.55};
double par_GMp_error[5] = {0.0, 0.04, 0.11, 0.28, 1.2};

double par_GMn[5] = {1.0, 2.33, 14.72, 24.20, 84.1};
double par_GMn_error[5] = {0.0, 1.4, 1.7, 9.8, 41};

double par_GEn[2] = {1.70, 3.30};
double par_GEn_error[2] = {0.04, 0.32};

TCanvas *c_kelly_GMn;

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
		return G_calc;		
	}
	if( normed ){
		return G_calc/G_D;				
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

		division = numerator/denominator;
		division_error = CalculateErrorMultiplicationDivision( numerator, numerator_error, denominator, denominator_error, division);
		division_error = division_error/division;		
	}

	if( normed ){
		division_error = division_error/G_D;
	}

	// cout << "division error: " << division_error << endl;
	return division_error;
}

//-------------------------------------
//----------------- GEp ------------
Double_t kelly_GEp_parameterization( Double_t *x, Double_t *par_G ){
	double GEp = calc_kelly_parameterization( 1, x[0], true );

	return GEp;
}
Double_t kelly_GEp_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GEp = calc_kelly_parameterization( 1, x[0], true );
	double error = calc_kelly_errors( 1, x[0], true );

	return GEp + GEp*fabs(error);
}
Double_t kelly_GEp_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GEp = calc_kelly_parameterization( 1, x[0], true );
	double error = calc_kelly_errors( 1, x[0], true );

	return GEp - GEp*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GMp ------------
Double_t kelly_GMp_parameterization( Double_t *x, Double_t *par_G ){
	double GMp = calc_kelly_parameterization( 2, x[0], true );

	return GMp;
}
Double_t kelly_GMp_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GMp = calc_kelly_parameterization( 2, x[0], true );
	double error = calc_kelly_errors( 2, x[0], false );

	return GMp + x[0]*GMp*fabs(error);
}
Double_t kelly_GMp_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GMp = calc_kelly_parameterization( 2, x[0], true );
	double error = calc_kelly_errors( 2, x[0], false );

	return GMp - x[0]*GMp*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GEn ------------
Double_t kelly_GEn_parameterization( Double_t *x, Double_t *par_G ){
	double GEn = calc_kelly_parameterization( 3, x[0], true );

	return GEn;
}
Double_t kelly_GEn_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GEn = calc_kelly_parameterization( 3, x[0], true );
	double error = calc_kelly_errors( 3, x[0], false );

	return GEn + GEn*fabs(error);
}
Double_t kelly_GEn_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GEn = calc_kelly_parameterization( 3, x[0], true );
	double error = calc_kelly_errors( 3, x[0], false );

	return GEn - GEn*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GMn ------------
Double_t kelly_GMn_parameterization( Double_t *x, Double_t *par_G ){
	double GMn = calc_kelly_parameterization( 4, x[0], true );

	return GMn;
}
Double_t kelly_GMn_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GMn = calc_kelly_parameterization( 4, x[0], true );
	double error = sqrt(x[0])*GMn*calc_kelly_errors( 4, x[0], false );
	
	return GMn + fabs(error);
}
Double_t kelly_GMn_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GMn = calc_kelly_parameterization( 4, x[0], true );
	double error = sqrt(x[0])*GMn*calc_kelly_errors( 4, x[0], false );
	
	return GMn - fabs(error);
}
//-------------------------------------

TF1 *tf1_kelly_GEp, *tf1_kelly_GEp_error_high_fit, *tf1_kelly_GEp_error_low_fit;
TF1 *tf1_kelly_GMp, *tf1_kelly_GMp_error_high_fit, *tf1_kelly_GMp_error_low_fit;
TF1 *tf1_kelly_GEn, *tf1_kelly_GEn_error_high_fit, *tf1_kelly_GEn_error_low_fit;
TF1 *tf1_kelly_GMn, *tf1_kelly_GMn_error_high_fit, *tf1_kelly_GMn_error_low_fit;

//-------------------------------
//FF 0: All, 1: GEp, 2: GMp, 3: GEn, 4: GMn
void plot_kelly_FFs( int FF = 0, bool plot_jboyd = false, bool plot_jboyd_GMn4 = false, bool plot_jboyd_approx = false ){

	if( FF == 0 || FF == 1 ){
		double GEp_min_x = 0.003;
		double GEp_max_x = 10;
		double GEp_min_y = -0.1;
		double GEp_max_y = 1.2;
		double GEp_fit_Npx = 5000;

		TCanvas *c_kelly_GEp = new TCanvas("c_kelly_GEp", "c_kelly_GEp", 600, 500);
		c_kelly_GEp->SetGrid();
		TH1D *h_GEp = new TH1D("h_GEp", "GEp - Kelly Parameterization", GEp_fit_Npx, GEp_min_x, GEp_max_x);
		h_GEp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEp->GetXaxis()->SetTitleOffset(1.25f);

		h_GEp->GetYaxis()->SetTitle("G_{E}^{p}/G_{D}");
		h_GEp->GetYaxis()->SetTitleOffset(1.25f);
		h_GEp->SetLineColor(0);
		h_GEp->GetYaxis()->SetRangeUser(GEp_min_y, GEp_max_y);
		h_GEp->GetXaxis()->SetRangeUser(GEp_min_x, GEp_max_x);
		h_GEp->SetStats(0);
		h_GEp->Draw();

		c_kelly_GEp->SetLogx();
		c_kelly_GEp->Range(GEp_min_x, GEp_min_y, GEp_max_x, GEp_max_y);
		tf1_kelly_GEp = new TF1("tf1_kelly_GEp", kelly_GEp_parameterization, GEp_min_x, GEp_max_x, 0);
		tf1_kelly_GEp->SetNpx(GEp_fit_Npx);
		tf1_kelly_GEp->SetLineColor(1);
		tf1_kelly_GEp->SetLineWidth(2);
		tf1_kelly_GEp->Draw("same");

		tf1_kelly_GEp_error_high_fit = new TF1("tf1_kelly_GEp_error_high_fit", kelly_GEp_error_high_fit, GEp_min_x, GEp_max_x, 0);
		tf1_kelly_GEp_error_high_fit->SetNpx(GEp_fit_Npx);
		tf1_kelly_GEp_error_high_fit->SetLineColor(1);
		tf1_kelly_GEp_error_high_fit->SetLineWidth(1);
		tf1_kelly_GEp_error_high_fit->SetLineStyle(7);
		tf1_kelly_GEp_error_high_fit->Draw("same");

		tf1_kelly_GEp_error_low_fit = new TF1("tf1_kelly_GEp_error_low", kelly_GEp_error_low_fit, GEp_min_x, GEp_max_x, 0);
		tf1_kelly_GEp_error_low_fit->SetNpx(GEp_fit_Npx);
		tf1_kelly_GEp_error_low_fit->SetLineColor(1);
		tf1_kelly_GEp_error_low_fit->SetLineWidth(1);
		tf1_kelly_GEp_error_low_fit->SetLineStyle(7);
		tf1_kelly_GEp_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd(c_kelly_GEp, h_GEp, tf1_kelly_GEp_error_high_fit, tf1_kelly_GEp_error_low_fit, tf1_kelly_GEp, GEp_min_y, GEp_max_y);
	
	}

	if( FF == 0 || FF == 2 ){
		double GMp_min_x = 0.01;
		double GMp_max_x = 60;
		double GMp_min_y = 0.6;
		double GMp_max_y = 1.15;
		double GMp_fit_Npx = 5000;

		TCanvas *c_kelly_GMp = new TCanvas("c_kelly_GMp", "c_kelly_GMp", 600, 500);
		c_kelly_GMp->SetGrid();
		TH1D *h_GMp = new TH1D("h_GMp", "GMp - Kelly Parameterization", GMp_fit_Npx, GMp_min_x, GMp_max_x);
		h_GMp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMp->GetXaxis()->SetTitleOffset(1.25f);

		h_GMp->GetYaxis()->SetTitle("G_{M}^{p}/(#mu_{n} G_{D})");
		h_GMp->GetYaxis()->SetTitleOffset(1.25f);
		h_GMp->SetLineColor(0);
		h_GMp->GetYaxis()->SetRangeUser(GMp_min_y, GMp_max_y);
		h_GMp->GetXaxis()->SetRangeUser(GMp_min_x, GMp_max_x);
		h_GMp->SetStats(0);
		h_GMp->Draw();

		c_kelly_GMp->SetLogx();
		c_kelly_GMp->Range(GMp_min_x, GMp_min_y, GMp_max_x, GMp_max_y);
		tf1_kelly_GMp = new TF1("tf1_kelly_GMp", kelly_GMp_parameterization, GMp_min_x, GMp_max_x, 0);
		tf1_kelly_GMp->SetNpx(GMp_fit_Npx);
		tf1_kelly_GMp->SetLineColor(1);
		tf1_kelly_GMp->SetLineWidth(2);
		tf1_kelly_GMp->Draw("same");

		tf1_kelly_GMp_error_high_fit = new TF1("tf1_kelly_GMp_error_high_fit", kelly_GMp_error_high_fit, GMp_min_x, GMp_max_x, 0);
		tf1_kelly_GMp_error_high_fit->SetNpx(GMp_fit_Npx);
		tf1_kelly_GMp_error_high_fit->SetLineColor(1);
		tf1_kelly_GMp_error_high_fit->SetLineWidth(1);
		tf1_kelly_GMp_error_high_fit->SetLineStyle(7);
		tf1_kelly_GMp_error_high_fit->Draw("same");

		tf1_kelly_GMp_error_low_fit = new TF1("tf1_kelly_GMp_error_low", kelly_GMp_error_low_fit, GMp_min_x, GMp_max_x, 0);
		tf1_kelly_GMp_error_low_fit->SetNpx(GMp_fit_Npx);
		tf1_kelly_GMp_error_low_fit->SetLineColor(1);
		tf1_kelly_GMp_error_low_fit->SetLineWidth(1);
		tf1_kelly_GMp_error_low_fit->SetLineStyle(7);
		tf1_kelly_GMp_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd(c_kelly_GMp, h_GMp, tf1_kelly_GMp_error_high_fit, tf1_kelly_GMp_error_low_fit, tf1_kelly_GMp, GMp_min_y, GMp_max_y);
	
	}

	if( FF == 0 || FF == 3 ){
		double GEn_min_x = 0.003;
		double GEn_max_x = 10;
		double GEn_min_y = -0.1;
		double GEn_max_y = 0.55;
		double GEn_fit_Npx = 5000;

		TCanvas *c_kelly_GEn = new TCanvas("c_kelly_GEn", "c_kelly_GEn", 600, 500);
		c_kelly_GEn->SetGrid();
		TH1D *h_GEn = new TH1D("h_GEn", "GEn - Kelly Parameterization", GEn_fit_Npx, GEn_min_x, GEn_max_x);
		h_GEn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEn->GetXaxis()->SetTitleOffset(1.25f);

		h_GEn->GetYaxis()->SetTitle("G_{E}^{n}/G_{D}");
		h_GEn->GetYaxis()->SetTitleOffset(1.25f);
		h_GEn->SetLineColor(0);
		h_GEn->GetYaxis()->SetRangeUser(GEn_min_y, GEn_max_y);
		h_GEn->GetXaxis()->SetRangeUser(GEn_min_x, GEn_max_x);
		h_GEn->SetStats(0);
		h_GEn->Draw();

		c_kelly_GEn->SetLogx();
		c_kelly_GEn->Range(GEn_min_x, GEn_min_y, GEn_max_x, GEn_max_y);
		tf1_kelly_GEn = new TF1("tf1_kelly_GEn", kelly_GEn_parameterization, GEn_min_x, GEn_max_x, 0);
		tf1_kelly_GEn->SetNpx(GEn_fit_Npx);
		tf1_kelly_GEn->SetLineColor(1);
		tf1_kelly_GEn->SetLineWidth(2);
		tf1_kelly_GEn->Draw("same");

		tf1_kelly_GEn_error_high_fit = new TF1("tf1_kelly_GEn_error_high_fit", kelly_GEn_error_high_fit, GEn_min_x, GEn_max_x, 0);
		tf1_kelly_GEn_error_high_fit->SetNpx(GEn_fit_Npx);
		tf1_kelly_GEn_error_high_fit->SetLineColor(1);
		tf1_kelly_GEn_error_high_fit->SetLineWidth(1);
		tf1_kelly_GEn_error_high_fit->SetLineStyle(7);
		tf1_kelly_GEn_error_high_fit->Draw("same");

		tf1_kelly_GEn_error_low_fit = new TF1("tf1_kelly_GEn_error_low", kelly_GEn_error_low_fit, GEn_min_x, GEn_max_x, 0);
		tf1_kelly_GEn_error_low_fit->SetNpx(GEn_fit_Npx);
		tf1_kelly_GEn_error_low_fit->SetLineColor(1);
		tf1_kelly_GEn_error_low_fit->SetLineWidth(1);
		tf1_kelly_GEn_error_low_fit->SetLineStyle(7);
		tf1_kelly_GEn_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd(c_kelly_GEn, h_GEn, tf1_kelly_GEn_error_high_fit, tf1_kelly_GEn_error_low_fit, tf1_kelly_GEn, GEn_min_y, GEn_max_y);
	
	}

	if( FF == 0 || FF == 4 ){
		double GMn_min_x = 0.01;
		double GMn_max_x = 60;
		double GMn_min_y = 0.4;
		double GMn_max_y = 1.15;
		double GMn_fit_Npx = 5000;

		c_kelly_GMn = new TCanvas("c_kelly_GMn", "c_kelly_GMn", 600, 500);
		c_kelly_GMn->SetGrid();
		TH1D *h_GMn = new TH1D("h_GMn", "GMn - Kelly Parameterization", GMn_fit_Npx, GMn_min_x, GMn_max_x);
		h_GMn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMn->GetXaxis()->SetTitleOffset(1.25f);

		h_GMn->GetYaxis()->SetTitle("G_{M}^{n}/(#mu_{n} G_{D})");
		h_GMn->GetYaxis()->SetTitleOffset(1.25f);
		h_GMn->SetLineColor(0);
		h_GMn->GetYaxis()->SetRangeUser(GMn_min_y, GMn_max_y);
		h_GMn->GetXaxis()->SetRangeUser(GMn_min_x, GMn_max_x);
		h_GMn->SetStats(0);
		h_GMn->Draw();

		c_kelly_GMn->SetLogx();
		c_kelly_GMn->Range(GMn_min_x, GMn_min_y, GMn_max_x, GMn_max_y);
		tf1_kelly_GMn = new TF1("tf1_kelly_GMn", kelly_GMn_parameterization, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn->SetLineColor(1);
		tf1_kelly_GMn->SetLineWidth(2);
		tf1_kelly_GMn->Draw("same");

		tf1_kelly_GMn_error_high_fit = new TF1("tf1_kelly_GMn_error_high_fit", kelly_GMn_error_high_fit, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn_error_high_fit->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn_error_high_fit->SetLineColor(1);
		tf1_kelly_GMn_error_high_fit->SetLineWidth(1);
		tf1_kelly_GMn_error_high_fit->SetLineStyle(7);
		tf1_kelly_GMn_error_high_fit->Draw("same");

		tf1_kelly_GMn_error_low_fit = new TF1("tf1_kelly_GMn_error_low", kelly_GMn_error_low_fit, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn_error_low_fit->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn_error_low_fit->SetLineColor(1);
		tf1_kelly_GMn_error_low_fit->SetLineWidth(1);
		tf1_kelly_GMn_error_low_fit->SetLineStyle(7);
		tf1_kelly_GMn_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd(c_kelly_GMn, h_GMn, tf1_kelly_GMn_error_high_fit, tf1_kelly_GMn_error_low_fit, tf1_kelly_GMn, GMn_min_y, GMn_max_y);
	
	}

	if( plot_jboyd ){
		if( plot_jboyd_GMn4 ){
			plot_jboyd_data( true );
		}
		else{
			plot_jboyd_data();		
		}


		TString drawStyle = "*same";

		if( !plot_jboyd_approx ){
			tge_GMn_SBS8_jboyd_stat->Draw("*same");
			tge_GMn_SBS9_jboyd_stat->Draw(drawStyle.Data());
			tge_GMn_SBS8_jboyd->Draw(drawStyle.Data());
			tge_GMn_SBS9_jboyd->Draw(drawStyle.Data());

			if( plot_jboyd_GMn4 ){
				tge_GMn_SBS4_jboyd->Draw(drawStyle.Data());				
			}

			TLegend *tleg_GMn_jboyd;
			tleg_GMn_jboyd = new TLegend(0.125, 0.20, 0.80, 0.35);
			tleg_GMn_jboyd->SetTextSize(0.03);

			if( plot_jboyd_GMn4 ){
				tleg_GMn_jboyd->AddEntry((TObject*)0, "", "");
				tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 G_{M}^{n} = %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[2], jboyd_Q2_vec[2]));
				tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 syst. err: %.4f", jboyd_GMn_syst_error_vec[2]));
				tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd_stat, Form("SBS4 stat. err: %.4f", jboyd_GMn_stat_error_vec[2] - jboyd_GMn_syst_error_vec[2]));
				tleg_GMn_jboyd->AddEntry((TObject*)0, "", "");
			}	

			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd, Form("SBS8 G_{M}^{n} = %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[0], jboyd_Q2_vec[0]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd, Form("SBS8 syst. err: %.4f", jboyd_GMn_syst_error_vec[0]));	
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd_stat, Form("SBS8 stat. err: %.4f", jboyd_GMn_stat_error_vec[0] - jboyd_GMn_syst_error_vec[0]));	
			tleg_GMn_jboyd->AddEntry((TObject*)0, "", "");

			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd, Form("SBS9 G_{M}^{n} = %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[1], jboyd_Q2_vec[1]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd, Form("SBS9 syst. err: %.4f ", jboyd_GMn_syst_error_vec[1]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd_stat, Form("SBS9 stat. err: %.4f", jboyd_GMn_stat_error_vec[1] - jboyd_GMn_syst_error_vec[1]));

			tleg_GMn_jboyd->Draw(drawStyle.Data());			
		}

		if( plot_jboyd_approx ){
			tge_GMn_SBS8_approx_jboyd->Draw(drawStyle.Data());
			tge_GMn_SBS9_approx_jboyd->Draw(drawStyle.Data());

			TLegend *tleg_GMn_approx_jboyd;
			tleg_GMn_approx_jboyd = new TLegend(0.125, 0.20, 0.55, 0.35);
			tleg_GMn_approx_jboyd->SetTextSize(0.03);
			tleg_GMn_approx_jboyd->AddEntry(tge_GMn_SBS8_approx_jboyd, Form("SBS8 G_{M}^{n} ~ %.5f +/- %.5f", jboyd_GMn_approx_vec[0], jboyd_GMn_approx_error_vec[0] ));
			tleg_GMn_approx_jboyd->AddEntry(tge_GMn_SBS9_approx_jboyd, Form("SBS9 G_{M}^{n} ~ %.5f +/- %.5f", jboyd_GMn_approx_vec[1], jboyd_GMn_approx_error_vec[1] ));
			tleg_GMn_approx_jboyd->Draw(drawStyle.Data());			
		
		}

		c_kelly_GMn->cd();
		tge_GMn_SBS4_jboyd_stat->Draw("*same");
		tge_GMn_SBS4_jboyd->Draw("*same");
		tge_GMn_SBS4_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS4_jboyd->Draw("same");
		tg_GMn_SBS4_jboyd->Draw("*same");
		tg_GMn_SBS4_jboyd->SetLineColor(1);
		tg_GMn_SBS4_jboyd->SetMarkerColor(1);
		tg_GMn_SBS4_jboyd->SetMarkerStyle(kOpenCircle);
		tg_GMn_SBS4_jboyd->Draw("same");

		tge_GMn_SBS8_jboyd_stat->Draw("*same");
		tge_GMn_SBS8_jboyd->Draw("*same");
		tge_GMn_SBS8_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS8_jboyd->Draw("same");
		tg_GMn_SBS8_jboyd->Draw("*same");
		tg_GMn_SBS8_jboyd->SetLineColor(1);
		tg_GMn_SBS8_jboyd->SetMarkerColor(1);
		tg_GMn_SBS8_jboyd->SetMarkerStyle(kOpenCircle);
		tg_GMn_SBS8_jboyd->Draw("same");

		tge_GMn_SBS9_jboyd_stat->Draw("*same");
		tge_GMn_SBS9_jboyd->Draw("*same");
		tge_GMn_SBS9_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS9_jboyd->Draw("same");
		tg_GMn_SBS9_jboyd->Draw("*same");
		tg_GMn_SBS9_jboyd->SetLineColor(1);
		tg_GMn_SBS9_jboyd->SetMarkerColor(1);
		tg_GMn_SBS9_jboyd->SetMarkerStyle(kOpenCircle);
		tg_GMn_SBS9_jboyd->Draw("same");

	}

}

void general_plot_kelly(TCanvas *c_to_plot_on, TLegend *tleg_to_add = NULL, int setLineColor = 2, int setLineStyle = 9, int setLineWidth = 2, bool plot_with_error = false){
	
	double GMn_min_x = 0.01;
	double GMn_max_x = 60;
	double GMn_min_y = 0.4;
	double GMn_max_y = 1.15;
	double GMn_fit_Npx = 5000;

	c_to_plot_on->cd();
	c_to_plot_on->SetGrid();
	TH1D *h_GMn = new TH1D("h_GMn", "GMn - Kelly Parameterization", GMn_fit_Npx, GMn_min_x, GMn_max_x);
	// h_GMn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
	// h_GMn->GetXaxis()->SetTitleOffset(1.25f);

	// h_GMn->GetYaxis()->SetTitle("G_{M}^{n}/(#mu_{n} G_{D})");
	// h_GMn->GetYaxis()->SetTitleOffset(1.25f);
	// h_GMn->SetLineColor(0);
	// h_GMn->GetYaxis()->SetRangeUser(GMn_min_y, GMn_max_y);
	// h_GMn->GetXaxis()->SetRangeUser(GMn_min_x, GMn_max_x);
	// h_GMn->SetStats(0);
	h_GMn->Draw("same");

	// c_to_plot_on->Range(GMn_min_x, GMn_min_y, GMn_max_x, GMn_max_y);
	tf1_kelly_GMn = new TF1("tf1_kelly_GMn", kelly_GMn_parameterization, GMn_min_x, GMn_max_x, 0);
	tf1_kelly_GMn->SetNpx(GMn_fit_Npx);
	tf1_kelly_GMn->SetLineColor(setLineColor);
	tf1_kelly_GMn->SetLineStyle(setLineStyle);
	tf1_kelly_GMn->SetLineWidth(setLineWidth);
	tf1_kelly_GMn->Draw("same");

	if( plot_with_error ){
		tf1_kelly_GMn_error_high_fit = new TF1("tf1_kelly_GMn_error_high_fit", kelly_GMn_error_high_fit, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn_error_high_fit->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn_error_high_fit->SetLineColor(1);
		tf1_kelly_GMn_error_high_fit->SetLineWidth(1);
		tf1_kelly_GMn_error_high_fit->SetLineStyle(7);
		tf1_kelly_GMn_error_high_fit->Draw("same");

		tf1_kelly_GMn_error_low_fit = new TF1("tf1_kelly_GMn_error_low", kelly_GMn_error_low_fit, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn_error_low_fit->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn_error_low_fit->SetLineColor(1);
		tf1_kelly_GMn_error_low_fit->SetLineWidth(1);
		tf1_kelly_GMn_error_low_fit->SetLineStyle(7);
		tf1_kelly_GMn_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd(c_to_plot_on, h_GMn, tf1_kelly_GMn_error_high_fit, tf1_kelly_GMn_error_low_fit, tf1_kelly_GMn, GMn_min_y, GMn_max_y);	
	}
	if( tleg_to_add != NULL ){
		TString label_gmn_nominal = "Kelly";
		if( doesLegendEntryExist( tleg_to_add, label_gmn_nominal) ){

		}
		else{
			tleg_to_add->AddEntry(tf1_kelly_GMn, label_gmn_nominal.Data());
			c_to_plot_on->Update();
			c_to_plot_on->Modified();			
		}

	}
	// if( tleg_to_add != NULL ){
	// 	tleg_to_add->AddEntry(tf1_kelly_GMn, "Kelly");
	// 	c_to_plot_on->Update();
	// 	c_to_plot_on->Modified();
	// }


}

#endif