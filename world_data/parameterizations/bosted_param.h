#ifndef BOSTED_PARAM_H
#define BOSTED_PARAM_H

const int num_pars = 5;

double par_G_bosted[num_pars], par_G_bosted_error[num_pars];

double par_GEp_bosted[num_pars] = {0.14, 3.01, 0.02, 1.20, 0.32};
double par_GEp_bosted_error[num_pars] = {0.0, 0.0, 0.0, 0.0, 0.0};

double par_GMp_bosted[num_pars] = {0.35, 2.44, 0.50, 1.04, 0.34};
double par_GMp_bosted_error[num_pars] = {0.0, 0.0, 0.0, 0.0, 0.0};

double par_GMn_bosted[num_pars] = {-1.74, 9.29, -7.63, 4.63, 0.0};
double par_GMn_bosted_error[num_pars] = {0.0, 0.0, 0.0, 0.0, 0.0};

double par_GEn_bosted[num_pars] = {1.25, 18.3, 0.0, 0.0, 0.0};
double par_GEn_bosted_error[num_pars] = {0.0, 0.0, 0.0, 0.0, 0.0};

Double_t calc_bosted_parameterization( int FF, double q2, bool normed = true ){
	double q = sqrt(q2);
	double mu_n = 1.9130427;
	double mu_p = 2.792847351;
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

	double par_G_select[num_pars];


//----GEp
	if( FF == 1 ){
		for( int i = 0; i < num_pars; i++ ){
			par_G_select[i] = par_GEp_bosted[i];
		}
	}
//----GMp
	if( FF == 2 ){
		for( int i = 0; i < num_pars; i++ ){
			par_G_select[i] = par_GMp_bosted[i];
		}
	}
//----GEn
	if( FF == 3 ){
		for( int i = 0; i < 2; i++ ){
			par_G_select[i] = par_GEn_bosted[i];
		}
		for( int i = 2; i < num_pars; i++ ){
			par_G_select[i] = 0.0;
		}
	}
//----GMn
	if( FF == 4 ){
		for( int i = 0; i < num_pars; i++ ){
			par_G_select[i] = par_GMn_bosted[i];
		}
	}

//----Calculation
	a_sum = 1.0;

	for( int i = 0; i < num_pars; i++ ){
		b_sum += par_G_select[i]*pow( q, i+1);
	}
	b_sum += 1.0;

	//--- For everything but GEn
	if( FF != 3){
		G_calc = a_sum/b_sum;	
	}

	//--- For GEn
	if( FF == 3 ){
		G_calc = (-par_GEn_bosted[0]*mu_n*tau*G_D)/(1 + par_GEn_bosted[1]*tau); 
	}


//-----Return normalized value or not
	if( !normed ){
		return G_calc;		
	}
	if( normed ){
		return G_calc/G_D;				
	}

}

Double_t calc_bosted_errors( int FF, double q2, bool normed = true ){
	
	double q = sqrt(q2);

	double delta2 = 0.71;
	double G_D = pow( (1 + (q2/delta2)), -2);

	double mu_n = 1.9130427;
	double mu_p = 2.792847351;

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

	double par_G_select[num_pars], par_G_select_error[num_pars];

//----GEp
	if( FF == 1 ){
		for( int i = 0; i < num_pars; i++ ){
			par_G_select[i] = par_GEp_bosted[i];
			par_G_select_error[i] = par_GEp_bosted_error[i];
		}
	}
//----GMp
	if( FF == 2 ){
		for( int i = 0; i < num_pars; i++ ){
			par_G_select[i] = par_GMp_bosted[i];
			par_G_select_error[i] = par_GMp_bosted_error[i];
		}
	}
//----GEn
	if( FF == 3 ){
		for( int i = 0; i < 2; i++ ){
			par_G_select[i] = par_GEn_bosted[i];
			par_G_select_error[i] = par_GEn_bosted_error[i];
		}
		for( int i = 2; i < num_pars; i++ ){
			par_G_select[i] = 0.0;
			par_G_select_error[i] = 0.0;
		}
	}
//----GMn
	if( FF == 4 ){
		for( int i = 0; i < num_pars; i++ ){
			par_G_select[i] = par_GMn_bosted[i];
			par_G_select_error[i] = par_GMn_bosted_error[i];
		}
	}

	if( FF != 3 ){
		numerator = 1.0;
		numerator_error = 0.0;

		for( int i = 0; i < num_pars; i++ ){
			denominator += par_G_select[i]*pow(q , i+1 );
			denominator_error += pow(par_G_select_error[i], 2);
		}

		denominator += 1.0;
		denominator_error = sqrt(denominator_error);

		division = numerator/denominator;
		division_error = CalculateErrorMultiplicationDivision( numerator, numerator_error, denominator, denominator_error, division);
	}

	if( FF == 3 ){
		numerator = par_GEn_bosted[0]*mu_n*tau*G_D;
		numerator_error = par_GEn_bosted_error[0];

		denominator = 1+ par_GEn_bosted[1]*tau;
		denominator_error = par_GEn_bosted_error[1];

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
Double_t bosted_GEp_parameterization( Double_t *x, Double_t *par_G ){
	double GEp = calc_bosted_parameterization( 1, x[0], true );

	return GEp;
}
Double_t bosted_GEp_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GEp = calc_bosted_parameterization( 1, x[0], true );
	double error = calc_bosted_errors( 1, x[0], true );

	return GEp + GEp*fabs(error);
}
Double_t bosted_GEp_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GEp = calc_bosted_parameterization( 1, x[0], true );
	double error = calc_bosted_errors( 1, x[0], true );

	return GEp - GEp*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GMp ------------
Double_t bosted_GMp_parameterization( Double_t *x, Double_t *par_G ){
	double GMp = calc_bosted_parameterization( 2, x[0], true );

	return GMp;
}
Double_t bosted_GMp_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GMp = calc_bosted_parameterization( 2, x[0], true );
	double error = calc_bosted_errors( 2, x[0], false );

	return GMp + x[0]*GMp*fabs(error);
}
Double_t bosted_GMp_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GMp = calc_bosted_parameterization( 2, x[0], true );
	double error = calc_bosted_errors( 2, x[0], false );

	return GMp - x[0]*GMp*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GEn ------------
Double_t bosted_GEn_parameterization( Double_t *x, Double_t *par_G ){
	double GEn = calc_bosted_parameterization( 3, x[0], true );

	return GEn;
}
Double_t bosted_GEn_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GEn = calc_bosted_parameterization( 3, x[0], true );
	double error = calc_bosted_errors( 3, x[0], false );

	return GEn + GEn*fabs(error);
}
Double_t bosted_GEn_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GEn = calc_bosted_parameterization( 3, x[0], true );
	double error = calc_bosted_errors( 3, x[0], false );

	return GEn - GEn*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GMn ------------
Double_t bosted_GMn_parameterization( Double_t *x, Double_t *par_G ){
	double GMn = calc_bosted_parameterization( 4, x[0], true );

	return GMn;
}
Double_t bosted_GMn_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GMn = calc_bosted_parameterization( 4, x[0], true );
	double error = sqrt(x[0])*GMn*calc_bosted_errors( 4, x[0], false );
	
	return GMn + fabs(error);
}
Double_t bosted_GMn_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GMn = calc_bosted_parameterization( 4, x[0], true );
	double error = sqrt(x[0])*GMn*calc_bosted_errors( 4, x[0], false );
	
	return GMn - fabs(error);
}
//-------------------------------------

TF1 *tf1_bosted_GEp, *tf1_bosted_GEp_error_high_fit, *tf1_bosted_GEp_error_low_fit;
TF1 *tf1_bosted_GMp, *tf1_bosted_GMp_error_high_fit, *tf1_bosted_GMp_error_low_fit;
TF1 *tf1_bosted_GEn, *tf1_bosted_GEn_error_high_fit, *tf1_bosted_GEn_error_low_fit;
TF1 *tf1_bosted_GMn, *tf1_bosted_GMn_error_high_fit, *tf1_bosted_GMn_error_low_fit;

//-------------------------------
//FF 0: All, 1: GEp, 2: GMp, 3: GEn, 4: GMn
void plot_bosted_FFs( int FF = 0, bool plot_jboyd = false, bool bool_plot_jboyd_GMn4 = false ){

	if( FF == 0 || FF == 1 ){
		double GEp_min_x = 0.1;
		double GEp_max_x = 20;
		double GEp_min_y = 0.5;
		double GEp_max_y = 1.5;
		double GEp_fit_Npx = 5000;

		TCanvas *c_bosted_GEp = new TCanvas("c_bosted_GEp", "c_bosted_GEp", 600, 500);
		TH1D *h_GEp = new TH1D("h_GEp", "GEp - Bosted Parameterization", GEp_fit_Npx, GEp_min_x, GEp_max_x);
		h_GEp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEp->GetXaxis()->SetTitleOffset(1.25f);

		h_GEp->GetYaxis()->SetTitle("G_{E}^{p}/G_{D}");
		h_GEp->GetYaxis()->SetTitleOffset(1.25f);
		h_GEp->SetLineColor(0);
		h_GEp->GetYaxis()->SetRangeUser(GEp_min_y, GEp_max_y);
		h_GEp->GetXaxis()->SetRangeUser(GEp_min_x, GEp_max_x);
		h_GEp->SetStats(0);
		h_GEp->Draw();

		c_bosted_GEp->SetLogx();
		c_bosted_GEp->Range(GEp_min_x, GEp_min_y, GEp_max_x, GEp_max_y);
		tf1_bosted_GEp = new TF1("tf1_bosted_GEp", bosted_GEp_parameterization, GEp_min_x, GEp_max_x, 0);
		tf1_bosted_GEp->SetNpx(GEp_fit_Npx);
		tf1_bosted_GEp->SetLineColor(1);
		tf1_bosted_GEp->SetLineWidth(2);
		tf1_bosted_GEp->Draw("same");

		tf1_bosted_GEp_error_high_fit = new TF1("tf1_bosted_GEp_error_high_fit", bosted_GEp_error_high_fit, GEp_min_x, GEp_max_x, 0);
		tf1_bosted_GEp_error_high_fit->SetNpx(GEp_fit_Npx);
		tf1_bosted_GEp_error_high_fit->SetLineColor(1);
		tf1_bosted_GEp_error_high_fit->SetLineWidth(1);
		tf1_bosted_GEp_error_high_fit->SetLineStyle(7);
		tf1_bosted_GEp_error_high_fit->Draw("same");

		tf1_bosted_GEp_error_low_fit = new TF1("tf1_bosted_GEp_error_low", bosted_GEp_error_low_fit, GEp_min_x, GEp_max_x, 0);
		tf1_bosted_GEp_error_low_fit->SetNpx(GEp_fit_Npx);
		tf1_bosted_GEp_error_low_fit->SetLineColor(1);
		tf1_bosted_GEp_error_low_fit->SetLineWidth(1);
		tf1_bosted_GEp_error_low_fit->SetLineStyle(7);
		tf1_bosted_GEp_error_low_fit->Draw("same");

		// fshade_btw_TF1s_w_3rd(c_bosted_GEp, h_GEp, tf1_bosted_GEp_error_high_fit, tf1_bosted_GEp_error_low_fit, tf1_bosted_GEp, GEp_min_y, GEp_max_y);
	
	}

	if( FF == 0 || FF == 2 ){
		double GMp_min_x = 0.1;
		double GMp_max_x = 30;
		double GMp_min_y = 0.6;
		double GMp_max_y = 1.1;
		double GMp_fit_Npx = 5000;

		TCanvas *c_bosted_GMp = new TCanvas("c_bosted_GMp", "c_bosted_GMp", 600, 500);
		TH1D *h_GMp = new TH1D("h_GMp", "GMp - Bosted Parameterization", GMp_fit_Npx, GMp_min_x, GMp_max_x);
		h_GMp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMp->GetXaxis()->SetTitleOffset(1.25f);

		h_GMp->GetYaxis()->SetTitle("G_{M}^{p}/G_{D}");
		h_GMp->GetYaxis()->SetTitleOffset(1.25f);
		h_GMp->SetLineColor(0);
		h_GMp->GetYaxis()->SetRangeUser(GMp_min_y, GMp_max_y);
		h_GMp->GetXaxis()->SetRangeUser(GMp_min_x, GMp_max_x);
		h_GMp->SetStats(0);
		h_GMp->Draw();

		c_bosted_GMp->SetLogx();
		c_bosted_GMp->Range(GMp_min_x, GMp_min_y, GMp_max_x, GMp_max_y);
		tf1_bosted_GMp = new TF1("tf1_bosted_GMp", bosted_GMp_parameterization, GMp_min_x, GMp_max_x, 0);
		tf1_bosted_GMp->SetNpx(GMp_fit_Npx);
		tf1_bosted_GMp->SetLineColor(1);
		tf1_bosted_GMp->SetLineWidth(2);
		tf1_bosted_GMp->Draw("same");

		tf1_bosted_GMp_error_high_fit = new TF1("tf1_bosted_GMp_error_high_fit", bosted_GMp_error_high_fit, GMp_min_x, GMp_max_x, 0);
		tf1_bosted_GMp_error_high_fit->SetNpx(GMp_fit_Npx);
		tf1_bosted_GMp_error_high_fit->SetLineColor(1);
		tf1_bosted_GMp_error_high_fit->SetLineWidth(1);
		tf1_bosted_GMp_error_high_fit->SetLineStyle(7);
		tf1_bosted_GMp_error_high_fit->Draw("same");

		tf1_bosted_GMp_error_low_fit = new TF1("tf1_bosted_GMp_error_low", bosted_GMp_error_low_fit, GMp_min_x, GMp_max_x, 0);
		tf1_bosted_GMp_error_low_fit->SetNpx(GMp_fit_Npx);
		tf1_bosted_GMp_error_low_fit->SetLineColor(1);
		tf1_bosted_GMp_error_low_fit->SetLineWidth(1);
		tf1_bosted_GMp_error_low_fit->SetLineStyle(7);
		tf1_bosted_GMp_error_low_fit->Draw("same");

		// fshade_btw_TF1s_w_3rd(c_bosted_GMp, h_GMp, tf1_bosted_GMp_error_high_fit, tf1_bosted_GMp_error_low_fit, tf1_bosted_GMp, GMp_min_y, GMp_max_y);
	
	}

	if( FF == 0 || FF == 3 ){
		double GEn_min_x = 0.2;
		double GEn_max_x = 5;
		double GEn_min_y = -0.2;
		double GEn_max_y = 0.6;
		double GEn_fit_Npx = 5000;

		TCanvas *c_bosted_GEn = new TCanvas("c_bosted_GEn", "c_bosted_GEn", 600, 500);
		TH1D *h_GEn = new TH1D("h_GEn", "GEn - Bosted Parameterization", GEn_fit_Npx, GEn_min_x, GEn_max_x);
		h_GEn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEn->GetXaxis()->SetTitleOffset(1.25f);

		h_GEn->GetYaxis()->SetTitle("G_{E}^{n}/G_{D}");
		h_GEn->GetYaxis()->SetTitleOffset(1.25f);
		h_GEn->SetLineColor(0);
		h_GEn->GetYaxis()->SetRangeUser(GEn_min_y, GEn_max_y);
		h_GEn->GetXaxis()->SetRangeUser(GEn_min_x, GEn_max_x);
		h_GEn->SetStats(0);
		h_GEn->Draw();

		c_bosted_GEn->SetLogx();
		c_bosted_GEn->Range(GEn_min_x, GEn_min_y, GEn_max_x, GEn_max_y);
		tf1_bosted_GEn = new TF1("tf1_bosted_GEn", bosted_GEn_parameterization, GEn_min_x, GEn_max_x, 0);
		tf1_bosted_GEn->SetNpx(GEn_fit_Npx);
		tf1_bosted_GEn->SetLineColor(1);
		tf1_bosted_GEn->SetLineWidth(2);
		tf1_bosted_GEn->Draw("same");

		tf1_bosted_GEn_error_high_fit = new TF1("tf1_bosted_GEn_error_high_fit", bosted_GEn_error_high_fit, GEn_min_x, GEn_max_x, 0);
		tf1_bosted_GEn_error_high_fit->SetNpx(GEn_fit_Npx);
		tf1_bosted_GEn_error_high_fit->SetLineColor(1);
		tf1_bosted_GEn_error_high_fit->SetLineWidth(1);
		tf1_bosted_GEn_error_high_fit->SetLineStyle(7);
		tf1_bosted_GEn_error_high_fit->Draw("same");

		tf1_bosted_GEn_error_low_fit = new TF1("tf1_bosted_GEn_error_low", bosted_GEn_error_low_fit, GEn_min_x, GEn_max_x, 0);
		tf1_bosted_GEn_error_low_fit->SetNpx(GEn_fit_Npx);
		tf1_bosted_GEn_error_low_fit->SetLineColor(1);
		tf1_bosted_GEn_error_low_fit->SetLineWidth(1);
		tf1_bosted_GEn_error_low_fit->SetLineStyle(7);
		tf1_bosted_GEn_error_low_fit->Draw("same");

		// fshade_btw_TF1s_w_3rd(c_bosted_GEn, h_GEn, tf1_bosted_GEn_error_high_fit, tf1_bosted_GEn_error_low_fit, tf1_bosted_GEn, GEn_min_y, GEn_max_y);
	
	}

	if( FF == 0 || FF == 4 ){
		double GMn_min_x = 0.1;
		double GMn_max_x = 10;
		double GMn_min_y = 0.6;
		double GMn_max_y = 1.2;
		double GMn_fit_Npx = 5000;

		TCanvas *c_bosted_GMn = new TCanvas("c_bosted_GMn", "c_bosted_GMn", 600, 500);
		TH1D *h_GMn = new TH1D("h_GMn", "GMn - Bosted Parameterization", GMn_fit_Npx, GMn_min_x, GMn_max_x);
		h_GMn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMn->GetXaxis()->SetTitleOffset(1.25f);

		h_GMn->GetYaxis()->SetTitle("G_{M}^{n}/G_{D}");
		h_GMn->GetYaxis()->SetTitleOffset(1.25f);
		h_GMn->SetLineColor(0);
		h_GMn->GetYaxis()->SetRangeUser(GMn_min_y, GMn_max_y);
		h_GMn->GetXaxis()->SetRangeUser(GMn_min_x, GMn_max_x);
		h_GMn->SetStats(0);
		h_GMn->Draw();

		c_bosted_GMn->SetLogx();
		c_bosted_GMn->Range(GMn_min_x, GMn_min_y, GMn_max_x, GMn_max_y);
		tf1_bosted_GMn = new TF1("tf1_bosted_GMn", bosted_GMn_parameterization, GMn_min_x, GMn_max_x, 0);
		tf1_bosted_GMn->SetNpx(GMn_fit_Npx);
		tf1_bosted_GMn->SetLineColor(1);
		tf1_bosted_GMn->SetLineWidth(2);
		tf1_bosted_GMn->Draw("same");

		tf1_bosted_GMn_error_high_fit = new TF1("tf1_bosted_GMn_error_high_fit", bosted_GMn_error_high_fit, GMn_min_x, GMn_max_x, 0);
		tf1_bosted_GMn_error_high_fit->SetNpx(GMn_fit_Npx);
		tf1_bosted_GMn_error_high_fit->SetLineColor(1);
		tf1_bosted_GMn_error_high_fit->SetLineWidth(1);
		tf1_bosted_GMn_error_high_fit->SetLineStyle(7);
		tf1_bosted_GMn_error_high_fit->Draw("same");

		tf1_bosted_GMn_error_low_fit = new TF1("tf1_bosted_GMn_error_low", bosted_GMn_error_low_fit, GMn_min_x, GMn_max_x, 0);
		tf1_bosted_GMn_error_low_fit->SetNpx(GMn_fit_Npx);
		tf1_bosted_GMn_error_low_fit->SetLineColor(1);
		tf1_bosted_GMn_error_low_fit->SetLineWidth(1);
		tf1_bosted_GMn_error_low_fit->SetLineStyle(7);
		tf1_bosted_GMn_error_low_fit->Draw("same");

		// fshade_btw_TF1s_w_3rd(c_bosted_GMn, h_GMn, tf1_bosted_GMn_error_high_fit, tf1_bosted_GMn_error_low_fit, tf1_bosted_GMn, GMn_min_y, GMn_max_y);
	
	}

	if( plot_jboyd ){

		TString drawStyle = "*same";
		if( bool_plot_jboyd_GMn4 ){
			plot_jboyd_data( bool_plot_jboyd_GMn4 );
		}
		else{
			plot_jboyd_data( false );
		}
		
		tge_GMn_SBS8_jboyd_stat->Draw("same");
		tge_GMn_SBS9_jboyd_stat->Draw(drawStyle.Data());
		tge_GMn_SBS8_jboyd->Draw(drawStyle.Data());
		tge_GMn_SBS9_jboyd->Draw(drawStyle.Data());

		if( bool_plot_jboyd_GMn4 ){
			tge_GMn_SBS4_jboyd->Draw(drawStyle.Data());
		}

		TLegend *tleg_GMn_jboyd;
		tleg_GMn_jboyd = new TLegend(0.125, 0.20, 0.80, 0.35);
		tleg_GMn_jboyd->SetTextSize(0.03);
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd, Form("SBS8 G_{M}^{n} = %.4f +/- %.4f +/- %.4f", jboyd_GMn_vec[0], jboyd_GMn_syst_error_vec[0], jboyd_GMn_stat_error_vec[0] - jboyd_GMn_syst_error_vec[0]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd, Form("SBS9 G_{M}^{n} = %.4f +/- %.4f +/- %.4f", jboyd_GMn_vec[1], jboyd_GMn_syst_error_vec[1], jboyd_GMn_stat_error_vec[1] - jboyd_GMn_syst_error_vec[1]));

		if( bool_plot_jboyd_GMn4 ){
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 G_{M}^{n} = %.4f +/- %.4f +/- %.4f", jboyd_GMn_vec[2], jboyd_GMn_syst_error_vec[2], jboyd_GMn_stat_error_vec[2] - jboyd_GMn_syst_error_vec[2]));
		}
		tleg_GMn_jboyd->Draw(drawStyle.Data());			


	}

}

void general_plot_bosted(TCanvas *c_to_plot_on, TLegend *tleg_to_add = NULL, int setLineColor = 1, int setLineStyle = 9, int setLineWidth = 2){

		double GMn_min_x = 0.1;
		double GMn_max_x = 10;
		double GMn_min_y = 0.6;
		double GMn_max_y = 1.2;
		double GMn_fit_Npx = 5000;

		c_to_plot_on->cd();
		TH1D *h_GMn = new TH1D("h_GMn", "GMn - Bosted Parameterization", GMn_fit_Npx, GMn_min_x, GMn_max_x);
		// h_GMn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		// h_GMn->GetXaxis()->SetTitleOffset(1.25f);

		// h_GMn->GetYaxis()->SetTitle("G_{M}^{n}/G_{D}");
		// h_GMn->GetYaxis()->SetTitleOffset(1.25f);
		// h_GMn->SetLineColor(0);
		// h_GMn->GetYaxis()->SetRangeUser(GMn_min_y, GMn_max_y);
		// h_GMn->GetXaxis()->SetRangeUser(GMn_min_x, GMn_max_x);
		// h_GMn->SetStats(0);
		h_GMn->Draw("same");

		// c_to_plot_on->SetLogx();
		// c_to_plot_on->Range(GMn_min_x, GMn_min_y, GMn_max_x, GMn_max_y);
		tf1_bosted_GMn = new TF1("tf1_bosted_GMn", bosted_GMn_parameterization, GMn_min_x, GMn_max_x, 0);
		tf1_bosted_GMn->SetNpx(GMn_fit_Npx);
		tf1_bosted_GMn->SetLineColor(setLineColor);
		tf1_bosted_GMn->SetLineStyle(setLineStyle);
		tf1_bosted_GMn->SetLineWidth(setLineWidth);
		tf1_bosted_GMn->Draw("same");

		// tf1_bosted_GMn_error_high_fit = new TF1("tf1_bosted_GMn_error_high_fit", bosted_GMn_error_high_fit, GMn_min_x, GMn_max_x, 0);
		// tf1_bosted_GMn_error_high_fit->SetNpx(GMn_fit_Npx);
		// tf1_bosted_GMn_error_high_fit->SetLineColor(1);
		// tf1_bosted_GMn_error_high_fit->SetLineWidth(1);
		// tf1_bosted_GMn_error_high_fit->SetLineStyle(7);
		// tf1_bosted_GMn_error_high_fit->Draw("same");

		// tf1_bosted_GMn_error_low_fit = new TF1("tf1_bosted_GMn_error_low", bosted_GMn_error_low_fit, GMn_min_x, GMn_max_x, 0);
		// tf1_bosted_GMn_error_low_fit->SetNpx(GMn_fit_Npx);
		// tf1_bosted_GMn_error_low_fit->SetLineColor(1);
		// tf1_bosted_GMn_error_low_fit->SetLineWidth(1);
		// tf1_bosted_GMn_error_low_fit->SetLineStyle(7);
		// tf1_bosted_GMn_error_low_fit->Draw("same");

		// fshade_btw_TF1s_w_3rd(c_bosted_GMn, h_GMn, tf1_bosted_GMn_error_high_fit, tf1_bosted_GMn_error_low_fit, tf1_bosted_GMn, GMn_min_y, GMn_max_y);

		if( tleg_to_add != NULL ){
			tleg_to_add->AddEntry(tf1_bosted_GMn, "Bosted");
			c_to_plot_on->Update();
			c_to_plot_on->Modified();
		}


}

#endif