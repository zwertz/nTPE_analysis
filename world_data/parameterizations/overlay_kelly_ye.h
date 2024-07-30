#ifndef OVERLAY_KELLY_YE_H
#define OVERLAY_KELLY_YE_H

#include "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/world_data/parameterizations/GetFFfunctions.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/plot_pimping.h"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/world_data/jboyd_data_points.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/calc_errors.h"


double par_G[4], par_G_error[4];

double par_GEp[5] = {1.0, -0.24, 10.98, 12.82, 21.97};
double par_GEp_error[5] = {0.0, 0.12, 0.19, 1.1, 6.8};

double par_GMp[5] = {1.0, 0.12, 10.97, 18.86, 6.55};
double par_GMp_error[5] = {0.0, 0.04, 0.11, 0.28, 1.2};

double par_GMn[5] = {1.0, 2.33, 14.72, 24.20, 84.1};
double par_GMn_error[5] = {0.0, 1.4, 1.7, 9.8, 41};

double par_GEn[2] = {1.70, 3.30};
double par_GEn_error[2] = {0.04, 0.32};

TCanvas *c_overlay_kelly_ye_GEp, *c_overlay_kelly_ye_GMp, *c_overlay_kelly_ye_GEn, *c_overlay_kelly_ye_GMn;

TF1 *tf1_ye_arrington_GEp_fit, *tf1_ye_arrington_GEp_fit_error_high, *tf1_ye_arrington_GEp_fit_error_low;
TF1 *tf1_ye_arrington_GMp_fit, *tf1_ye_arrington_GMp_fit_error_high, *tf1_ye_arrington_GMp_fit_error_low;
TF1 *tf1_ye_arrington_GEn_fit, *tf1_ye_arrington_GEn_fit_error_high, *tf1_ye_arrington_GEn_fit_error_low;
TF1 *tf1_ye_arrington_GMn_fit, *tf1_ye_arrington_GMn_fit_error_high, *tf1_ye_arrington_GMn_fit_error_low;

TGraph *tg_ye_arrington;
TGraphErrors *tg_GMn_exp_all;
TPolyLine *GMn_shaded_area;
auto *mg_ye_arrington = new TMultiGraph();

TGraph *graph = new TGraph();

vector<double> ye_GMn_vec;
vector<double> ye_Q2_vec;
vector<double> GMn_calc_vec, G_D_calc_vec, GMn_GD_calc_vec;
vector<double> q2_vec, z_vec, z_numer_vec, z_denom_vec;

double par_aN_vec[11] = {0.257758326959, -1.079540642058, 1.182183812195, 0.711015085833, -1.348080936796, -1.662444025208,
	2.624354426029, 1.751234494568, -4.922300878888, 3.197892727312, -0.712072389946};

double par_aN_err_vec[11];

double par_low_Q2[2] = {-0.20765505, 0.99767103};
double par_mid_Q2[15] = {-2.06920873e0, 6.43156400e-2, -3.55593786e-1, 4.14897660e-1, 1.95746824e0, 2.70525700e-1, -1.52685784e0,
	-4.43527359e-1, 5.16884065e-1, 2.07915837e-1, -7.48665703e-2, -4.25411431e-2, 1.54965016e-3, 3.25322279e-3, 4.20819518e-4
  };
double par_high_Q2[3] = {0.50859057, 1.96863291, 0.23213950};


//---------------GEp------------------
Double_t GEp_parameterization(Double_t *x, Double_t *par_aN_vec ){
	double GEp_calc = GetFFvalue(1, x[0]);
	return GEp_calc;
}

Double_t GEp_error_high_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GEp_calc = GetFFvalue(1, x[0]);
	double error = GetFFerror(1, x[0]);

	return GEp_calc + error;
}

Double_t GEp_error_low_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GEp_calc = GetFFvalue(1, x[0]);
	double error = GetFFerror(1, x[0]);

	return GEp_calc - error;
}
//-------------------------------
//---------------GMp------------------
Double_t GMp_parameterization(Double_t *x, Double_t *par_aN_vec ){
	double GMp_calc = GetFFvalue(2, x[0]);
	return GMp_calc;
}

Double_t GMp_error_high_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GMp_calc = GetFFvalue(2, x[0]);
	double error = GetFFerror(2, x[0]);

	return GMp_calc + error;
}

Double_t GMp_error_low_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GMp_calc = GetFFvalue(2, x[0]);
	double error = GetFFerror(2, x[0]);

	return GMp_calc - error;
}
//-------------------------------
//---------------GEn------------------
Double_t GEn_parameterization(Double_t *x, Double_t *par_aN_vec ){
	double GEn_calc = GetFFvalue(3, x[0]);
	return GEn_calc;
}

Double_t GEn_error_high_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GEn_calc = GetFFvalue(3, x[0]);
	double error = GetFFerror(3, x[0]);

	return GEn_calc + error;
}

Double_t GEn_error_low_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GEn_calc = GetFFvalue(3, x[0]);
	double error = GetFFerror(3, x[0]);

	return GEn_calc - error;
}
//-------------------------------
//---------- GMn -------------------------
Double_t GMn_parameterization(Double_t *x, Double_t *par_aN_vec ){
	double GMn_calc = GetFFvalue(4, x[0]);
	return GMn_calc;
}

Double_t GMn_error_high_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GMn_calc = GetFFvalue(4, x[0]);
	double error = GetFFerror(4, x[0]);

	return GMn_calc + error;
}

Double_t GMn_error_low_fit(Double_t *x, Double_t *par_aN_vec_err ){
	double GMn_calc = GetFFvalue(4, x[0]);
	double error = GetFFerror(4, x[0]);

	return GMn_calc - error;
}

//-------------------------------
//FF 0: All, 1: GEp, 2: GMp, 3: GEn, 4: GMn
void plot_ye_arrington( int FF = 0, bool plot_jboyd = false, bool bool_plot_jboyd_GMn4 = false ){

	// int ye_line_color = 632;
	// int ye_shade_color = 622;

	int ye_line_color = 600;
	int ye_shade_color = 590;

	tf1_ye_arrington_GMp_fit = new TF1("tf1_ye_arrington_GMp_fit", GMp_parameterization, 0, 20, 11);
	tf1_ye_arrington_GEn_fit = new TF1("tf1_ye_arrington_GEn_fit", GEn_parameterization, 0, 20, 11);

	if( FF == 0 || FF == 1 ){
		double GEp_min_x = 0.01;
		double GEp_max_x = 100;
		double GEp_min_y = -1.0;
		double GEp_max_y = 1.1;
		double GEp_fit_Npx = 5000;

		c_overlay_kelly_ye_GEp->cd();
		TH1D *h_GEp = new TH1D("h_GEp", "GEp - Ye Parameterization", GEp_fit_Npx, GEp_min_x, GEp_max_x);
		h_GEp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEp->GetXaxis()->SetTitleOffset(1.25f);

		h_GEp->GetYaxis()->SetTitle("G_{E}^{p}/G_{D}");
		h_GEp->GetYaxis()->SetTitleOffset(1.25f);
		h_GEp->SetLineColor(0);
		h_GEp->GetYaxis()->SetRangeUser(GEp_min_y, GEp_max_y);
		h_GEp->GetXaxis()->SetRangeUser(GEp_min_x, GEp_max_x);
		h_GEp->SetStats(0);
		h_GEp->Draw("same");

		c_overlay_kelly_ye_GEp->SetLogx();
		c_overlay_kelly_ye_GEp->Range(GEp_min_x, GEp_min_y, GEp_max_x, GEp_max_y);
		tf1_ye_arrington_GEp_fit = new TF1("tf1_ye_arrington_GEp_fit", GEp_parameterization, GEp_min_x, GEp_max_x, 0);
		tf1_ye_arrington_GEp_fit->SetNpx(GEp_fit_Npx);
		tf1_ye_arrington_GEp_fit->SetLineColor(ye_line_color);
		tf1_ye_arrington_GEp_fit->SetLineWidth(2);
		tf1_ye_arrington_GEp_fit->Draw("same");

		tf1_ye_arrington_GEp_fit_error_high = new TF1("tf1_ye_arrington_GEp_fit_error_high", GEp_error_high_fit, GEp_min_x, GEp_max_x, 0);
		tf1_ye_arrington_GEp_fit_error_high->SetNpx(GEp_fit_Npx);
		tf1_ye_arrington_GEp_fit_error_high->SetLineColor(ye_line_color);
		tf1_ye_arrington_GEp_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GEp_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GEp_fit_error_high->Draw("same");

		tf1_ye_arrington_GEp_fit_error_low = new TF1("tf1_ye_arrington_GEp_fit_error_low", GEp_error_low_fit, GEp_min_x, GEp_max_x, 0);
		tf1_ye_arrington_GEp_fit_error_low->SetNpx(GEp_fit_Npx);
		tf1_ye_arrington_GEp_fit_error_low->SetLineColor(ye_line_color);
		tf1_ye_arrington_GEp_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GEp_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GEp_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GEp, h_GEp, tf1_ye_arrington_GEp_fit_error_high, tf1_ye_arrington_GEp_fit_error_low, tf1_ye_arrington_GEp_fit, GEp_min_y, GEp_max_y, ye_line_color, ye_shade_color);

	}

	if( FF == 0 || FF == 2 ){
		double GMp_min_x = 0.01;
		double GMp_max_x = 100;
		double GMp_min_y = 0.4;
		double GMp_max_y = 1.15;
		double GMp_fit_Npx = 5000;

		c_overlay_kelly_ye_GMp->cd();
		TH1D *h_GMp = new TH1D("h_GMp", "GMp - Ye Parameterization", GMp_fit_Npx, GMp_min_x, GMp_max_x);
		h_GMp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMp->GetXaxis()->SetTitleOffset(1.25f);

		h_GMp->GetYaxis()->SetTitle("G_{M}^{p}/(#mu_{p} G_{D})");
		h_GMp->GetYaxis()->SetTitleOffset(1.25f);
		h_GMp->SetLineColor(0);
		h_GMp->GetYaxis()->SetRangeUser(GMp_min_y, GMp_max_y);
		h_GMp->GetXaxis()->SetRangeUser(GMp_min_x, GMp_max_x);
		h_GMp->SetStats(0);
		h_GMp->Draw("same");

		c_overlay_kelly_ye_GMp->SetLogx();
		c_overlay_kelly_ye_GMp->Range(GMp_min_x, GMp_min_y, GMp_max_x, GMp_max_y);
		tf1_ye_arrington_GMp_fit = new TF1("tf1_ye_arrington_GMp_fit", GMp_parameterization, GMp_min_x, GMp_max_x, 0);
		tf1_ye_arrington_GMp_fit->SetNpx(GMp_fit_Npx);
		tf1_ye_arrington_GMp_fit->SetLineColor(ye_line_color);
		tf1_ye_arrington_GMp_fit->SetLineWidth(2);
		tf1_ye_arrington_GMp_fit->Draw("same");

		tf1_ye_arrington_GMp_fit_error_high = new TF1("tf1_ye_arrington_GMp_fit_error_high", GMp_error_high_fit, GMp_min_x, GMp_max_x, 0);
		tf1_ye_arrington_GMp_fit_error_high->SetNpx(GMp_fit_Npx);
		tf1_ye_arrington_GMp_fit_error_high->SetLineColor(ye_line_color);
		tf1_ye_arrington_GMp_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GMp_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GMp_fit_error_high->Draw("same");

		tf1_ye_arrington_GMp_fit_error_low = new TF1("tf1_ye_arrington_GMp_fit_error_low", GMp_error_low_fit, GMp_min_x, GMp_max_x, 0);
		tf1_ye_arrington_GMp_fit_error_low->SetNpx(GMp_fit_Npx);
		tf1_ye_arrington_GMp_fit_error_low->SetLineColor(ye_line_color);
		tf1_ye_arrington_GMp_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GMp_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GMp_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GMp, h_GMp, tf1_ye_arrington_GMp_fit_error_high, tf1_ye_arrington_GMp_fit_error_low, tf1_ye_arrington_GMp_fit, GMp_min_y, GMp_max_y, ye_line_color, ye_shade_color);

	}

	if( FF == 0 || FF == 3 ){
		double GEn_min_x = 0.01;
		double GEn_max_x = 20;
		double GEn_min_y = -0.025;
		double GEn_max_y = 0.8;
		double GEn_fit_Npx = 5000;

		c_overlay_kelly_ye_GEn->cd();
		TH1D *h_GEn = new TH1D("h_GEn", "GEn - Ye Parameterization", GEn_fit_Npx, GEn_min_x, GEn_max_x);
		h_GEn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEn->GetXaxis()->SetTitleOffset(1.25f);

		h_GEn->GetYaxis()->SetTitle("G_{E}^{n}/G_{D}");
		h_GEn->GetYaxis()->SetTitleOffset(1.25f);
		h_GEn->SetLineColor(0);
		h_GEn->GetYaxis()->SetRangeUser(GEn_min_y, GEn_max_y);
		h_GEn->GetXaxis()->SetRangeUser(GEn_min_x, GEn_max_x);
		h_GEn->SetStats(0);
		h_GEn->Draw("same");

		c_overlay_kelly_ye_GEn->SetLogx();
		c_overlay_kelly_ye_GEn->Range(GEn_min_x, GEn_min_y, GEn_max_x, GEn_max_y);
		tf1_ye_arrington_GEn_fit = new TF1("tf1_ye_arrington_GEn_fit", GEn_parameterization, GEn_min_x, GEn_max_x, 0);
		tf1_ye_arrington_GEn_fit->SetNpx(GEn_fit_Npx);
		tf1_ye_arrington_GEn_fit->SetLineColor(ye_line_color);
		tf1_ye_arrington_GEn_fit->SetLineWidth(2);
		tf1_ye_arrington_GEn_fit->Draw("same");

		tf1_ye_arrington_GEn_fit_error_high = new TF1("tf1_ye_arrington_GEn_fit_error_high", GEn_error_high_fit, GEn_min_x, GEn_max_x, 0);
		tf1_ye_arrington_GEn_fit_error_high->SetNpx(GEn_fit_Npx);
		tf1_ye_arrington_GEn_fit_error_high->SetLineColor(ye_line_color);
		tf1_ye_arrington_GEn_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GEn_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GEn_fit_error_high->Draw("same");

		tf1_ye_arrington_GEn_fit_error_low = new TF1("tf1_ye_arrington_GEn_fit_error_low", GEn_error_low_fit, GEn_min_x, GEn_max_x, 0);
		tf1_ye_arrington_GEn_fit_error_low->SetNpx(GEn_fit_Npx);
		tf1_ye_arrington_GEn_fit_error_low->SetLineColor(ye_line_color);
		tf1_ye_arrington_GEn_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GEn_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GEn_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GEn, h_GEn, tf1_ye_arrington_GEn_fit_error_high, tf1_ye_arrington_GEn_fit_error_low, tf1_ye_arrington_GEn_fit, GEn_min_y, GEn_max_y, ye_line_color, ye_shade_color);

	}

	if( FF == 0 || FF == 4 ){
		double GMn_min_x = 0.01;
		double GMn_max_x = 20;
		double GMn_min_y = 0.175;
		double GMn_max_y = 1.2;
		double GMn_fit_Npx = 5000;

		c_overlay_kelly_ye_GMn->SetGrid();	
		c_overlay_kelly_ye_GMn->cd();
		TH1D *h_GMn = new TH1D("h_GMn", "GMn - Kelly and Ye Parameterizations", GMn_fit_Npx, GMn_min_x, GMn_max_x);
		h_GMn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMn->GetXaxis()->SetTitleOffset(1.25f);

		h_GMn->GetYaxis()->SetTitle("G_{M}^{n}/(#mu_{n} G_{D})");
		h_GMn->GetYaxis()->SetTitleOffset(1.25f);
		h_GMn->SetLineColor(0);
		h_GMn->GetYaxis()->SetRangeUser(GMn_min_y, GMn_max_y);
		h_GMn->GetXaxis()->SetRangeUser(GMn_min_x, GMn_max_x);
		h_GMn->SetStats(0);
		h_GMn->Draw("same");

		c_overlay_kelly_ye_GMn->SetLogx();
		c_overlay_kelly_ye_GMn->Range(GMn_min_x, GMn_min_y, GMn_max_x, GMn_max_y);
		tf1_ye_arrington_GMn_fit = new TF1("tf1_ye_arrington_GMn_fit", GMn_parameterization, GMn_min_x, GMn_max_x, 0);
		tf1_ye_arrington_GMn_fit->SetNpx(GMn_fit_Npx);
		tf1_ye_arrington_GMn_fit->SetLineColor(ye_line_color);
		tf1_ye_arrington_GMn_fit->SetLineWidth(2);
		tf1_ye_arrington_GMn_fit->Draw("same");

		tf1_ye_arrington_GMn_fit_error_high = new TF1("tf1_ye_arrington_GMn_fit_error_high", GMn_error_high_fit, GMn_min_x, GMn_max_x, 0);
		tf1_ye_arrington_GMn_fit_error_high->SetNpx(GMn_fit_Npx);
		tf1_ye_arrington_GMn_fit_error_high->SetLineColor(ye_line_color);
		tf1_ye_arrington_GMn_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GMn_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GMn_fit_error_high->Draw("same");

		tf1_ye_arrington_GMn_fit_error_low = new TF1("tf1_ye_arrington_GMn_fit_error_low", GMn_error_low_fit, GMn_min_x, GMn_max_x, 0);
		tf1_ye_arrington_GMn_fit_error_low->SetNpx(GMn_fit_Npx);
		tf1_ye_arrington_GMn_fit_error_low->SetLineColor(ye_line_color);
		tf1_ye_arrington_GMn_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GMn_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GMn_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GMn, h_GMn, tf1_ye_arrington_GMn_fit_error_high, tf1_ye_arrington_GMn_fit_error_low, tf1_ye_arrington_GMn_fit, GMn_min_y, GMn_max_y, ye_line_color, ye_shade_color);

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
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd, Form("SBS8 G_{M}^{n} = %.4f +/- %.4f +/- %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[0], jboyd_GMn_syst_error_vec[0], jboyd_GMn_stat_error_vec[0] - jboyd_GMn_syst_error_vec[0], jboyd_Q2_vec[0]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd, Form("SBS9 G_{M}^{n} = %.4f +/- %.4f +/- %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[1], jboyd_GMn_syst_error_vec[1], jboyd_GMn_stat_error_vec[1] - jboyd_GMn_syst_error_vec[1], jboyd_Q2_vec[1]));

		if( bool_plot_jboyd_GMn4 ){
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 G_{M}^{n} = %.4f +/- %.4f +/- %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[2], jboyd_GMn_syst_error_vec[2], jboyd_GMn_stat_error_vec[2] - jboyd_GMn_syst_error_vec[2], jboyd_Q2_vec[2]));
		}

		tleg_GMn_jboyd->Draw(drawStyle.Data());			

		c_overlay_kelly_ye_GMn->cd();
		tge_GMn_SBS4_jboyd_stat->Draw("*same");
		tge_GMn_SBS4_jboyd->Draw("*same");
		tge_GMn_SBS4_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS4_jboyd->Draw("same");
		tg_GMn_SBS4_jboyd->Draw("*same");
		tg_GMn_SBS4_jboyd->SetLineColor(ye_line_color);
		tg_GMn_SBS4_jboyd->SetMarkerColor(1);
		tg_GMn_SBS4_jboyd->SetMarkerStyle(kFullCircle);
		tg_GMn_SBS4_jboyd->Draw("same");

		tge_GMn_SBS8_jboyd_stat->Draw("*same");
		tge_GMn_SBS8_jboyd->Draw("*same");
		tge_GMn_SBS8_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS8_jboyd->Draw("same");
		tg_GMn_SBS8_jboyd->Draw("*same");
		tg_GMn_SBS8_jboyd->SetLineColor(ye_line_color);
		tg_GMn_SBS8_jboyd->SetMarkerColor(1);
		tg_GMn_SBS8_jboyd->SetMarkerStyle(kFullCircle);
		tg_GMn_SBS8_jboyd->Draw("same");

		tge_GMn_SBS9_jboyd_stat->Draw("*same");
		tge_GMn_SBS9_jboyd->Draw("*same");
		tge_GMn_SBS9_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS9_jboyd->Draw("same");
		tg_GMn_SBS9_jboyd->Draw("*same");
		tg_GMn_SBS9_jboyd->SetLineColor(ye_line_color);
		tg_GMn_SBS9_jboyd->SetMarkerColor(1);
		tg_GMn_SBS9_jboyd->SetMarkerStyle(kFullCircle);
		tg_GMn_SBS9_jboyd->Draw("same");

	}

}
Double_t calc_overlay_kelly_ye_parameterization( int FF, double q2, bool normed = true ){
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

Double_t calc_overlay_kelly_ye_errors( int FF, double q2, bool normed = true ){
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
	double GEp = calc_overlay_kelly_ye_parameterization( 1, x[0], true );

	return GEp;
}
Double_t kelly_GEp_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GEp = calc_overlay_kelly_ye_parameterization( 1, x[0], true );
	double error = calc_overlay_kelly_ye_errors( 1, x[0], true );

	return GEp + GEp*fabs(error);
}
Double_t kelly_GEp_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GEp = calc_overlay_kelly_ye_parameterization( 1, x[0], true );
	double error = calc_overlay_kelly_ye_errors( 1, x[0], true );

	return GEp - GEp*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GMp ------------
Double_t kelly_GMp_parameterization( Double_t *x, Double_t *par_G ){
	double GMp = calc_overlay_kelly_ye_parameterization( 2, x[0], true );

	return GMp;
}
Double_t kelly_GMp_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GMp = calc_overlay_kelly_ye_parameterization( 2, x[0], true );
	double error = calc_overlay_kelly_ye_errors( 2, x[0], false );

	return GMp + x[0]*GMp*fabs(error);
}
Double_t kelly_GMp_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GMp = calc_overlay_kelly_ye_parameterization( 2, x[0], true );
	double error = calc_overlay_kelly_ye_errors( 2, x[0], false );

	return GMp - x[0]*GMp*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GEn ------------
Double_t kelly_GEn_parameterization( Double_t *x, Double_t *par_G ){
	double GEn = calc_overlay_kelly_ye_parameterization( 3, x[0], true );

	return GEn;
}
Double_t kelly_GEn_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GEn = calc_overlay_kelly_ye_parameterization( 3, x[0], true );
	double error = calc_overlay_kelly_ye_errors( 3, x[0], false );

	return GEn + GEn*fabs(error);
}
Double_t kelly_GEn_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GEn = calc_overlay_kelly_ye_parameterization( 3, x[0], true );
	double error = calc_overlay_kelly_ye_errors( 3, x[0], false );

	return GEn - GEn*fabs(error);
}
//-------------------------------------

//-------------------------------------
//----------------- GMn ------------
Double_t kelly_GMn_parameterization( Double_t *x, Double_t *par_G ){
	double GMn = calc_overlay_kelly_ye_parameterization( 4, x[0], true );

	return GMn;
}
Double_t kelly_GMn_error_high_fit( Double_t *x, Double_t *par_G_error ){
	double GMn = calc_overlay_kelly_ye_parameterization( 4, x[0], true );
	double error = sqrt(x[0])*GMn*calc_overlay_kelly_ye_errors( 4, x[0], false );
	
	return GMn + fabs(error);
}
Double_t kelly_GMn_error_low_fit( Double_t *x, Double_t *par_G_error ){
	double GMn = calc_overlay_kelly_ye_parameterization( 4, x[0], true );
	double error = sqrt(x[0])*GMn*calc_overlay_kelly_ye_errors( 4, x[0], false );
	
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

	// int kelly_line_color = 600;
	// int kelly_shade_color = 590;

	int kelly_line_color = 632;
	int kelly_shade_color = 622;

	if( FF == 0 || FF == 1 ){
		double GEp_min_x = 0.003;
		double GEp_max_x = 10;
		double GEp_min_y = -0.1;
		double GEp_max_y = 1.2;
		double GEp_fit_Npx = 5000;

		c_overlay_kelly_ye_GEp->cd();
		TH1D *h_GEp = new TH1D("h_GEp", "GEp - Kelly Parameterization", GEp_fit_Npx, GEp_min_x, GEp_max_x);
		h_GEp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEp->GetXaxis()->SetTitleOffset(1.25f);

		h_GEp->GetYaxis()->SetTitle("G_{E}^{p}/G_{D}");
		h_GEp->GetYaxis()->SetTitleOffset(1.25f);
		h_GEp->SetLineColor(0);
		h_GEp->GetYaxis()->SetRangeUser(GEp_min_y, GEp_max_y);
		h_GEp->GetXaxis()->SetRangeUser(GEp_min_x, GEp_max_x);
		h_GEp->SetStats(0);
		h_GEp->Draw("same");

		c_overlay_kelly_ye_GEp->SetLogx();
		c_overlay_kelly_ye_GEp->Range(GEp_min_x, GEp_min_y, GEp_max_x, GEp_max_y);
		tf1_kelly_GEp = new TF1("tf1_kelly_GEp", kelly_GEp_parameterization, GEp_min_x, GEp_max_x, 0);
		tf1_kelly_GEp->SetNpx(GEp_fit_Npx);
		tf1_kelly_GEp->SetLineColor(kelly_line_color);
		tf1_kelly_GEp->SetLineWidth(2);
		tf1_kelly_GEp->Draw("same");

		tf1_kelly_GEp_error_high_fit = new TF1("tf1_kelly_GEp_error_high_fit", kelly_GEp_error_high_fit, GEp_min_x, GEp_max_x, 0);
		tf1_kelly_GEp_error_high_fit->SetNpx(GEp_fit_Npx);
		tf1_kelly_GEp_error_high_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GEp_error_high_fit->SetLineWidth(1);
		tf1_kelly_GEp_error_high_fit->SetLineStyle(7);
		tf1_kelly_GEp_error_high_fit->Draw("same");

		tf1_kelly_GEp_error_low_fit = new TF1("tf1_kelly_GEp_error_low", kelly_GEp_error_low_fit, GEp_min_x, GEp_max_x, 0);
		tf1_kelly_GEp_error_low_fit->SetNpx(GEp_fit_Npx);
		tf1_kelly_GEp_error_low_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GEp_error_low_fit->SetLineWidth(1);
		tf1_kelly_GEp_error_low_fit->SetLineStyle(7);
		tf1_kelly_GEp_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GEp, h_GEp, tf1_kelly_GEp_error_high_fit, tf1_kelly_GEp_error_low_fit, tf1_kelly_GEp, GEp_min_y, GEp_max_y, kelly_line_color, kelly_shade_color);
	
	}

	if( FF == 0 || FF == 2 ){
		double GMp_min_x = 0.01;
		double GMp_max_x = 60;
		double GMp_min_y = 0.6;
		double GMp_max_y = 1.15;
		double GMp_fit_Npx = 5000;

		c_overlay_kelly_ye_GMp->cd();
		TH1D *h_GMp = new TH1D("h_GMp", "GMp - Kelly Parameterization", GMp_fit_Npx, GMp_min_x, GMp_max_x);
		h_GMp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMp->GetXaxis()->SetTitleOffset(1.25f);

		h_GMp->GetYaxis()->SetTitle("G_{M}^{p}/(#mu_{p} G_{D})");
		h_GMp->GetYaxis()->SetTitleOffset(1.25f);
		h_GMp->SetLineColor(0);
		h_GMp->GetYaxis()->SetRangeUser(GMp_min_y, GMp_max_y);
		h_GMp->GetXaxis()->SetRangeUser(GMp_min_x, GMp_max_x);
		h_GMp->SetStats(0);
		h_GMp->Draw("same");

		c_overlay_kelly_ye_GMp->SetLogx();
		c_overlay_kelly_ye_GMp->Range(GMp_min_x, GMp_min_y, GMp_max_x, GMp_max_y);
		tf1_kelly_GMp = new TF1("tf1_kelly_GMp", kelly_GMp_parameterization, GMp_min_x, GMp_max_x, 0);
		tf1_kelly_GMp->SetNpx(GMp_fit_Npx);
		tf1_kelly_GMp->SetLineColor(kelly_line_color);
		tf1_kelly_GMp->SetLineWidth(2);
		tf1_kelly_GMp->Draw("same");

		tf1_kelly_GMp_error_high_fit = new TF1("tf1_kelly_GMp_error_high_fit", kelly_GMp_error_high_fit, GMp_min_x, GMp_max_x, 0);
		tf1_kelly_GMp_error_high_fit->SetNpx(GMp_fit_Npx);
		tf1_kelly_GMp_error_high_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GMp_error_high_fit->SetLineWidth(1);
		tf1_kelly_GMp_error_high_fit->SetLineStyle(7);
		tf1_kelly_GMp_error_high_fit->Draw("same");

		tf1_kelly_GMp_error_low_fit = new TF1("tf1_kelly_GMp_error_low", kelly_GMp_error_low_fit, GMp_min_x, GMp_max_x, 0);
		tf1_kelly_GMp_error_low_fit->SetNpx(GMp_fit_Npx);
		tf1_kelly_GMp_error_low_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GMp_error_low_fit->SetLineWidth(1);
		tf1_kelly_GMp_error_low_fit->SetLineStyle(7);
		tf1_kelly_GMp_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GMp, h_GMp, tf1_kelly_GMp_error_high_fit, tf1_kelly_GMp_error_low_fit, tf1_kelly_GMp, GMp_min_y, GMp_max_y, kelly_line_color, kelly_shade_color);
	
	}

	if( FF == 0 || FF == 3 ){
		double GEn_min_x = 0.003;
		double GEn_max_x = 10;
		double GEn_min_y = -0.1;
		double GEn_max_y = 0.55;
		double GEn_fit_Npx = 5000;

		c_overlay_kelly_ye_GEn->cd();
		TH1D *h_GEn = new TH1D("h_GEn", "GEn - Kelly Parameterization", GEn_fit_Npx, GEn_min_x, GEn_max_x);
		h_GEn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEn->GetXaxis()->SetTitleOffset(1.25f);

		h_GEn->GetYaxis()->SetTitle("G_{E}^{n}/G_{D}");
		h_GEn->GetYaxis()->SetTitleOffset(1.25f);
		h_GEn->SetLineColor(0);
		h_GEn->GetYaxis()->SetRangeUser(GEn_min_y, GEn_max_y);
		h_GEn->GetXaxis()->SetRangeUser(GEn_min_x, GEn_max_x);
		h_GEn->SetStats(0);
		h_GEn->Draw("same");

		c_overlay_kelly_ye_GEn->SetLogx();
		c_overlay_kelly_ye_GEn->Range(GEn_min_x, GEn_min_y, GEn_max_x, GEn_max_y);
		tf1_kelly_GEn = new TF1("tf1_kelly_GEn", kelly_GEn_parameterization, GEn_min_x, GEn_max_x, 0);
		tf1_kelly_GEn->SetNpx(GEn_fit_Npx);
		tf1_kelly_GEn->SetLineColor(kelly_line_color);
		tf1_kelly_GEn->SetLineWidth(2);
		tf1_kelly_GEn->Draw("same");

		tf1_kelly_GEn_error_high_fit = new TF1("tf1_kelly_GEn_error_high_fit", kelly_GEn_error_high_fit, GEn_min_x, GEn_max_x, 0);
		tf1_kelly_GEn_error_high_fit->SetNpx(GEn_fit_Npx);
		tf1_kelly_GEn_error_high_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GEn_error_high_fit->SetLineWidth(1);
		tf1_kelly_GEn_error_high_fit->SetLineStyle(7);
		tf1_kelly_GEn_error_high_fit->Draw("same");

		tf1_kelly_GEn_error_low_fit = new TF1("tf1_kelly_GEn_error_low", kelly_GEn_error_low_fit, GEn_min_x, GEn_max_x, 0);
		tf1_kelly_GEn_error_low_fit->SetNpx(GEn_fit_Npx);
		tf1_kelly_GEn_error_low_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GEn_error_low_fit->SetLineWidth(1);
		tf1_kelly_GEn_error_low_fit->SetLineStyle(7);
		tf1_kelly_GEn_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GEn, h_GEn, tf1_kelly_GEn_error_high_fit, tf1_kelly_GEn_error_low_fit, tf1_kelly_GEn, GEn_min_y, GEn_max_y, kelly_line_color, kelly_shade_color);
	
	}

	if( FF == 0 || FF == 4 ){
		double GMn_min_x = 0.01;
		double GMn_max_x = 20;
		double GMn_min_y = 0.4;
		double GMn_max_y = 1.15;
		double GMn_fit_Npx = 5000;

		c_overlay_kelly_ye_GMn->SetGrid();	
		c_overlay_kelly_ye_GMn->cd();
		TH1D *h_GMn = new TH1D("h_GMn", "GMn - Kelly  and Ye Parameterizations", GMn_fit_Npx, GMn_min_x, GMn_max_x);
		h_GMn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMn->GetXaxis()->SetTitleOffset(1.25f);

		h_GMn->GetYaxis()->SetTitle("G_{M}^{n}/(#mu_{n} G_{D})");
		h_GMn->GetYaxis()->SetTitleOffset(1.25f);
		h_GMn->SetLineColor(0);
		h_GMn->GetYaxis()->SetRangeUser(GMn_min_y, GMn_max_y);
		h_GMn->GetXaxis()->SetRangeUser(GMn_min_x, GMn_max_x);
		h_GMn->SetStats(0);
		h_GMn->Draw("same");

		c_overlay_kelly_ye_GMn->SetLogx();
		c_overlay_kelly_ye_GMn->Range(GMn_min_x, GMn_min_y, GMn_max_x, GMn_max_y);
		tf1_kelly_GMn = new TF1("tf1_kelly_GMn", kelly_GMn_parameterization, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn->SetLineColor(kelly_line_color);
		tf1_kelly_GMn->SetLineWidth(2);
		tf1_kelly_GMn->Draw("same");

		tf1_kelly_GMn_error_high_fit = new TF1("tf1_kelly_GMn_error_high_fit", kelly_GMn_error_high_fit, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn_error_high_fit->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn_error_high_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GMn_error_high_fit->SetLineWidth(1);
		tf1_kelly_GMn_error_high_fit->SetLineStyle(7);
		tf1_kelly_GMn_error_high_fit->Draw("same");

		tf1_kelly_GMn_error_low_fit = new TF1("tf1_kelly_GMn_error_low", kelly_GMn_error_low_fit, GMn_min_x, GMn_max_x, 0);
		tf1_kelly_GMn_error_low_fit->SetNpx(GMn_fit_Npx);
		tf1_kelly_GMn_error_low_fit->SetLineColor(kelly_line_color);
		tf1_kelly_GMn_error_low_fit->SetLineWidth(1);
		tf1_kelly_GMn_error_low_fit->SetLineStyle(7);
		tf1_kelly_GMn_error_low_fit->Draw("same");

		fshade_btw_TF1s_w_3rd_SAME(c_overlay_kelly_ye_GMn, h_GMn, tf1_kelly_GMn_error_high_fit, tf1_kelly_GMn_error_low_fit, tf1_kelly_GMn, GMn_min_y, GMn_max_y, kelly_line_color, kelly_shade_color);
	
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
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd, Form("SBS8 G_{M}^{n} = %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[0], jboyd_Q2_vec[0]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd, Form("SBS8 syst. err: %.4f", jboyd_GMn_syst_error_vec[0]));	
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd_stat, Form("SBS8 stat. err: %.4f", jboyd_GMn_stat_error_vec[0] - jboyd_GMn_syst_error_vec[0]));	
			tleg_GMn_jboyd->AddEntry((TObject*)0, "", "");

			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd, Form("SBS9 G_{M}^{n} = %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[1], jboyd_Q2_vec[1]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd, Form("SBS9 syst. err: %.4f ", jboyd_GMn_syst_error_vec[1]));
			tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd_stat, Form("SBS9 stat. err: %.4f", jboyd_GMn_stat_error_vec[1] - jboyd_GMn_syst_error_vec[1]));

			if( plot_jboyd_GMn4 ){
				tleg_GMn_jboyd->AddEntry((TObject*)0, "", "");
				tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 G_{M}^{n} = %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[2], jboyd_Q2_vec[2]));
				tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 syst. err: %.4f", jboyd_GMn_syst_error_vec[2]));
				tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 syst. err: %.4f", jboyd_GMn_stat_error_vec[2] - jboyd_GMn_syst_error_vec[2]));
			}	

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

		c_overlay_kelly_ye_GMn->cd();
		tge_GMn_SBS4_jboyd_stat->Draw("*same");
		tge_GMn_SBS4_jboyd->Draw("*same");
		tge_GMn_SBS4_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS4_jboyd->Draw("same");
		tg_GMn_SBS4_jboyd->Draw("*same");
		tg_GMn_SBS4_jboyd->SetLineColor(kelly_line_color);
		tg_GMn_SBS4_jboyd->SetMarkerColor(1);
		tg_GMn_SBS4_jboyd->SetMarkerStyle(kOpenCircle);
		tg_GMn_SBS4_jboyd->Draw("same");

		tge_GMn_SBS8_jboyd_stat->Draw("*same");
		tge_GMn_SBS8_jboyd->Draw("*same");
		tge_GMn_SBS8_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS8_jboyd->Draw("same");
		tg_GMn_SBS8_jboyd->Draw("*same");
		tg_GMn_SBS8_jboyd->SetLineColor(kelly_line_color);
		tg_GMn_SBS8_jboyd->SetMarkerColor(1);
		tg_GMn_SBS8_jboyd->SetMarkerStyle(kOpenCircle);
		tg_GMn_SBS8_jboyd->Draw("same");

		tge_GMn_SBS9_jboyd_stat->Draw("*same");
		tge_GMn_SBS9_jboyd->Draw("*same");
		tge_GMn_SBS9_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS9_jboyd->Draw("same");
		tg_GMn_SBS9_jboyd->Draw("*same");
		tg_GMn_SBS9_jboyd->SetLineColor(kelly_line_color);
		tg_GMn_SBS9_jboyd->SetMarkerColor(1);
		tg_GMn_SBS9_jboyd->SetMarkerStyle(kOpenCircle);
		tg_GMn_SBS9_jboyd->Draw("same");

	}

}

void plot_GMp(){
	c_overlay_kelly_ye_GMp = new TCanvas("c_overlay_kelly_ye_GMp", "c_overlay_kelly_ye_GMp", 600, 500);
	plot_kelly_FFs(2);
	plot_ye_arrington(2);
}

void plot_GEp(){
	c_overlay_kelly_ye_GEp = new TCanvas("c_overlay_kelly_ye_GEp", "c_overlay_kelly_ye_GEp", 600, 500);
	plot_kelly_FFs(1);
	plot_ye_arrington(1);
}

void plot_GMn(){
	c_overlay_kelly_ye_GMn = new TCanvas("c_overlay_kelly_ye_GMn", "c_overlay_kelly_ye_GMn", 600, 500);
	plot_kelly_FFs(4);
	plot_ye_arrington(4);

	TLegend *tl_GMn = new TLegend(0.125, 0.78, 0.25, 0.885);
	tl_GMn->AddEntry(tf1_kelly_GMn, "Kelly");
	tl_GMn->AddEntry(tf1_ye_arrington_GMn_fit, "Ye");
	tl_GMn->Draw("same");
}

void plot_GEn(){
	c_overlay_kelly_ye_GEn = new TCanvas("c_overlay_kelly_ye_GEn", "c_overlay_kelly_ye_GEn", 600, 500);
	plot_kelly_FFs(3);
	plot_ye_arrington(3);
}	

void plot_SBS4_GMn_line(bool plot_theory = true, bool plot_data = 1){
	c_overlay_kelly_ye_GMn->cd();
	TLine *tl_SBS4_GMn_theory = new TLine(SBS4_Q2_theory, 0.4, SBS4_Q2_theory, 1.5);
	tl_SBS4_GMn_theory->SetLineColor(kSpring + 1);
	tl_SBS4_GMn_theory->SetLineStyle(5);

	TLine *tl_SBS4_GMn_data = new TLine(jboyd_Q2_SBS4, 0.4, jboyd_Q2_SBS4, 1.5);
	tl_SBS4_GMn_data->SetLineColor(kOrange -2);
	tl_SBS4_GMn_data->SetLineStyle(6);

	if( plot_theory ){
		tl_SBS4_GMn_theory->Draw("same");
	}
	if( plot_data ){
		tl_SBS4_GMn_data->Draw("same");
	}
}
void plot_SBS8_GMn_line(bool plot_theory = true, bool plot_data = 1){
	c_overlay_kelly_ye_GMn->cd();
	TLine *tl_SBS8_GMn_theory = new TLine(SBS8_Q2_theory, 0.4, SBS8_Q2_theory, 1.5);
	tl_SBS8_GMn_theory->SetLineColor(kOrange + 8);
	tl_SBS8_GMn_theory->SetLineStyle(5);

	TLine *tl_SBS8_GMn_data = new TLine(jboyd_Q2_SBS8, 0.4, jboyd_Q2_SBS8, 1.5);
	tl_SBS8_GMn_data->SetLineColor(kViolet - 6);
	tl_SBS8_GMn_data->SetLineStyle(6);

	if( plot_theory ){
		tl_SBS8_GMn_theory->Draw("same");
	}
	if( plot_data ){
		tl_SBS8_GMn_data->Draw("same");
	}
}
void plot_SBS9_GMn_line(bool plot_theory = true, bool plot_data = 1){
	c_overlay_kelly_ye_GMn->cd();
	TLine *tl_SBS9_GMn_theory = new TLine(SBS9_Q2_theory, 0.4, SBS9_Q2_theory, 1.125);
	tl_SBS9_GMn_theory->SetLineColor(kAzure + 8);
	tl_SBS9_GMn_theory->SetLineStyle(5);

	TLine *tl_SBS9_GMn_data = new TLine(jboyd_Q2_SBS9, 0.4, jboyd_Q2_SBS9, 1.125);
	tl_SBS9_GMn_data->SetLineColor(kBlack + 1);
	tl_SBS9_GMn_data->SetLineStyle(6);

	if( plot_theory ){
		tl_SBS9_GMn_theory->Draw("same");
	}
	if( plot_data ){
		tl_SBS9_GMn_data->Draw("same");
	}
}

void draw_text_box(bool SBS4 = true, bool SBS8 = true,  bool SBS9 = true ){
	TPaveText *tpt_GMn = new TPaveText(0.125, 0.15, 0.6, 0.65, "NDC");
	c_overlay_kelly_ye_GMn->cd();
	if( SBS4 ){
		tpt_GMn->AddText("SBS4: ");
		tpt_GMn->AddText(Form("Nominal Q^{2}: %0.4f -- Data Q^{2}: %0.4f", SBS4_Q2_theory, jboyd_Q2_SBS4));
		tpt_GMn->AddText("- - - - - ");
		tpt_GMn->AddText(Form("Nominal Kelly: %0.4f -- Nominal Ye: %0.4f", tf1_kelly_GMn->Eval(SBS4_Q2_theory), tf1_ye_arrington_GMn_fit->Eval(SBS4_Q2_theory)));
		tpt_GMn->AddText("-----------------------------");
	}

	if( SBS8 ){
		tpt_GMn->AddText("SBS8: ");
		tpt_GMn->AddText(Form("Nominal Q^{2}: %0.4f -- Data Q^{2}: %0.4f", SBS8_Q2_theory, jboyd_Q2_SBS8));
		tpt_GMn->AddText("- - - - - ");
		tpt_GMn->AddText(Form("Nominal Kelly: %0.4f -- Nominal Ye: %0.4f", tf1_kelly_GMn->Eval(SBS8_Q2_theory), tf1_ye_arrington_GMn_fit->Eval(SBS8_Q2_theory)));
		tpt_GMn->AddText("-----------------------------");
	}

	if( SBS9 ){
		tpt_GMn->AddText("SBS9: ");
		tpt_GMn->AddText(Form("Nominal Q^{2}: %0.4f -- Data Q^{2}: %0.4f", SBS9_Q2_theory, jboyd_Q2_SBS9));
		tpt_GMn->AddText("- - - - - ");
		tpt_GMn->AddText(Form("Nominal Kelly: %0.4f -- Nominal Ye: %0.4f", tf1_kelly_GMn->Eval(SBS9_Q2_theory), tf1_ye_arrington_GMn_fit->Eval(SBS9_Q2_theory)));
		tpt_GMn->AddText("-----------------------------");
	}
	tpt_GMn->Draw("same");
}

void plot_jboyd_points(bool SBS4 = true, bool SBS8 = true, bool SBS9 = true){
		if( SBS4 ){
			plot_jboyd_data( true );
		}
		else{
			plot_jboyd_data();		
		}


		TString drawStyle = "*same";


		tge_GMn_SBS8_jboyd_stat->Draw("*same");
		tge_GMn_SBS9_jboyd_stat->Draw(drawStyle.Data());
		tge_GMn_SBS8_jboyd->Draw(drawStyle.Data());
		tge_GMn_SBS9_jboyd->Draw(drawStyle.Data());

		if( SBS4 ){
			tge_GMn_SBS4_jboyd->Draw(drawStyle.Data());				
		}


		TLegend *tleg_GMn_jboyd;
		tleg_GMn_jboyd = new TLegend(0.125, 0.15, 0.6, 0.65);
		tleg_GMn_jboyd->SetTextSize(0.03);
		if( SBS4 ){
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


		c_overlay_kelly_ye_GMn->cd();
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

void overlay_kelly_ye(){
	
}

#endif