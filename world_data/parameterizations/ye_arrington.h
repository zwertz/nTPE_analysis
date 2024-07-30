#ifndef YE_ARRINGTON_H
#define YE_ARRINGTON_H

#include "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/world_data/parameterizations/GetFFfunctions.h"
#include "/work/halla/sbs/jboyd/include/utility_functions.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/plot_pimping.h"
TF1 *tf1_ye_arrington_GEp_fit, *tf1_ye_arrington_GEp_fit_error_high, *tf1_ye_arrington_GEp_fit_error_low;
TF1 *tf1_ye_arrington_GMp_fit, *tf1_ye_arrington_GMp_fit_error_high, *tf1_ye_arrington_GMp_fit_error_low;
TF1 *tf1_ye_arrington_GEn_fit, *tf1_ye_arrington_GEn_fit_error_high, *tf1_ye_arrington_GEn_fit_error_low;
TF1 *tf1_ye_arrington_GMn_fit, *tf1_ye_arrington_GMn_fit_error_high, *tf1_ye_arrington_GMn_fit_error_low;

TGraph *tg_ye_arrington;
TGraphErrors *tg_GMn_exp_all;
TPolyLine *GMn_shaded_area;
auto *mg_ye_arrington = new TMultiGraph();

TGraph *graph = new TGraph();
TCanvas *c_ye_arrington_GMn;
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

	tf1_ye_arrington_GMp_fit = new TF1("tf1_ye_arrington_GMp_fit", GMp_parameterization, 0, 20, 11);
	tf1_ye_arrington_GEn_fit = new TF1("tf1_ye_arrington_GEn_fit", GEn_parameterization, 0, 20, 11);

	if( FF == 0 || FF == 1 ){
		double GEp_min_x = 0.01;
		double GEp_max_x = 100;
		double GEp_min_y = -1.0;
		double GEp_max_y = 1.1;
		double GEp_fit_Npx = 5000;

		TCanvas *c_ye_arrington_GEp = new TCanvas("c_ye_arrington_GEp", "c_ye_arrington_GEp", 600, 500);
		TH1D *h_GEp = new TH1D("h_GEp", "GEp - Ye Parameterization", GEp_fit_Npx, GEp_min_x, GEp_max_x);
		h_GEp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEp->GetXaxis()->SetTitleOffset(1.25f);

		h_GEp->GetYaxis()->SetTitle("G_{E}^{p}/G_{D}");
		h_GEp->GetYaxis()->SetTitleOffset(1.25f);
		h_GEp->SetLineColor(0);
		h_GEp->GetYaxis()->SetRangeUser(GEp_min_y, GEp_max_y);
		h_GEp->GetXaxis()->SetRangeUser(GEp_min_x, GEp_max_x);
		h_GEp->SetStats(0);
		h_GEp->Draw();

		c_ye_arrington_GEp->SetLogx();
		c_ye_arrington_GEp->Range(GEp_min_x, GEp_min_y, GEp_max_x, GEp_max_y);
		tf1_ye_arrington_GEp_fit = new TF1("tf1_ye_arrington_GEp_fit", GEp_parameterization, GEp_min_x, GEp_max_x, 0);
		tf1_ye_arrington_GEp_fit->SetNpx(GEp_fit_Npx);
		tf1_ye_arrington_GEp_fit->SetLineColor(1);
		tf1_ye_arrington_GEp_fit->SetLineWidth(2);
		tf1_ye_arrington_GEp_fit->Draw("same");

		tf1_ye_arrington_GEp_fit_error_high = new TF1("tf1_ye_arrington_GEp_fit_error_high", GEp_error_high_fit, GEp_min_x, GEp_max_x, 0);
		tf1_ye_arrington_GEp_fit_error_high->SetNpx(GEp_fit_Npx);
		tf1_ye_arrington_GEp_fit_error_high->SetLineColor(1);
		tf1_ye_arrington_GEp_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GEp_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GEp_fit_error_high->Draw("same");

		tf1_ye_arrington_GEp_fit_error_low = new TF1("tf1_ye_arrington_GEp_fit_error_low", GEp_error_low_fit, GEp_min_x, GEp_max_x, 0);
		tf1_ye_arrington_GEp_fit_error_low->SetNpx(GEp_fit_Npx);
		tf1_ye_arrington_GEp_fit_error_low->SetLineColor(1);
		tf1_ye_arrington_GEp_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GEp_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GEp_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd(c_ye_arrington_GEp, h_GEp, tf1_ye_arrington_GEp_fit_error_high, tf1_ye_arrington_GEp_fit_error_low, tf1_ye_arrington_GEp_fit, GEp_min_y, GEp_max_y);

	}

	if( FF == 0 || FF == 2 ){
		double GMp_min_x = 0.01;
		double GMp_max_x = 100;
		double GMp_min_y = 0.4;
		double GMp_max_y = 1.15;
		double GMp_fit_Npx = 5000;

		TCanvas *c_ye_arrington_GMp = new TCanvas("c_ye_arrington_GMp", "c_ye_arrington_GMp", 600, 500);
		TH1D *h_GMp = new TH1D("h_GMp", "GMp - Ye Parameterization", GMp_fit_Npx, GMp_min_x, GMp_max_x);
		h_GMp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMp->GetXaxis()->SetTitleOffset(1.25f);

		h_GMp->GetYaxis()->SetTitle("G_{M}^{p}/(#mu_{p} G_{D}");
		h_GMp->GetYaxis()->SetTitleOffset(1.25f);
		h_GMp->SetLineColor(0);
		h_GMp->GetYaxis()->SetRangeUser(GMp_min_y, GMp_max_y);
		h_GMp->GetXaxis()->SetRangeUser(GMp_min_x, GMp_max_x);
		h_GMp->SetStats(0);
		h_GMp->Draw();

		c_ye_arrington_GMp->SetLogx();
		c_ye_arrington_GMp->Range(GMp_min_x, GMp_min_y, GMp_max_x, GMp_max_y);
		tf1_ye_arrington_GMp_fit = new TF1("tf1_ye_arrington_GMp_fit", GMp_parameterization, GMp_min_x, GMp_max_x, 0);
		tf1_ye_arrington_GMp_fit->SetNpx(GMp_fit_Npx);
		tf1_ye_arrington_GMp_fit->SetLineColor(1);
		tf1_ye_arrington_GMp_fit->SetLineWidth(2);
		tf1_ye_arrington_GMp_fit->Draw("same");

		tf1_ye_arrington_GMp_fit_error_high = new TF1("tf1_ye_arrington_GMp_fit_error_high", GMp_error_high_fit, GMp_min_x, GMp_max_x, 0);
		tf1_ye_arrington_GMp_fit_error_high->SetNpx(GMp_fit_Npx);
		tf1_ye_arrington_GMp_fit_error_high->SetLineColor(1);
		tf1_ye_arrington_GMp_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GMp_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GMp_fit_error_high->Draw("same");

		tf1_ye_arrington_GMp_fit_error_low = new TF1("tf1_ye_arrington_GMp_fit_error_low", GMp_error_low_fit, GMp_min_x, GMp_max_x, 0);
		tf1_ye_arrington_GMp_fit_error_low->SetNpx(GMp_fit_Npx);
		tf1_ye_arrington_GMp_fit_error_low->SetLineColor(1);
		tf1_ye_arrington_GMp_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GMp_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GMp_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd(c_ye_arrington_GMp, h_GMp, tf1_ye_arrington_GMp_fit_error_high, tf1_ye_arrington_GMp_fit_error_low, tf1_ye_arrington_GMp_fit, GMp_min_y, GMp_max_y);

	}

	if( FF == 0 || FF == 3 ){
		double GEn_min_x = 0.01;
		double GEn_max_x = 20;
		double GEn_min_y = -0.025;
		double GEn_max_y = 0.8;
		double GEn_fit_Npx = 5000;

		TCanvas *c_ye_arrington_GEn = new TCanvas("c_ye_arrington_GEn", "c_ye_arrington_GEn", 600, 500);
		TH1D *h_GEn = new TH1D("h_GEn", "GEn - Ye Parameterization", GEn_fit_Npx, GEn_min_x, GEn_max_x);
		h_GEn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GEn->GetXaxis()->SetTitleOffset(1.25f);

		h_GEn->GetYaxis()->SetTitle("G_{E}^{n}/G_{D}");
		h_GEn->GetYaxis()->SetTitleOffset(1.25f);
		h_GEn->SetLineColor(0);
		h_GEn->GetYaxis()->SetRangeUser(GEn_min_y, GEn_max_y);
		h_GEn->GetXaxis()->SetRangeUser(GEn_min_x, GEn_max_x);
		h_GEn->SetStats(0);
		h_GEn->Draw();

		c_ye_arrington_GEn->SetLogx();
		c_ye_arrington_GEn->Range(GEn_min_x, GEn_min_y, GEn_max_x, GEn_max_y);
		tf1_ye_arrington_GEn_fit = new TF1("tf1_ye_arrington_GEn_fit", GEn_parameterization, GEn_min_x, GEn_max_x, 0);
		tf1_ye_arrington_GEn_fit->SetNpx(GEn_fit_Npx);
		tf1_ye_arrington_GEn_fit->SetLineColor(1);
		tf1_ye_arrington_GEn_fit->SetLineWidth(2);
		tf1_ye_arrington_GEn_fit->Draw("same");

		tf1_ye_arrington_GEn_fit_error_high = new TF1("tf1_ye_arrington_GEn_fit_error_high", GEn_error_high_fit, GEn_min_x, GEn_max_x, 0);
		tf1_ye_arrington_GEn_fit_error_high->SetNpx(GEn_fit_Npx);
		tf1_ye_arrington_GEn_fit_error_high->SetLineColor(1);
		tf1_ye_arrington_GEn_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GEn_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GEn_fit_error_high->Draw("same");

		tf1_ye_arrington_GEn_fit_error_low = new TF1("tf1_ye_arrington_GEn_fit_error_low", GEn_error_low_fit, GEn_min_x, GEn_max_x, 0);
		tf1_ye_arrington_GEn_fit_error_low->SetNpx(GEn_fit_Npx);
		tf1_ye_arrington_GEn_fit_error_low->SetLineColor(1);
		tf1_ye_arrington_GEn_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GEn_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GEn_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd(c_ye_arrington_GEn, h_GEn, tf1_ye_arrington_GEn_fit_error_high, tf1_ye_arrington_GEn_fit_error_low, tf1_ye_arrington_GEn_fit, GEn_min_y, GEn_max_y);

	}

	if( FF == 0 || FF == 4 ){
		double GMn_min_x = 0.01;
		double GMn_max_x = 20;
		double GMn_min_y = 0.175;
		double GMn_max_y = 1.2;
		double GMn_fit_Npx = 5000;

		c_ye_arrington_GMn = new TCanvas("c_ye_arrington_GMn", "c_ye_arrington_GMn", 600, 500);
		TH1D *h_GMn = new TH1D("h_GMn", "GMn - Ye Parameterization", GMn_fit_Npx, GMn_min_x, GMn_max_x);
		h_GMn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
		h_GMn->GetXaxis()->SetTitleOffset(1.25f);

		h_GMn->GetYaxis()->SetTitle("G_{M}^{n}/(#mu_{n} G_{D})");
		h_GMn->GetYaxis()->SetTitleOffset(1.25f);
		h_GMn->SetLineColor(0);
		h_GMn->GetYaxis()->SetRangeUser(GMn_min_y, GMn_max_y);
		h_GMn->GetXaxis()->SetRangeUser(GMn_min_x, GMn_max_x);
		h_GMn->SetStats(0);
		h_GMn->Draw();

		c_ye_arrington_GMn->SetLogx();
		c_ye_arrington_GMn->Range(GMn_min_x, GMn_min_y, GMn_max_x, GMn_max_y);
		tf1_ye_arrington_GMn_fit = new TF1("tf1_ye_arrington_GMn_fit", GMn_parameterization, GMn_min_x, GMn_max_x, 0);
		tf1_ye_arrington_GMn_fit->SetNpx(GMn_fit_Npx);
		tf1_ye_arrington_GMn_fit->SetLineColor(1);
		tf1_ye_arrington_GMn_fit->SetLineWidth(2);
		tf1_ye_arrington_GMn_fit->Draw("same");

		tf1_ye_arrington_GMn_fit_error_high = new TF1("tf1_ye_arrington_GMn_fit_error_high", GMn_error_high_fit, GMn_min_x, GMn_max_x, 0);
		tf1_ye_arrington_GMn_fit_error_high->SetNpx(GMn_fit_Npx);
		tf1_ye_arrington_GMn_fit_error_high->SetLineColor(1);
		tf1_ye_arrington_GMn_fit_error_high->SetLineWidth(1);
		tf1_ye_arrington_GMn_fit_error_high->SetLineStyle(7);
		tf1_ye_arrington_GMn_fit_error_high->Draw("same");

		tf1_ye_arrington_GMn_fit_error_low = new TF1("tf1_ye_arrington_GMn_fit_error_low", GMn_error_low_fit, GMn_min_x, GMn_max_x, 0);
		tf1_ye_arrington_GMn_fit_error_low->SetNpx(GMn_fit_Npx);
		tf1_ye_arrington_GMn_fit_error_low->SetLineColor(1);
		tf1_ye_arrington_GMn_fit_error_low->SetLineWidth(1);
		tf1_ye_arrington_GMn_fit_error_low->SetLineStyle(7);
		tf1_ye_arrington_GMn_fit_error_low->Draw("same");

		fshade_btw_TF1s_w_3rd(c_ye_arrington_GMn, h_GMn, tf1_ye_arrington_GMn_fit_error_high, tf1_ye_arrington_GMn_fit_error_low, tf1_ye_arrington_GMn_fit, GMn_min_y, GMn_max_y);

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

		// if( bool_plot_jboyd_GMn4 ){
		// 	tleg_GMn_jboyd->AddEntry(tge_GMn_SBS4_jboyd, Form("SBS4 G_{M}^{n} = %.4f +/- %.4f +/- %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[2], jboyd_GMn_syst_error_vec[2], jboyd_GMn_stat_error_vec[2] - jboyd_GMn_syst_error_vec[2], jboyd_Q2_vec[2]));
		// 	tleg_GMn_jboyd->AddEntry((TObject*)0, "", "");
		// }

		// tleg_GMn_jboyd->AddEntry(tge_GMn_SBS8_jboyd, Form("SBS8 G_{M}^{n} = %.4f +/- %.4f +/- %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[0], jboyd_GMn_syst_error_vec[0], jboyd_GMn_stat_error_vec[0] - jboyd_GMn_syst_error_vec[0], jboyd_Q2_vec[0]));
		// tleg_GMn_jboyd->AddEntry((TObject*)0, "", "");
		// tleg_GMn_jboyd->AddEntry(tge_GMn_SBS9_jboyd, Form("SBS9 G_{M}^{n} = %.4f +/- %.4f +/- %.4f, Q^{2} = %0.3f", jboyd_GMn_vec[1], jboyd_GMn_syst_error_vec[1], jboyd_GMn_stat_error_vec[1] - jboyd_GMn_syst_error_vec[1], jboyd_Q2_vec[1]));

		if( bool_plot_jboyd_GMn4 ){
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


		tleg_GMn_jboyd->Draw(drawStyle.Data());			

		c_ye_arrington_GMn->cd();
		tge_GMn_SBS4_jboyd_stat->Draw("*same");
		tge_GMn_SBS4_jboyd->Draw("*same");
		tge_GMn_SBS4_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS4_jboyd->Draw("same");
		tg_GMn_SBS4_jboyd->Draw("*same");
		tg_GMn_SBS4_jboyd->SetLineColor(1);
		tg_GMn_SBS4_jboyd->SetMarkerColor(1);
		tg_GMn_SBS4_jboyd->SetMarkerStyle(kFullCircle);
		tg_GMn_SBS4_jboyd->Draw("same");

		tge_GMn_SBS8_jboyd_stat->Draw("*same");
		tge_GMn_SBS8_jboyd->Draw("*same");
		tge_GMn_SBS8_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS8_jboyd->Draw("same");
		tg_GMn_SBS8_jboyd->Draw("*same");
		tg_GMn_SBS8_jboyd->SetLineColor(1);
		tg_GMn_SBS8_jboyd->SetMarkerColor(1);
		tg_GMn_SBS8_jboyd->SetMarkerStyle(kFullCircle);
		tg_GMn_SBS8_jboyd->Draw("same");

		tge_GMn_SBS9_jboyd_stat->Draw("*same");
		tge_GMn_SBS9_jboyd->Draw("*same");
		tge_GMn_SBS9_jboyd->SetMarkerStyle(kFullDotLarge);
		tge_GMn_SBS9_jboyd->Draw("same");
		tg_GMn_SBS9_jboyd->Draw("*same");
		tg_GMn_SBS9_jboyd->SetLineColor(1);
		tg_GMn_SBS9_jboyd->SetMarkerColor(1);
		tg_GMn_SBS9_jboyd->SetMarkerStyle(kFullCircle);
		tg_GMn_SBS9_jboyd->Draw("same");

	}

}

void general_plot_ye( TCanvas *c_to_plot_on, TLegend *tleg_to_add = NULL, int setLineColor = 2, int setLineStyle = 9, int setLineWidth = 2, bool plot_with_error = false){
	
	tf1_ye_arrington_GMp_fit = new TF1("tf1_ye_arrington_GMp_fit", GMp_parameterization, 0, 20, 11);
	tf1_ye_arrington_GEn_fit = new TF1("tf1_ye_arrington_GEn_fit", GEn_parameterization, 0, 20, 11);

	double GMn_min_x = 0.01;
	double GMn_max_x = 20;
	double GMn_min_y = 0.175;
	double GMn_max_y = 1.2;
	double GMn_fit_Npx = 5000;

	c_to_plot_on->cd();
	c_to_plot_on->SetGrid();

	TH1D *h_GMn = new TH1D("h_GMn", "GMn - Ye Parameterization", GMn_fit_Npx, GMn_min_x, GMn_max_x);
	h_GMn->Draw("same");

	tf1_ye_arrington_GMn_fit = new TF1("tf1_ye_arrington_GMn_fit", GMn_parameterization, GMn_min_x, GMn_max_x, 0);
	tf1_ye_arrington_GMn_fit->SetNpx(GMn_fit_Npx);
	tf1_ye_arrington_GMn_fit->SetLineColor(setLineColor);
	tf1_ye_arrington_GMn_fit->SetLineStyle(setLineStyle);
	tf1_ye_arrington_GMn_fit->SetLineWidth(2);
	tf1_ye_arrington_GMn_fit->Draw("same");

	if( tleg_to_add != NULL ){
		TString label_gmn_nominal = "Ye18";
		if( doesLegendEntryExist( tleg_to_add, label_gmn_nominal) ){

		}
		else{
			tleg_to_add->AddEntry(tf1_ye_arrington_GMn_fit, label_gmn_nominal.Data());
			c_to_plot_on->Update();
			c_to_plot_on->Modified();			
		}

	}

	// if( tleg_to_add != NULL ){
	// 	tleg_to_add->AddEntry(tf1_ye_arrington_GMn_fit, "Ye18");
	// 	c_to_plot_on->Update();
	// 	c_to_plot_on->Modified();
	// }

}

#endif