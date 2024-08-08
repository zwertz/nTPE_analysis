#ifndef NTPE_PROPOSAL_FIG2_LT_PT_OVERLAY_H
#define NTPE_PROPOSAL_FIG2_LT_PT_OVERLAY_H

#include "/w/halla-scshelf2102/sbs/ewertz/nTPE_analysis/world_data/jboyd_data_points.h"


auto mg_nTPE_LT_PT_overlay = new TMultiGraph();
TGraphErrors *tg_nTPE_proposal_fig2_LT;
TGraphErrors *tg_nTPE_proposal_fig2_PT;

TGraphErrors *tg_jboyd_LT_SBS8;
TGraphErrors *tg_jboyd_LT_SBS9;
TLegend *tl_jboyd_LT = new TLegend(0.50, 0.15, 0.80, 0.35, "Experimental RS values:");

vector<double> jboyd_nTPE_LTRS_SBS8;
vector<double> jboyd_nTPE_LTRS_SBS9;

double SBS8_RS, SBS9_RS;
double SBS8_RS_err, SBS9_RS_err;
double SBS8_tau_n, SBS9_tau_n;
double SBS8_sigma_L, SBS8_sigma_T, SBS9_sigma_L, SBS9_sigma_T;

double SBS8_nTPE_GEn, SBS9_nTPE_GEn;

vector<vector<double>> nTPE_proposal_fig2_red_LT_PT = {
	{0.9033816425120773, 0.17794612794612796, 0.0, 0.0},
	{2.9021739130434785, 0.2583333333333333, 0.0, 0.0},
	{5.903381642512077, 0.2574915824915825, 0.0, 0.0},
};

vector<vector<double>> nTPE_proposal_fig2_green_LT_PT = {
	{1.1026570048309179, 0.17752525252525253, 0.0, 0.0},
	{3.101449275362319, 0.2617003367003367, 0.0, 0.0},
	{6.102657004830918, 0.2692760942760943, 0.0, 0.0},
};

vector<vector<double>> nTPE_proposal_fig2_blue_LT_PT = {
	{1.000000000000000, 0.18047138047138045, 0.0, 0.0},
	{3.004830917874396, 0.2684343434343434, 0.0, 0.0},
	{6.000000000000000, 0.2814814814814815, 0.0, 0.0},
};

void plot_nTPE_proposal_fig2_LT_PT_overlay(bool plot_jboyd = false, bool plot_exp_info = false ){

	SBS8_tau_n = jboyd_Q2_SBS8/(4.0*physics_constants::M_n*physics_constants::M_n);
	SBS9_tau_n = jboyd_Q2_SBS9/(4.0*physics_constants::M_n*physics_constants::M_n);

	SBS8_nTPE_GEn = 0.4137855; //ye: 0.515279, kelly: 0.4137855;
	SBS9_nTPE_GEn = 0.4142309; //ye: 0.516186, kelly: 0.4142309;

	SBS8_sigma_L = SBS8_nTPE_GEn*SBS8_nTPE_GEn;
	SBS9_sigma_L = SBS9_nTPE_GEn*SBS9_nTPE_GEn;

	SBS8_sigma_T = SBS8_tau_n*jboyd_GMn_SBS8*jboyd_GMn_SBS8;
	SBS9_sigma_T = SBS9_tau_n*jboyd_GMn_SBS9*jboyd_GMn_SBS9;

	// SBS8_RS	= 2.0*SBS8_sigma_L/SBS8_sigma_T;
	SBS8_RS = 0.417606402;
	// SBS8_RS = 0.3435;
	SBS9_RS	= 2.0*SBS9_sigma_L/SBS9_sigma_T;

	SBS8_RS_err = 0.051561254;
	SBS9_RS_err	= 0.07737;

	vector<double> jboyd_nTPE_LTRS_SBS8 = {jboyd_Q2_SBS8, SBS8_RS, 0.0, SBS8_RS_err};
	vector<double> jboyd_nTPE_LTRS_SBS9 = {jboyd_Q2_SBS9, SBS9_RS, 0.0, SBS9_RS_err};

	TCanvas *c_ntpe_LT_PT_overlay = new TCanvas("c_ntpe_LT_PT_overlay", "c_ntpe_LT_PT_overlay", 700, 500);
	const int n_LT = nTPE_proposal_fig2_green.size();
	const int n_PT = nTPE_proposal_fig2_green_LT_PT.size();
	const int n_jboyd_LTRS = 1;

	double eps_LT[n_LT], eps_LT_err[n_LT], ratio_LT[n_LT], ratio_LT_err[n_LT];
	double eps_PT[n_PT], eps_PT_err[n_PT], ratio_PT[n_PT], ratio_PT_err[n_PT];

	if( plot_jboyd ){
		tg_jboyd_LT_SBS8 = new TGraphErrors(1, &jboyd_nTPE_LTRS_SBS8[0], &jboyd_nTPE_LTRS_SBS8[1], &jboyd_nTPE_LTRS_SBS8[2], &jboyd_nTPE_LTRS_SBS8[3]);
		tg_jboyd_LT_SBS8->SetLineColor(kMagenta);
		tg_jboyd_LT_SBS8->SetMarkerColor(kMagenta);

		// tg_jboyd_LT_SBS9 = new TGraphErrors(1, &jboyd_nTPE_LTRS_SBS9[0], &jboyd_nTPE_LTRS_SBS9[1], &jboyd_nTPE_LTRS_SBS9[2], &jboyd_nTPE_LTRS_SBS9[3]);
		// tg_jboyd_LT_SBS9->SetLineColor(kViolet -1);
		// tg_jboyd_LT_SBS9->SetMarkerColor(kViolet -1);

		// tl_jboyd_LT->AddEntry(tg_jboyd_LT_SBS8, Form("SBS8: %.4f +/- %0.4f", jboyd_nTPE_LTRS_SBS8[1], SBS8_RS_err));
		// tl_jboyd_LT->AddEntry(tg_jboyd_LT_SBS9, Form("SBS9: %.4f +/- %0.4f", jboyd_nTPE_LTRS_SBS9[1], SBS9_RS_err));
		tl_jboyd_LT->AddEntry(tg_jboyd_LT_SBS8, Form("SBS 8&9: %.4f +/- %0.4f", jboyd_nTPE_LTRS_SBS8[1], SBS8_RS_err));

		tl_jboyd_LT->SetTextSize(0.04f);
	}

	for(int i = 0; i < n_LT; i++ ){
		eps_LT[i] = nTPE_proposal_fig2_green[i][0];
		ratio_LT[i] = nTPE_proposal_fig2_green[i][1];
		eps_LT_err[i] = nTPE_proposal_fig2_green[i][2];
		ratio_LT_err[i] = nTPE_proposal_fig2_green[i][3];
	}

	for(int i = 0; i < n_PT; i++ ){
		eps_PT[i] = nTPE_proposal_fig2_green_LT_PT[i][0];
		ratio_PT[i] = nTPE_proposal_fig2_green_LT_PT[i][1];
		eps_PT_err[i] = nTPE_proposal_fig2_green_LT_PT[i][2];
		ratio_PT_err[i] = nTPE_proposal_fig2_green_LT_PT[i][3];
	}

	tg_nTPE_proposal_fig2_LT = new TGraphErrors(n_LT, eps_LT, ratio_LT, eps_LT_err, ratio_LT_err);
	tg_nTPE_proposal_fig2_LT->SetMarkerColor(kRed);
	tg_nTPE_proposal_fig2_LT->SetLineColor(kRed);
	tg_nTPE_proposal_fig2_LT->SetLineStyle(7);
	tg_nTPE_proposal_fig2_LT->SetMarkerStyle(21);
	tg_nTPE_proposal_fig2_LT->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_LT->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_LT->SetTitle("L-T Separation - #mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_LT->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_LT->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_LT; i++ ){
		tg_nTPE_proposal_fig2_LT->SetPoint(i, eps_LT[i], ratio_LT[i]);
		tg_nTPE_proposal_fig2_LT->SetPointError(i, eps_LT_err[i], ratio_LT_err[i]);
	}

	tg_nTPE_proposal_fig2_PT = new TGraphErrors(n_PT, eps_PT, ratio_PT, eps_PT_err, ratio_PT_err);
	tg_nTPE_proposal_fig2_PT->SetMarkerColor(kBlue);
	tg_nTPE_proposal_fig2_PT->SetLineColor(kBlue);
	tg_nTPE_proposal_fig2_PT->SetLineStyle(3);
	tg_nTPE_proposal_fig2_PT->SetMarkerStyle(8);
	tg_nTPE_proposal_fig2_PT->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_PT->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_PT->SetTitle("Polarization Transfer - #mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_PT->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_PT->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_PT; i++ ){
		tg_nTPE_proposal_fig2_PT->SetPoint(i, eps_PT[i], ratio_PT[i]);
		tg_nTPE_proposal_fig2_PT->SetPointError(i, eps_PT_err[i], ratio_PT_err[i]);
	}

	mg_nTPE_LT_PT_overlay->SetMinimum(0.0);
	mg_nTPE_LT_PT_overlay->SetMaximum(0.5);
	mg_nTPE_LT_PT_overlay->SetTitle("LT & PT Methods - Neutron FF Ratio vs Momentum Transfer, Q^{2}; Q^{2} (GeV/c)^{2}; #propto_{n} G_{E}^{n}/G_{M}^{n}");
	mg_nTPE_LT_PT_overlay->Add(tg_nTPE_proposal_fig2_LT);
	mg_nTPE_LT_PT_overlay->Add(tg_nTPE_proposal_fig2_PT);

	mg_nTPE_LT_PT_overlay->Draw("APL");

	TLegend *tl_nTPE_LT_PT_overlay = new TLegend(0.12, 0.70, 0.52, 0.90);
	tl_nTPE_LT_PT_overlay->AddEntry(tg_nTPE_proposal_fig2_LT, "LT Separation");
	tl_nTPE_LT_PT_overlay->AddEntry(tg_nTPE_proposal_fig2_PT, "Polarization Transfer");
	tl_nTPE_LT_PT_overlay->Draw("same");

	if( plot_jboyd ){
		tg_jboyd_LT_SBS8->Draw("sameP*");
		// tg_jboyd_LT_SBS9->Draw("sameP*");
		tl_jboyd_LT->Draw("same");
	}

	if( plot_exp_info ){

		TF1 *tf_nTPE_LT_simple_linear = new TF1("tf_nTPE_LT_simple_linear", "pol1", eps_LT[1], eps_LT[2] );
		tg_nTPE_proposal_fig2_LT->Fit("tf_nTPE_LT_simple_linear", "R+0");

		TF1 *tf_nTPE_PT_simple_linear = new TF1("tf_nTPE_PT_simple_linear", "pol1", eps_PT[1], eps_PT[2] );
		tg_nTPE_proposal_fig2_PT->Fit("tf_nTPE_PT_simple_linear", "R+0");

		TLine *tline_LT_PT_overlay = new TLine(4.5, 0, 4.5, 0.5);
		tline_LT_PT_overlay->SetLineStyle(10);
		tline_LT_PT_overlay->Draw("same");

		TLatex txt_overlay_Q2;
		txt_overlay_Q2.SetTextSize(0.035f);
		txt_overlay_Q2.DrawLatex(4.54, 0.47, Form("ntpe, Q^{2} = 4.5 (GeV/c)^{2}"));
		txt_overlay_Q2.Draw("same");

		TLatex txt_overlay_LT;
		txt_overlay_LT.SetTextColor(kRed);
		txt_overlay_LT.SetTextSize(0.035f);
		txt_overlay_LT.DrawLatex(4.54, 0.35, Form("LT Separation --> %0.3f", tf_nTPE_LT_simple_linear->Eval(4.5)));
		txt_overlay_LT.Draw("same");

		TLatex txt_overlay_PT;
		txt_overlay_PT.SetTextColor(kBlue);
		txt_overlay_PT.SetTextSize(0.035f);
		txt_overlay_PT.DrawLatex(4.54, 0.23, Form("Pol. Transfer --> %0.3f", tf_nTPE_PT_simple_linear->Eval(4.5)));
		txt_overlay_PT.Draw("same");

		// tf_nTPE_LT_simple_linear->Draw("same");
		// tf_nTPE_PT_simple_linear->Draw("same");

	}

//-----------------------------------------------------------------------------
//Plot it with polynomial fits

	TCanvas *c_ntpe_LT_PT_overlay_fits = new TCanvas("c_ntpe_LT_PT_overlay_fits", "c_ntpe_LT_PT_overlay_fits", 700, 500);

	TF1 *tf_nTPE_proposal_fig2_LT = new TF1("tf_nTPE_proposal_fig2_LT", "pol2", 0, 6.5);
	tg_nTPE_proposal_fig2_LT->Fit("tf_nTPE_proposal_fig2_LT", "R+0");
	tf_nTPE_proposal_fig2_LT->SetLineStyle(7);
	tf_nTPE_proposal_fig2_LT->SetLineColor(kRed);

	TF1 *tf_nTPE_proposal_fig2_PT = new TF1("tf_nTPE_proposal_fig2_PT", "pol2", 0, 6.5);
	tg_nTPE_proposal_fig2_PT->Fit("tf_nTPE_proposal_fig2_PT", "R+0");
	tf_nTPE_proposal_fig2_PT->SetLineStyle(3);
	tf_nTPE_proposal_fig2_PT->SetLineColor(kBlue);	

	mg_nTPE_LT_PT_overlay->Draw("AP");
	tf_nTPE_proposal_fig2_LT->Draw("same");
	tf_nTPE_proposal_fig2_PT->Draw("same");

	tl_nTPE_LT_PT_overlay->Draw("same");

	if( plot_jboyd ){
		tg_jboyd_LT_SBS8->Draw("sameP*");
		// tg_jboyd_LT_SBS9->Draw("sameP*");
		tl_jboyd_LT->Draw("same");
	}

	if( plot_exp_info ){
		TLine *tline_LT_PT_overlay = new TLine(4.5, 0, 4.5, 0.5);
		tline_LT_PT_overlay->SetLineStyle(10);
		tline_LT_PT_overlay->Draw("same");

		TLatex txt_overlay_Q2;
		txt_overlay_Q2.SetTextSize(0.035f);
		txt_overlay_Q2.DrawLatex(4.54, 0.47, Form("ntpe, Q^{2} = 4.5 (GeV/c)^{2}"));
		txt_overlay_Q2.Draw("same");

		TLatex txt_overlay_LT;
		txt_overlay_LT.SetTextColor(kBlue);
		txt_overlay_LT.SetTextSize(0.035f);
		txt_overlay_LT.DrawLatex(4.54, 0.35, Form("LT Separation --> %0.3f", tf_nTPE_proposal_fig2_LT->Eval(4.5)));
		txt_overlay_LT.Draw("same");

		TLatex txt_overlay_PT;
		txt_overlay_PT.SetTextColor(kViolet);
		txt_overlay_PT.SetTextSize(0.035f);
		txt_overlay_PT.DrawLatex(4.54, 0.23, Form("Pol. Transfer --> %0.3f", tf_nTPE_proposal_fig2_PT->Eval(4.5)));
		txt_overlay_PT.Draw("same");

	}
}

#endif
