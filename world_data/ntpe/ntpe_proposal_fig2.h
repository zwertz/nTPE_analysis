#ifndef NTPE_PROPOSAL_FIG2_H
#define NTPE_PROPOSAL_FIG2_H

vector<vector<double>> nTPE_proposal_fig2_red = {
	{0.9551724137931035, 0.19609164420485176, 0.0, 0.0},
	{2.9492610837438424, 0.33409703504043126, 0.0, 0.0},
	{5.951231527093596, 0.4553908355795149, 0.0, 0.0},
};

vector<vector<double>> nTPE_proposal_fig2_green = {
	{1.0527093596059114, 0.19609164420485176, 0.0, 0.0},
	{3.0522167487684726, 0.32601078167115904, 0.0, 0.0},
	{6.054187192118226, 0.41549865229110516, 0.0, 0.0},
};

vector<vector<double>> nTPE_proposal_fig2_blue = {
	{1.0039408866995072, 0.1804582210242588, 0.0, 0.0},
	{2.998029556650246, 0.2688679245283019, 0.0, 0.0},
	{6.000000000000000, 0.281266846361186, 0.0, 0.0},
};
auto *mg_nTPE_proposal_fig2_LT = new TMultiGraph();

void plot_nTPE_proposal_fig2(bool plot_jboyd = true, bool overlay_LT_PT = true ){

	TCanvas *c_GMN_proposal_data = new TCanvas("c_GMN_proposal_data", "GMN_proposal_data", 700, 500);
	const int n_red = nTPE_proposal_fig2_red.size();
	const int n_green = nTPE_proposal_fig2_green.size();
	const int n_blue = nTPE_proposal_fig2_blue.size();

	double eps_red[n_red], eps_red_err[n_red], ratio_red[n_red], ratio_red_err[n_red];
	double eps_green[n_green], eps_green_err[n_green], ratio_green[n_green], ratio_green_err[n_green];
	double eps_blue[n_blue], eps_blue_err[n_blue], ratio_blue[n_blue], ratio_blue_err[n_blue];

	for(int i = 0; i < n_red; i++ ){
		eps_red[i] = nTPE_proposal_fig2_red[i][0];
		ratio_red[i] = nTPE_proposal_fig2_red[i][1];
		eps_red_err[i] = nTPE_proposal_fig2_red[i][2];
		ratio_red_err[i] = nTPE_proposal_fig2_red[i][3];
	}

	for(int i = 0; i < n_green; i++ ){
		eps_green[i] = nTPE_proposal_fig2_green[i][0];
		ratio_green[i] = nTPE_proposal_fig2_green[i][1];
		eps_green_err[i] = nTPE_proposal_fig2_green[i][2];
		ratio_green_err[i] = nTPE_proposal_fig2_green[i][3];
	}

	for(int i = 0; i < n_blue; i++ ){
		eps_blue[i] = nTPE_proposal_fig2_blue[i][0];
		ratio_blue[i] = nTPE_proposal_fig2_blue[i][1];
		eps_blue_err[i] = nTPE_proposal_fig2_blue[i][2];
		ratio_blue_err[i] = nTPE_proposal_fig2_blue[i][3];
	}

	TGraphErrors *tg_nTPE_proposal_fig2_red = new TGraphErrors(n_red, eps_red, ratio_red, eps_red_err, ratio_red_err);
	tg_nTPE_proposal_fig2_red->SetMarkerColor(kRed);
	tg_nTPE_proposal_fig2_red->SetLineColor(kRed);
	tg_nTPE_proposal_fig2_red->SetMarkerStyle(21);
	tg_nTPE_proposal_fig2_red->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_red->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_red->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_red->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_red->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_red; i++ ){
		tg_nTPE_proposal_fig2_red->SetPoint(i, eps_red[i], ratio_red[i]);
		tg_nTPE_proposal_fig2_red->SetPointError(i, eps_red_err[i], ratio_red_err[i]);
	}

	TGraphErrors *tg_nTPE_proposal_fig2_green = new TGraphErrors(n_green, eps_green, ratio_green, eps_green_err, ratio_green_err);
	tg_nTPE_proposal_fig2_green->SetMarkerColor(kGreen);
	tg_nTPE_proposal_fig2_green->SetLineColor(kGreen);
	tg_nTPE_proposal_fig2_green->SetMarkerStyle(8);
	tg_nTPE_proposal_fig2_green->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_green->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_green->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_green->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_green->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_green; i++ ){
		tg_nTPE_proposal_fig2_green->SetPoint(i, eps_green[i], ratio_green[i]);
		tg_nTPE_proposal_fig2_green->SetPointError(i, eps_green_err[i], ratio_green_err[i]);
	}

	TGraphErrors *tg_nTPE_proposal_fig2_blue = new TGraphErrors(n_blue, eps_blue, ratio_blue, eps_blue_err, ratio_blue_err);
	tg_nTPE_proposal_fig2_blue->SetMarkerColor(kBlue);
	tg_nTPE_proposal_fig2_blue->SetLineColor(kBlue);
	tg_nTPE_proposal_fig2_blue->SetMarkerStyle(4);
	tg_nTPE_proposal_fig2_blue->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_blue->SetLineWidth(6);
	tg_nTPE_proposal_fig2_blue->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_blue->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_blue->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_blue->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_blue; i++ ){
		tg_nTPE_proposal_fig2_blue->SetPoint(i, eps_blue[i], ratio_blue[i]);
		tg_nTPE_proposal_fig2_blue->SetPointError(i, eps_blue_err[i], ratio_blue_err[i]);
	}


	mg_nTPE_proposal_fig2_LT->SetMinimum(0.0);
	mg_nTPE_proposal_fig2_LT->SetMaximum(0.5);
	mg_nTPE_proposal_fig2_LT->SetTitle("LT Sep. NTPE - Neutron FF Ratio vs Momentum Transfer, Q^{2}; Q^{2} (GeV/c)^{2}; G_{E}^{n}/G_{M}^{n}");
	mg_nTPE_proposal_fig2_LT->Add(tg_nTPE_proposal_fig2_red);
	mg_nTPE_proposal_fig2_LT->Add(tg_nTPE_proposal_fig2_green);
	mg_nTPE_proposal_fig2_LT->Add(tg_nTPE_proposal_fig2_blue);


	if( plot_jboyd ){

		cout << "Plotting JBoyd's GMn points...." << endl;
		double jboyd_SBS8 = 0.349320738; //0.844, //0.859
		double jboyd_SBS9 = 0.348481676;
								//SBS9          SBS8
		double ratio_jboyd[2] = {0.348481676, 0.349320738}; //0.86709 //0.844 //0.8508
		double eps_jboyd[2] = {4.504703224, 4.503761453};
		const int n_jboyd = 2;
		double eps_jboyd_err[n_jboyd], ratio_jboyd_err[n_jboyd];
		TGraphErrors *tg_nTPE_jboyd = new TGraphErrors(n_jboyd, eps_jboyd, ratio_jboyd, eps_jboyd_err, ratio_jboyd_err);

		tg_nTPE_jboyd->SetPoint(0, eps_jboyd[0], ratio_jboyd[0]);
		tg_nTPE_jboyd->SetMarkerStyle(5);
		tg_nTPE_jboyd->SetMarkerSize(2);
		tg_nTPE_jboyd->SetMarkerColor(kMagenta);
		tg_nTPE_jboyd->SetLineColor(kMagenta);

		// plot_GMn_proposal_data();
		// c_GMN_proposal_data->cd()
		// TCanvas *c = new TCanvas("c", "c", 600, 500);
		// tg_GMn_jboyd->Draw("AP+same");
		mg_nTPE_proposal_fig2_LT->Add(tg_nTPE_jboyd);
	}

	mg_nTPE_proposal_fig2_LT->Draw("AP");
	mg_nTPE_proposal_fig2_LT->GetXaxis()->SetLimits(0, 6.5);
	c_GMN_proposal_data->Modified();

	TF1 *fit_red, *fit_green, *fit_blue;
	if( true ){
		fit_red = new TF1("fit_red", "pol2", 0, 6.5);
		fit_red->SetLineWidth(1);
		fit_red->SetLineStyle(6);
		tg_nTPE_proposal_fig2_red->Fit("fit_red", "R+");

		TLine *tl_Q2 = new TLine(4.5, 0, 4.5, 0.5);
		tl_Q2->SetLineColor(kBlue);
		tl_Q2->SetLineStyle(10);
		tl_Q2->Draw("same");

		double eval_red = fit_red->Eval(4.5);

		TLatex txt_red;
		txt_red.SetTextColor(kRed);
		txt_red.SetTextSize(0.035f);
		txt_red.DrawLatex(4.56, 0.18, Form("0.2 < #epsilon < 0.9 --> %0.3f", eval_red));

		fit_green = new TF1("fit_green", "pol2", 0, 6.5);
		fit_green->SetLineColor(kGreen);
		fit_green->SetLineWidth(1);
		fit_green->SetLineStyle(7);
		tg_nTPE_proposal_fig2_green->Fit("fit_green", "R+");

		double eval_green = fit_green->Eval(4.5);

		TLatex txt_green;
		txt_green.SetTextColor(kGreen);
		txt_green.SetTextSize(0.035f);
		txt_green.DrawLatex(4.56, 0.15, Form("0.5 < #epsilon < 0.8 --> %0.3f", eval_green));

		fit_blue = new TF1("fit_blue", "pol2", 0, 6.5);
		fit_blue->SetLineColor(kBlue);
		fit_blue->SetLineWidth(1);
		fit_blue->SetLineStyle(8);
		tg_nTPE_proposal_fig2_blue->Fit("fit_blue", "R+");

		double eval_blue = fit_blue->Eval(4.5);

		TLatex txt_blue;
		txt_blue.SetTextColor(kBlue);
		txt_blue.SetTextSize(0.035f);
		txt_blue.DrawLatex(4.56, 0.12, Form("OPE --> %0.3f", eval_blue));

	}

}

#endif
