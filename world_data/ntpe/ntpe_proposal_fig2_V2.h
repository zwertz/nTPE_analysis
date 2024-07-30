#ifndef NTPE_PROPOSAL_FIG2_V2
#define NTPE_PROPOSAL_FIG2_V2


vector<vector<double>> nTPE_proposal_fig2_red_V2 = {
	{0.9033816425120773, 0.17794612794612796, 0.0, 0.0},
	{2.9021739130434785, 0.2583333333333333, 0.0, 0.0},
	{5.903381642512077, 0.2574915824915825, 0.0, 0.0},
};

vector<vector<double>> nTPE_proposal_fig2_green_V2 = {
	{1.1026570048309179, 0.17752525252525253, 0.0, 0.0},
	{3.101449275362319, 0.2617003367003367, 0.0, 0.0},
	{6.102657004830918, 0.2692760942760943, 0.0, 0.0},
};

vector<vector<double>> nTPE_proposal_fig2_blue_V2 = {
	{1.000000000000000, 0.18047138047138045, 0.0, 0.0},
	{3.004830917874396, 0.2684343434343434, 0.0, 0.0},
	{6.000000000000000, 0.2814814814814815, 0.0, 0.0},
};


auto *mg_nTPE_proposal_fig2_PT = new TMultiGraph();
void plot_nTPE_proposal_fig2_V2(bool plot_jboyd = false, bool overlay_LT_PT = false ){

	TCanvas *c_GMN_proposal_data = new TCanvas("c_GMN_proposal_data", "GMN_proposal_data", 700, 500);
	const int n_red = nTPE_proposal_fig2_red_V2.size();
	const int n_green = nTPE_proposal_fig2_green_V2.size();
	const int n_blue = nTPE_proposal_fig2_blue_V2.size();

	double eps_red[n_red], eps_red_err[n_red], ratio_red[n_red], ratio_red_err[n_red];
	double eps_green[n_green], eps_green_err[n_green], ratio_green[n_green], ratio_green_err[n_green];
	double eps_blue[n_blue], eps_blue_err[n_blue], ratio_blue[n_blue], ratio_blue_err[n_blue];

	for(int i = 0; i < n_red; i++ ){
		eps_red[i] = nTPE_proposal_fig2_red_V2[i][0];
		ratio_red[i] = nTPE_proposal_fig2_red_V2[i][1];
		eps_red_err[i] = nTPE_proposal_fig2_red_V2[i][2];
		ratio_red_err[i] = nTPE_proposal_fig2_red_V2[i][3];
	}

	for(int i = 0; i < n_green; i++ ){
		eps_green[i] = nTPE_proposal_fig2_green_V2[i][0];
		ratio_green[i] = nTPE_proposal_fig2_green_V2[i][1];
		eps_green_err[i] = nTPE_proposal_fig2_green_V2[i][2];
		ratio_green_err[i] = nTPE_proposal_fig2_green_V2[i][3];
	}

	for(int i = 0; i < n_blue; i++ ){
		eps_blue[i] = nTPE_proposal_fig2_blue_V2[i][0];
		ratio_blue[i] = nTPE_proposal_fig2_blue_V2[i][1];
		eps_blue_err[i] = nTPE_proposal_fig2_blue_V2[i][2];
		ratio_blue_err[i] = nTPE_proposal_fig2_blue_V2[i][3];
	}

	TGraphErrors *tg_nTPE_proposal_fig2_red_V2 = new TGraphErrors(n_red, eps_red, ratio_red, eps_red_err, ratio_red_err);
	tg_nTPE_proposal_fig2_red_V2->SetMarkerColor(kRed);
	tg_nTPE_proposal_fig2_red_V2->SetLineColor(kRed);
	tg_nTPE_proposal_fig2_red_V2->SetMarkerStyle(21);
	tg_nTPE_proposal_fig2_red_V2->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_red_V2->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_red_V2->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_red_V2->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_red_V2->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_red; i++ ){
		tg_nTPE_proposal_fig2_red_V2->SetPoint(i, eps_red[i], ratio_red[i]);
		tg_nTPE_proposal_fig2_red_V2->SetPointError(i, eps_red_err[i], ratio_red_err[i]);
	}

	TGraphErrors *tg_nTPE_proposal_fig2_green_V2 = new TGraphErrors(n_green, eps_green, ratio_green, eps_green_err, ratio_green_err);
	tg_nTPE_proposal_fig2_green_V2->SetMarkerColor(kGreen);
	tg_nTPE_proposal_fig2_green_V2->SetLineColor(kGreen);
	tg_nTPE_proposal_fig2_green_V2->SetMarkerStyle(8);
	tg_nTPE_proposal_fig2_green_V2->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_green_V2->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_green_V2->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_green_V2->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_green_V2->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_green; i++ ){
		tg_nTPE_proposal_fig2_green_V2->SetPoint(i, eps_green[i], ratio_green[i]);
		tg_nTPE_proposal_fig2_green_V2->SetPointError(i, eps_green_err[i], ratio_green_err[i]);
	}

	TGraphErrors *tg_nTPE_proposal_fig2_blue_V2 = new TGraphErrors(n_blue, eps_blue, ratio_blue, eps_blue_err, ratio_blue_err);
	tg_nTPE_proposal_fig2_blue_V2->SetMarkerColor(kBlue);
	tg_nTPE_proposal_fig2_blue_V2->SetLineColor(kBlue);
	tg_nTPE_proposal_fig2_blue_V2->SetMarkerStyle(4);
	tg_nTPE_proposal_fig2_blue_V2->SetMarkerSize(1);
	tg_nTPE_proposal_fig2_blue_V2->SetLineWidth(6);
	tg_nTPE_proposal_fig2_blue_V2->SetMinimum(0.0);
	tg_nTPE_proposal_fig2_blue_V2->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n} vs Momentum Transfer, Q^{2}");
	tg_nTPE_proposal_fig2_blue_V2->GetXaxis()->SetTitle("#mu_{n} G_{E}^{n}/G_{M}^{n}");
	tg_nTPE_proposal_fig2_blue_V2->GetYaxis()->SetTitle("#epsilon");

	for(int i = 0; i < n_blue; i++ ){
		tg_nTPE_proposal_fig2_blue_V2->SetPoint(i, eps_blue[i], ratio_blue[i]);
		tg_nTPE_proposal_fig2_blue_V2->SetPointError(i, eps_blue_err[i], ratio_blue_err[i]);
	}

	mg_nTPE_proposal_fig2_PT->SetMinimum(0.0);
	mg_nTPE_proposal_fig2_PT->SetMaximum(0.35);
	mg_nTPE_proposal_fig2_PT->SetTitle("Pol. Trans. TPE - Neutron FF Ratio vs Momentum Transfer, Q^{2}; Q^{2} (GeV/c)^{2}; G_{E}^{n}/G_{M}^{n}");
	mg_nTPE_proposal_fig2_PT->Add(tg_nTPE_proposal_fig2_red_V2);
	mg_nTPE_proposal_fig2_PT->Add(tg_nTPE_proposal_fig2_green_V2);
	mg_nTPE_proposal_fig2_PT->Add(tg_nTPE_proposal_fig2_blue_V2);


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
		mg_nTPE_proposal_fig2_PT->Add(tg_nTPE_jboyd);
	}

	mg_nTPE_proposal_fig2_PT->Draw("AP");
	mg_nTPE_proposal_fig2_PT->GetXaxis()->SetLimits(0, 6.5);
	c_GMN_proposal_data->Modified();

	TF1 *fit_red, *fit_green, *fit_blue;
	if( true ){
		fit_red = new TF1("fit_red", "pol2", 0, 6.5);
		fit_red->SetLineWidth(1);
		fit_red->SetLineStyle(6);
		tg_nTPE_proposal_fig2_red_V2->Fit("fit_red", "R+");

		TLine *tl_Q2 = new TLine(4.5, 0, 4.5, 0.35);
		tl_Q2->SetLineColor(kBlue);
		tl_Q2->SetLineStyle(10);
		tl_Q2->Draw("same");

		double eval_red = fit_red->Eval(4.5);

		TLatex txt_red;
		txt_red.SetTextColor(kRed);
		txt_red.SetTextSize(0.035f);
		txt_red.DrawLatex(4.54, 0.16, Form("#epsilon = 0.8 --> %0.3f", eval_red));

		fit_green = new TF1("fit_green", "pol2", 0, 6.5);
		fit_green->SetLineColor(kGreen);
		fit_green->SetLineWidth(1);
		fit_green->SetLineStyle(7);
		tg_nTPE_proposal_fig2_green_V2->Fit("fit_green", "R+");

		double eval_green = fit_green->Eval(4.5);

		TLatex txt_green;
		txt_green.SetTextColor(kGreen);
		txt_green.SetTextSize(0.035f);
		txt_green.DrawLatex(4.54, 0.18, Form("#epsilon = 0.3 --> %0.3f", eval_green));

		fit_blue = new TF1("fit_blue", "pol2", 0, 6.5);
		fit_blue->SetLineColor(kBlue);
		fit_blue->SetLineWidth(1);
		fit_blue->SetLineStyle(8);
		tg_nTPE_proposal_fig2_blue_V2->Fit("fit_blue", "R+");

		double eval_blue = fit_blue->Eval(4.5);

		TLatex txt_blue;
		txt_blue.SetTextColor(kBlue);
		txt_blue.SetTextSize(0.035f);
		txt_blue.DrawLatex(4.54, 0.2, Form("OPE --> %0.3f", eval_blue));
	}

}

#endif