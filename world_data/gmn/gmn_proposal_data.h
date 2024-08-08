#ifndef GMN_PROPOSAL_DATA_H
#define GMN_PROPOSAL_DATA_H
#include "/w/halla-scshelf2102/sbs/ewertz/nTPE_analysis/world_data/jboyd_data_points.h"
TLegend *tleg_gmn_proposal_data;

//E_e, theta_e, Q^2, epsilon, sigma_r, sigma_r_erro
vector<vector<double>> kinem_and_reduced_cs_GMp12 = {
	{2.222,		42.001,		1.577,	0.701,	4.273E-2,	0.040E-2 },
	{2.222,		48.666,		1.858,	0.615,	2.983E-2,	0.057E-2 },
	{6.427,		24.250,		4.543,	0.826,	3.813E-3,	0.057E-3 },
	{6.427, 	30.909,		5.947,	0.709,	1.805E-3,	0.025E-3 },
	{6.427,		37.008,		6.993,	0.599,	1.113E-3,	0.016E-3 },
	{6.427,		44.500,		7.992,	0.478,	7.289E-4,	0.109E-4 },
	{8.518,		30.909,		9.002,	0.648,	5.163E-4,	0.078E-4 },
	{6.427,		55.900,		9.053,	0.332,	4.859E-4,	0.107E-4 },
	{8.518,		34.400,		9.807,	0.580,	3.923E-4,	0.059E-4 },
	{8.518,		42.001,		11.19,	0.448,	2.565E-4,	0.041E-4 },
	{8.518, 	48.666,		12.07,	0.356,	1.933E-4,	0.043E-4 },
	{8.518,		53.501,		12.57,	0.301,	1.664E-4,	0.053E-4 },
	{10.587,	48.666,		15.76,	0.309,	8.405E-4,	0.227E-4 }
};

//Green is SLAC [26]: S. Rock et. al., Phys. Rev. D 46 24 (1992).
vector<vector<double>> GMn_proposal_green = {
	{2.502673, 0.94289898, 0.0, 0.10295031/2},
	{3.957219251, 0.89370425, 0.0, 0.11304408/2},
	{6.010695187, 0.87525622, 0.0, 0.11708081/2},
	{8.021390374, 0.67847731, 0.0, 0.18773292/2},
	{10.01069519, 0.58828697, 0.0, 0.29673916/2},
};
// Green is SLAC [27]: A. Lung et al., Phys. Rev. Lett. 70, 718 (1993).
vector<vector<double>> GMn_proposal_magenta = {
	{1.689839572, 1.05153734, 0.0, 0.10295024/2},
	{2.459893048, 1.01669107, 0.0, 0.08680124/2},
	{3.251336898, 0.96544656, 0.0, 0.12111801/2},
	{4.021390374, 0.92240117, 0.0, 0.16149069/2},
};

//older.... [31] W. Bartel et al., Phys. Lett. B 30, 285 (1969); ibid. 39, 407 (1972); Nucl. Phys. B 58, 429 (1973).
vector<vector<double>> GMn_proposal_orange = {
	{1.0053475935828877, 1.0453879941434845, 0.0, 0.207919255/2},
	{1.1764705882352942, 1.0412884333821375, 0.0, 0.183695652/2},
	{1.5187165775401068, 1.0535871156661785, 0.0, 0.106987577/2},
	{1.7967914438502672, 1.2093704245973644, 0.0, 0.215993789/2},
	{2.6951871657754007, 0.8404099560761346, 0.0, 0.51878882/2},
};

vector<vector<double>> GMn_proposal_all_blacks = {
	{1.0481283422459893, 0.9756954612005855, 0.0, 0.15947205/2},
	{1.1122994652406417, 0.9101024890190336, 0.0, 0.078726708/2},
	{1.3048128342245988, 1.026939970717423, 0.0, 0.197826087/2},
	{1.6042780748663101, 0.936749633967789, 0.0, 0.276552795/2},
	{1.7540106951871657, 1.0474377745241579, 0.0, 0.127173904/2},
	{1.7326203208556148, 0.998243045387994, 0.0, 0.189347827/2},
};

//Blue is CLAS e5 [37] J, Lachniet thesis, Carnegie Mellon University, unpublished, June, 2005
// [39] J. Lachniet et al, nucl-ex 0811.1716, submitted to Phys. Rev. Lett.
vector<vector<double>> GMn_proposal_blue = {
	{0.5133689839572192, 0.9756954612005855, 0.0, 0.064596273/2},
	{0.6631016042780749, 0.9510980966325036, 0.0, 0.050465839/2},
	{0.8342245989304813, 0.973645680819912, 0.0, 0.048447205/2},
	{0.9839572192513368, 0.9961932650073205, 0.0, 0.044409938/2},
	{1.1336898395721924, 1.006442166910688, 0.0, 0.050465839/2},
	{1.2834224598930482, 1.03103953147877, 0.0, 0.042391304/2},
	{1.4331550802139037, 1.022840409956076, 0.0, 0.040372671/2},
	{1.5828877005347592, 1.0125915080527086, 0.0, 0.038354037/2},
	{1.7326203208556148, 1.0125915080527086, 0.0, 0.038354037/2},
	{1.8823529411764706, 1.0166910688140556, 0.0, 0.038354037/2},
	{2.032085561497326, 1.0002928257686674, 0.0, 0.042391304/2},
	{2.1818181818181817, 0.9900439238653, 0.0, 0.042391304/2},
	{2.3315508021390374, 1.002342606149341, 0.0, 0.044409938/2},
	{2.4812834224598928, 1.0207906295754026, 0.0, 0.046428571/2},
	{2.6524064171122994, 1.035139092240117, 0.0, 0.046428571/2},
	{2.8021390374331547, 1.0002928257686674, 0.0, 0.046428571/2},
	{2.9304812834224596, 0.9961932650073205, 0.0, 0.048447205/2},
	{3.1016042780748663, 1.002342606149341, 0.0, 0.048447205/2},
	{3.2513368983957216, 0.977745241581259, 0.0, 0.046428571/2},
	{3.4010695187165774, 0.994143484626647, 0.0, 0.050465839/2},
	{3.550802139037433, 0.9920937042459735, 0.0, 0.050465839/2},
	{3.7005347593582885, 0.973645680819912, 0.0, 0.052484472/2},
	{3.8502673796791442, 0.9797950219619325, 0.0, 0.054503106/2},
	{4.021390374331551, 1.0289897510980965, 0.0, 0.056521739/2},
	{4.149732620320855, 0.9183016105417275, 0.0, 0.058540373/2},
	{4.3208556149732615, 0.9080527086383601, 0.0, 0.064596273/2},
	{4.470588235294118, 0.9900439238653, 0.0, 0.082763975/2},
	{4.620320855614973, 1.1068814055636895, 0.0, 0.195807453/2},
	{4.770053475935828, 1.0904831625183014, 0.0, 0.258385093/2},
};

void plot_GMn_proposal_data(bool plot_jboyd = true, double custom_xmin = -999.9, double custom_xmax = -999.9, double custom_ymin = -999.9, double custom_ymax = -999.9 ){

	TCanvas *c_GMn_proposal_data = new TCanvas("c_GMn_proposal_data", "c_GMn_proposal_data", 600, 500);
	c_GMn_proposal_data->SetGrid();
	const int n_green = GMn_proposal_green.size();
	const int n_magenta = GMn_proposal_magenta.size();
	const int n_orange = GMn_proposal_orange.size();
	const int n_all_blacks = GMn_proposal_all_blacks.size();
	const int n_blue = GMn_proposal_blue.size();


	double Q2_green[n_green], Q2_green_err[n_green], GMn_green[n_green], GMn_green_err[n_green];
	double Q2_magenta[n_magenta], Q2_magenta_err[n_magenta], GMn_magenta[n_magenta], GMn_magenta_err[n_magenta];
	double Q2_orange[n_orange], Q2_orange_err[n_orange], GMn_orange[n_orange], GMn_orange_err[n_orange];
	double Q2_all_blacks[n_all_blacks], Q2_all_blacks_err[n_all_blacks], GMn_all_blacks[n_all_blacks], GMn_all_blacks_err[n_all_blacks];
	double Q2_blue[n_blue], Q2_blue_err[n_blue], GMn_blue[n_blue], GMn_blue_err[n_blue];

	for(int i = 0; i < n_green; i++ ){
		Q2_green[i] = GMn_proposal_green[i][0];
		GMn_green[i] = GMn_proposal_green[i][1];
		Q2_green_err[i] = GMn_proposal_green[i][2];
		GMn_green_err[i] = GMn_proposal_green[i][3];
	}

	for(int i = 0; i < n_magenta; i++ ){
		Q2_magenta[i] = GMn_proposal_magenta[i][0];
		GMn_magenta[i] = GMn_proposal_magenta[i][1];
		Q2_magenta_err[i] = GMn_proposal_magenta[i][2];
		GMn_magenta_err[i] = GMn_proposal_magenta[i][3];
	}

	for(int i = 0; i < n_orange; i++ ){
		Q2_orange[i] = GMn_proposal_orange[i][0];
		GMn_orange[i] = GMn_proposal_orange[i][1];
		Q2_orange_err[i] = GMn_proposal_orange[i][2];
		GMn_orange_err[i] = GMn_proposal_orange[i][3];
	}

	for(int i = 0; i < n_all_blacks; i++ ){
		Q2_all_blacks[i] = GMn_proposal_all_blacks[i][0];
		GMn_all_blacks[i] = GMn_proposal_all_blacks[i][1];
		Q2_all_blacks_err[i] = GMn_proposal_all_blacks[i][2];
		GMn_all_blacks_err[i] = GMn_proposal_all_blacks[i][3];
	}

	for(int i = 0; i < n_blue; i++ ){
		Q2_blue[i] = GMn_proposal_blue[i][0];
		GMn_blue[i] = GMn_proposal_blue[i][1];
		Q2_blue_err[i] = GMn_proposal_blue[i][2];
		GMn_blue_err[i] = GMn_proposal_blue[i][3];
	}

	TGraphErrors *tg_GMn_proposal_green = new TGraphErrors(n_green, Q2_green, GMn_green, Q2_green_err, GMn_green_err);
	tg_GMn_proposal_green->SetMarkerColor(kGreen);
	tg_GMn_proposal_green->SetLineColor(kGreen);
	tg_GMn_proposal_green->SetMarkerStyle(53);
	tg_GMn_proposal_green->SetMarkerSize(1);
	tg_GMn_proposal_green->SetMinimum(0.0);
	tg_GMn_proposal_green->SetTitle("G^{n}_{M} World Data (Normalized to G_{D})");
	tg_GMn_proposal_green->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMn_proposal_green->GetYaxis()->SetTitle("G_{M}^{n}/#mu_{n}G_{D}");

	for(int i = 0; i < n_green; i++ ){
		tg_GMn_proposal_green->SetPoint(i, Q2_green[i], GMn_green[i]);
		tg_GMn_proposal_green->SetPointError(i, Q2_green_err[i], GMn_green_err[i]);
	}

	TGraphErrors *tg_GMn_proposal_magenta = new TGraphErrors(n_magenta, Q2_magenta, GMn_magenta, Q2_magenta_err, GMn_magenta_err);
	tg_GMn_proposal_magenta->SetMarkerColor(kMagenta);
	tg_GMn_proposal_magenta->SetLineColor(kMagenta);
	tg_GMn_proposal_magenta->SetMarkerStyle(52);
	tg_GMn_proposal_magenta->SetMarkerSize(1);
	tg_GMn_proposal_magenta->SetMinimum(0.0);
	tg_GMn_proposal_magenta->SetTitle("Rosenbluth Separation GMn, Normalized G_{M}^{n}");
	tg_GMn_proposal_magenta->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMn_proposal_magenta->GetYaxis()->SetTitle("G_{M}^{n}/#mu_{n}G_{D}");

	for(int i = 0; i < n_magenta; i++ ){
		tg_GMn_proposal_magenta->SetPoint(i, Q2_magenta[i], GMn_magenta[i]);
		tg_GMn_proposal_magenta->SetPointError(i, Q2_magenta_err[i], GMn_magenta_err[i]);
	}

	TGraphErrors *tg_GMn_proposal_orange = new TGraphErrors(n_orange, Q2_orange, GMn_orange, Q2_orange_err, GMn_orange_err);
	tg_GMn_proposal_orange->SetMarkerColor(kOrange);
	tg_GMn_proposal_orange->SetLineColor(kOrange);
	tg_GMn_proposal_orange->SetMarkerStyle(54);
	tg_GMn_proposal_orange->SetMarkerSize(1);
	tg_GMn_proposal_orange->SetMinimum(0.0);
	tg_GMn_proposal_orange->SetTitle("Rosenbluth Separation GMn, Normalized G_{M}^{n}");
	tg_GMn_proposal_orange->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMn_proposal_orange->GetYaxis()->SetTitle("G_{M}^{n}/#mu_{n}G_{D}");

	for(int i = 0; i < n_orange; i++ ){
		tg_GMn_proposal_orange->SetPoint(i, Q2_orange[i], GMn_orange[i]);
		tg_GMn_proposal_orange->SetPointError(i, Q2_orange_err[i], GMn_orange_err[i]);
	}

	TGraphErrors *tg_GMn_proposal_all_blacks = new TGraphErrors(n_all_blacks, Q2_all_blacks, GMn_all_blacks, Q2_all_blacks_err, GMn_all_blacks_err);
	tg_GMn_proposal_all_blacks->SetMarkerColor(kBlack);
	tg_GMn_proposal_all_blacks->SetLineColor(kBlack);
	tg_GMn_proposal_all_blacks->SetMarkerStyle(22);
	tg_GMn_proposal_all_blacks->SetMarkerSize(1);
	tg_GMn_proposal_all_blacks->SetMinimum(0.0);
	tg_GMn_proposal_all_blacks->SetTitle("Rosenbluth Separation GMn, Normalized G_{M}^{n}");
	tg_GMn_proposal_all_blacks->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMn_proposal_all_blacks->GetYaxis()->SetTitle("G_{M}^{n}/#mu_{n}G_{D}");

	for(int i = 0; i < n_all_blacks; i++ ){
		tg_GMn_proposal_all_blacks->SetPoint(i, Q2_all_blacks[i], GMn_all_blacks[i]);
		tg_GMn_proposal_all_blacks->SetPointError(i, Q2_all_blacks_err[i], GMn_all_blacks_err[i]);
	}

	TGraphErrors *tg_GMn_proposal_blue = new TGraphErrors(n_blue, Q2_blue, GMn_blue, Q2_blue_err, GMn_blue_err);
	tg_GMn_proposal_blue->SetMarkerColor(4);
	tg_GMn_proposal_blue->SetLineColor(4);
	tg_GMn_proposal_blue->SetMarkerStyle(52);
	tg_GMn_proposal_blue->SetMarkerSize(1);
	tg_GMn_proposal_blue->SetMinimum(0.0);
	tg_GMn_proposal_blue->SetTitle("Rosenbluth Separation GMn, Normalized G_{M}^{n}");
	tg_GMn_proposal_blue->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMn_proposal_blue->GetYaxis()->SetTitle("G_{M}^{n}/#mu_{n}G_{D}");

	for(int i = 0; i < n_blue; i++ ){
		tg_GMn_proposal_blue->SetPoint(i, Q2_blue[i], GMn_blue[i]);
		tg_GMn_proposal_blue->SetPointError(i, Q2_blue_err[i], GMn_blue_err[i]);
	}

	auto *mg_GMn_proposal = new TMultiGraph();
	mg_GMn_proposal->SetMinimum(0.0);
	mg_GMn_proposal->SetMaximum(1.5);
	mg_GMn_proposal->SetTitle("Rosenbluth Separation GMn, Normalized G_{M}^{n}; Q^{2} (GeV/c)^{2}; G_{M}^{n}/#mu_{n}G_{D}");
	mg_GMn_proposal->Add(tg_GMn_proposal_all_blacks);
	mg_GMn_proposal->Add(tg_GMn_proposal_green);
	mg_GMn_proposal->Add(tg_GMn_proposal_magenta);
	mg_GMn_proposal->Add(tg_GMn_proposal_orange);
	mg_GMn_proposal->Add(tg_GMn_proposal_blue);



	TLegend *tleg_GMn_jboyd;
	if( plot_jboyd ){

		cout << "Plotting JBoyd's GMn points...." << endl;
		double jboyd_SBS8 = jboyd_GMn_SBS8; //0.8783 //0.844, //0.859
		double jboyd_SBS8_err = jboyd_GMn_SBS8_syst_error + jboyd_GMn_SBS8_stat_error; 
		double Q2_SBS8 = jboyd_Q2_SBS8;

		double jboyd_SBS9 = jboyd_GMn_SBS9;
		double jboyd_SBS9_err = jboyd_GMn_SBS9_syst_error + jboyd_GMn_SBS9_stat_error;
		double Q2_SBS9 = jboyd_Q2_SBS9;

		double GMn_jboyd[2] = {jboyd_SBS8, jboyd_SBS9}; //0.86709 //0.844 //0.8508
		double Q2_jboyd[2] = {Q2_SBS8, Q2_SBS9};
		const int n_jboyd = 2;
		double Q2_jboyd_err[2] = {0.0, 0.0};
		double GMn_jboyd_err[2] = {jboyd_SBS8_err, jboyd_SBS9_err};
		TGraphErrors *tg_GMn_jboyd = new TGraphErrors(n_jboyd, Q2_jboyd, GMn_jboyd, Q2_jboyd_err, GMn_jboyd_err);

		tg_GMn_jboyd->SetPoint(0, Q2_jboyd[0], GMn_jboyd[0]);
		tg_GMn_jboyd->SetMarkerStyle(1);
		tg_GMn_jboyd->SetMarkerSize(2);
		tg_GMn_jboyd->SetMarkerColor(kRed);
		tg_GMn_jboyd->SetLineColor(kRed);

///-------------
		double jboyd_SBS8_arr[1] = {jboyd_SBS8}; //0.8783 //0.844, //0.859
		double jboyd_SBS8_err_arr[1] = {jboyd_SBS8_err}; 
		double Q2_SBS8_arr[1] = {Q2_SBS8};
		double Q2_SBS8_err_arr[1] = {0.0};

		TGraphErrors *tg_GMn_jboyd_SBS8 = new TGraphErrors(1, Q2_SBS8_arr, jboyd_SBS8_arr, Q2_SBS8_err_arr, jboyd_SBS8_err_arr );
		tg_GMn_jboyd_SBS8->SetPoint(0, Q2_SBS8_arr[0], jboyd_SBS8_arr[0]);
		tg_GMn_jboyd_SBS8->SetMarkerStyle(23);
		tg_GMn_jboyd_SBS8->SetMarkerSize(1);
		tg_GMn_jboyd_SBS8->SetLineWidth(1);
		tg_GMn_jboyd_SBS8->SetMarkerColor(kRed);
		tg_GMn_jboyd_SBS8->SetLineColor(kRed);

		double jboyd_SBS9_arr[1] = {jboyd_SBS9};
		double jboyd_SBS9_err_arr[1] = {jboyd_SBS9_err};
		double Q2_SBS9_arr[1] = {Q2_SBS9};		
		double Q2_SBS9_err_arr[1] = {0.0};

		TGraphErrors *tg_GMn_jboyd_SBS9 = new TGraphErrors(1, Q2_SBS9_arr, jboyd_SBS9_arr, Q2_SBS9_err_arr, jboyd_SBS9_err_arr );
		tg_GMn_jboyd_SBS9->SetPoint(0, Q2_SBS9_arr[0], jboyd_SBS9_arr[0]);
		tg_GMn_jboyd_SBS9->SetMarkerStyle(22);
		tg_GMn_jboyd_SBS9->SetMarkerSize(1);
		tg_GMn_jboyd_SBS9->SetLineWidth(1);
		tg_GMn_jboyd_SBS9->SetMarkerColor(kCyan);
		tg_GMn_jboyd_SBS9->SetLineColor(kCyan);

		tleg_GMn_jboyd = new TLegend(0.55, 0.60, 0.85, 0.80);
		tleg_GMn_jboyd->AddEntry(tg_GMn_jboyd_SBS8, Form("SBS8 G_{M}^{n} = %.3f +/- %.3f", jboyd_SBS8_arr[0], jboyd_SBS8_err_arr[0] ));
		tleg_GMn_jboyd->AddEntry(tg_GMn_jboyd_SBS9, Form("SBS9 G_{M}^{n} = %.3f +/- %.3f", jboyd_SBS9_arr[0], jboyd_SBS9_err_arr[0] ));

		// plot_GMn_proposal_data();
		// c_GMN_proposal_data->cd()
		// TCanvas *c = new TCanvas("c", "c", 600, 500);
		// tg_GMn_jboyd->Draw("AP+same");


		// mg_GMn_proposal->Add(tg_GMn_jboyd);

		mg_GMn_proposal->Add(tg_GMn_jboyd_SBS8);
		mg_GMn_proposal->Add(tg_GMn_jboyd_SBS9);
	}

	mg_GMn_proposal->Draw("AP");
	if( custom_xmin != -999.9 && custom_xmax != -999.9 ){
		mg_GMn_proposal->GetXaxis()->SetLimits(custom_xmin, custom_xmax);		
	}
	if( custom_ymin != -999.9 && custom_ymax != -999.9 ){
		mg_GMn_proposal->GetYaxis()->SetLimits(custom_ymin, custom_ymax);		
	}
	else{
		mg_GMn_proposal->GetXaxis()->SetLimits(0, 18);		
		mg_GMn_proposal->GetYaxis()->SetRangeUser(0.4, 1.4);		
	}


	if( plot_jboyd ){
		tleg_GMn_jboyd->Draw("same");
	}

	tleg_gmn_proposal_data = new TLegend(0.15, 0.15, 0.40, 0.35);
	tleg_gmn_proposal_data->AddEntry(tg_GMn_proposal_green, "SLAC (Rock)");
	tleg_gmn_proposal_data->AddEntry(tg_GMn_proposal_magenta, "SLAC (Lung)");
	tleg_gmn_proposal_data->AddEntry(tg_GMn_proposal_orange, "Bartel");
	tleg_gmn_proposal_data->AddEntry(tg_GMn_proposal_blue, "CLAS e5");
	tleg_gmn_proposal_data->Draw("same");

	c_GMn_proposal_data->Modified();

}




void plot_gmn_nominal( TCanvas *c_to_plot_on, TLegend *tleg_to_add = NULL, double custom_GMn_const = -999.99 ){

	c_to_plot_on->cd();

	int n_gmn_nominal = 5;
	double GMn_const;

	if( custom_GMn_const != 999.99 ){
		GMn_const = custom_GMn_const;
	}
	else{
		GMn_const = 0.41;
	}

	vector<double> gmn_nominal_Q2= {3.0, 4.5, 7.5, 10, 13.6};

	int x_min = 0.0;
	int x_max = 18.0;
	int increments = 2.0;

	int n_bins = (x_max - x_min)*increments;

	TH1D *h_gmn_nominal = new TH1D("h_gmn_nominal", "h_gmn_nominal", n_bins, x_min, x_max);
	for( size_t i = 0; i < gmn_nominal_Q2.size(); i++ ){
		double Q2_bin = h_gmn_nominal->FindBin( gmn_nominal_Q2[i] );
		h_gmn_nominal->SetBinContent(Q2_bin, GMn_const);
	}
	h_gmn_nominal->SetMarkerColor(kRed);
	h_gmn_nominal->SetMarkerStyle(33);
	h_gmn_nominal->SetMarkerSize(2);
	h_gmn_nominal->Draw("P+same");

	if( tleg_to_add != NULL ){
		TString label_gmn_nominal = "GMn Meas. Pts.";
		if( doesLegendEntryExist( tleg_to_add, label_gmn_nominal) ){

		}
		else{
			tleg_to_add->AddEntry(h_gmn_nominal, label_gmn_nominal.Data());
			c_to_plot_on->Update();
			c_to_plot_on->Modified();			
		}

	}
}

void plot_gmn_nominal_tg( TCanvas *c_to_plot_on, TLegend *tleg_to_add = NULL, 	double GMn_const = 0.01 ){

	c_to_plot_on->cd();

	int n_gmn_nominal = 5;

	double gmn_nominal_Q2[5] = {3.0, 4.5, 7.5, 10, 13.6};
	double gmn_nominal_GMn[5] = {GMn_const, GMn_const, GMn_const, GMn_const, GMn_const};

	TGraph *tg_gmn_nominal = new TGraph(n_gmn_nominal, gmn_nominal_Q2, gmn_nominal_GMn);
	tg_gmn_nominal->SetMarkerColor(kRed);
	tg_gmn_nominal->SetMarkerStyle(48);
	tg_gmn_nominal->SetMarkerSize(2);
	tg_gmn_nominal->Draw("P+same");

	if( tleg_to_add != NULL ){
		tleg_to_add->AddEntry(tg_gmn_nominal, "GMn Meas. Pts.");
		c_to_plot_on->Update();
		c_to_plot_on->Modified();
	}
}


#endif
