#ifndef NTPE_DATA_TABLES_H
#define NTPE_DATA_TABLES_H

//Q2, Q2_err, GMp, GMp_err
vector<vector<double> >christy_ntpe_GMp = {
	{5.9942, 0.00167163, 1.00045, 0.0106811},
	{7.0199, 0.00103557, 0.96732, 0.0150966},
	{7.9432, 0.000709518, 0.943263, 0.0178099},
	{8.994, 0.000498122, 0.934086, 0.0156591},
	{9.8398, 0.000369459, 0.90902, 0.0289806},
	{12.249, 0.000180151, 0.858464, 0.019318},
	{15.721, 8.57e-05, 0.840226, 0.0250849}
};

void plot_christy_ntpe_GMp(){

	//JBoyd's extracted points:
	//Q2 = 4.407051, GMp = 1.032015
	//Q2 = 4.503205, GMp = 1.033596

	TCanvas *c_christy_ntpe_GMp = new TCanvas("c_christy_ntpe_GMp", "c_christy_ntpe_GMp", 600, 500);
	const int n = christy_ntpe_GMp.size();
	double Q2[n], Q2_err[n], GMp[n], GMp_err[n];

	for(int i = 0; i < n; i++ ){
		Q2[i] = christy_ntpe_GMp[i][0];
		Q2_err[i] = christy_ntpe_GMp[i][1];

		GMp[i] = christy_ntpe_GMp[i][2];
		GMp_err[i] = christy_ntpe_GMp[i][3];
	}

	TGraphErrors *tg_christy_ntpe_GMp = new TGraphErrors(n, Q2, GMp, Q2_err, GMp_err);
	gPad->Modified();
	tg_christy_ntpe_GMp->GetXaxis()->SetLimits(0.0, 18.0);
	// tg_christy_ntpe_GMp->SetMinimum(0.0);
	tg_christy_ntpe_GMp->SetTitle("Rosenbluth Separation nTPE, Normalized G_{M}^{p} [Christy, et al.]");
	tg_christy_ntpe_GMp->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_christy_ntpe_GMp->GetYaxis()->SetTitle("G_{M}^{p}/#mu_{p}G_{D}");

	for(int i = 0; i < n; i++ ){
		tg_christy_ntpe_GMp->SetPoint(i, Q2[i], GMp[i]);
		tg_christy_ntpe_GMp->SetPointError(i, Q2_err[i], GMp_err[i]);
	}
	
	// , "Rosenbluth Separation nTPE, Normalized G_{M}^{p}; Q^{2} (GeV/c)^{2}; G_{M}^{p}/#mu_{p}G_{D}

	tg_christy_ntpe_GMp->Draw("AP*");

}

//Q^2, sigma_t x 10^5, sigma_t_err, sigma_L x 10^5, sigma_L_err, G_{M}/(mu_{p} G_{D}) [OPE], Gm/mu_p*G_D error, mu_{p}G_{E}/G_{M}, mu_p* G_E/G_M _error
vector<vector<double>> rosenbluth_sep_ntpe = {
	{5.994,		167,	4,		7.1,	4.6,	1.000,	0.011,	0.75,	0.25},
	{7.020,		104,	3,		9.3,	5.3,	0.967,	0.015,	1.18,	0.35},
	{7.943,		71.0,	2.7,	4.1,	3.9,	0.943,	0.018,	1.0,	0.5},
	{8.993,		49.8,	1.7,	0.7,	3.0, 	0.934,	0.016,	0.5,	1.2},
	{9.840,		36.9,	2.4,	1.9,	3.5,	0.909,	0.029,	1.1,	1.0},
	{12.249,	18.0,	0.8,	1.2,	1.8,	0.858,	0.019,	1.3,	1.1},
	{15.721,	8.6,	0.5,	-0.2,	1.2,	0.840,	0.025,	-0.9,	2.8}
};

void plot_rosenbluth_sep_ntpe(){
	TCanvas *c_rosenbluth_sep_ntpe = new TCanvas("c_rosenbluth_sep_ntpe", "c_rosenbluth_sep_ntpe", 600, 500);
	const int n = rosenbluth_sep_ntpe.size();
	double Q2[n], Q2_err[n], GMp[n], GMp_err[n];

	for(int i = 0; i < n; i++ ){
		Q2[i] = rosenbluth_sep_ntpe[i][0];
		Q2_err[i] = rosenbluth_sep_ntpe[i][1]*1e-5;

		GMp[i] = rosenbluth_sep_ntpe[i][5];
		GMp_err[i] = rosenbluth_sep_ntpe[i][6];
	}

	TGraphErrors *tg_rosenbluth_sep_ntpe = new TGraphErrors(n, Q2, GMp, Q2_err, GMp_err);
	tg_rosenbluth_sep_ntpe->SetMinimum(0.0);
	tg_rosenbluth_sep_ntpe->SetTitle("Rosenbluth Separation nTPE, Normalized G_{M}^{p}");
	tg_rosenbluth_sep_ntpe->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_rosenbluth_sep_ntpe->GetYaxis()->SetTitle("G_{M}^{p}/#mu_{p}G_{D}");

	for(int i = 0; i < n; i++ ){
		tg_rosenbluth_sep_ntpe->SetPoint(i, Q2[i], GMp[i]);
		tg_rosenbluth_sep_ntpe->SetPointError(i, Q2_err[i], GMp_err[i]);
	}
	
	// , "Rosenbluth Separation nTPE, Normalized G_{M}^{p}; Q^{2} (GeV/c)^{2}; G_{M}^{p}/#mu_{p}G_{D}

	tg_rosenbluth_sep_ntpe->Draw("AP*");

}

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

vector<vector<double>> GMn_proposal_green = {
	{2.502673, 0.94289898, 0.0, 0.10295031/2},
	{3.957219251, 0.89370425, 0.0, 0.11304408/2},
	{6.010695187, 0.87525622, 0.0, 0.11708081/2},
	{8.021390374, 0.67847731, 0.0, 0.18773292/2},
	{10.01069519, 0.58828697, 0.0, 0.29673916/2},
};

vector<vector<double>> GMn_proposal_magenta = {
	{1.689839572, 1.05153734, 0.0, 0.10295024/2},
	{2.459893048, 1.01669107, 0.0, 0.08680124/2},
	{3.251336898, 0.96544656, 0.0, 0.12111801/2},
	{4.021390374, 0.92240117, 0.0, 0.16149069/2},
};

vector<vector<double>> GMn_proposal_yellow = {
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

void plot_GMn_proposal_data(bool plot_jboyd = false){

	TCanvas *c_GMN_proposal_data = new TCanvas("c_GMN_proposal_data", "GMN_proposal_data", 600, 500);
	const int n_green = GMn_proposal_green.size();
	const int n_magenta = GMn_proposal_magenta.size();
	const int n_yellow = GMn_proposal_yellow.size();
	const int n_all_blacks = GMn_proposal_all_blacks.size();
	const int n_blue = GMn_proposal_blue.size();


	double Q2_green[n_green], Q2_green_err[n_green], GMn_green[n_green], GMn_green_err[n_green];
	double Q2_magenta[n_magenta], Q2_magenta_err[n_magenta], GMn_magenta[n_magenta], GMn_magenta_err[n_magenta];
	double Q2_yellow[n_yellow], Q2_yellow_err[n_yellow], GMn_yellow[n_yellow], GMn_yellow_err[n_yellow];
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

	for(int i = 0; i < n_yellow; i++ ){
		Q2_yellow[i] = GMn_proposal_yellow[i][0];
		GMn_yellow[i] = GMn_proposal_yellow[i][1];
		Q2_yellow_err[i] = GMn_proposal_yellow[i][2];
		GMn_yellow_err[i] = GMn_proposal_yellow[i][3];
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
	tg_GMn_proposal_green->SetTitle("Rosenbluth Separation GMn, Normalized G_{M}^{n}");
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

	TGraphErrors *tg_GMn_proposal_yellow = new TGraphErrors(n_yellow, Q2_yellow, GMn_yellow, Q2_yellow_err, GMn_yellow_err);
	tg_GMn_proposal_yellow->SetMarkerColor(kYellow);
	tg_GMn_proposal_yellow->SetLineColor(kYellow);
	tg_GMn_proposal_yellow->SetMarkerStyle(54);
	tg_GMn_proposal_yellow->SetMarkerSize(1);
	tg_GMn_proposal_yellow->SetMinimum(0.0);
	tg_GMn_proposal_yellow->SetTitle("Rosenbluth Separation GMn, Normalized G_{M}^{n}");
	tg_GMn_proposal_yellow->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMn_proposal_yellow->GetYaxis()->SetTitle("G_{M}^{n}/#mu_{n}G_{D}");

	for(int i = 0; i < n_yellow; i++ ){
		tg_GMn_proposal_yellow->SetPoint(i, Q2_yellow[i], GMn_yellow[i]);
		tg_GMn_proposal_yellow->SetPointError(i, Q2_yellow_err[i], GMn_yellow_err[i]);
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
	mg_GMn_proposal->Add(tg_GMn_proposal_yellow);
	mg_GMn_proposal->Add(tg_GMn_proposal_blue);

	TLegend *tleg_GMn_jboyd;
	if( plot_jboyd ){

		cout << "Plotting JBoyd's GMn points...." << endl;
		double jboyd_SBS8 = 0.87007995; //0.8783 //0.844, //0.859
		double jboyd_SBS8_err = 0.1056551; 
		double Q2_SBS8 = 4.5038;

		double jboyd_SBS9 = 0.92007;
		double jboyd_SBS9_err = 0.115;
		double Q2_SBS9 = 4.5047;

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
	mg_GMn_proposal->GetXaxis()->SetLimits(0, 11);

	if( plot_jboyd ){
		tleg_GMn_jboyd->Draw("same");
	}

	c_GMN_proposal_data->Modified();

}

// void plot_GMn_jboyd(){

// 	double jboyd_SBS8 = 0.844;
// 	double jboyd_SBS9 = 0.0;

// 	double GMn_jboyd[1] = {0.844};
// 	double Q2_jboyd[1] = {4.5};
// 	const int n_jboyd = 1;
// 	double Q2_jboyd_err[n_jboyd], GMn_jboyd_err[n_jboyd];
// 	TGraphErrors *tg_GMn_jboyd = new TGraphErrors(n_jboyd, Q2_jboyd, GMn_jboyd, Q2_jboyd_err, GMn_jboyd_err);

// 	tg_GMn_jboyd->SetPoint(0, Q2_jboyd[0], GMn_jboyd[0]);
// 	tg_GMn_jboyd->SetMarkerStyle(39);
// 	tg_GMn_jboyd->SetMarkerColor(6);
// 	tg_GMn_jboyd->SetLineColor(6);

// 	plot_GMn_proposal_data();
// 	// c_GMN_proposal_data->cd()
// 	TCanvas *c = new TCanvas("c", "c", 600, 500);
// 	tg_GMn_jboyd->Draw();
// }

vector<vector<double>> GMp_PR12_07_108_hollow_square = {
	{2.8333333333333335, 1.0226252158894646, 0.0, 0.017962005},
	{3.6333333333333333, 1.0419689119170985, 0.0, 0.020478412},
	{5, 1.010880829015544, 0.0, 0.025906735},
	{5, 1.0046632124352333, 0.0, 0.020725385},
	{5, 0.9991364421416236, 0.0, 0.018307425},
	{7.3, 0.9487046632124353, 0.0, 0.019343695},
	{9.633333333333333, 0.8906735751295337, 0.0, 0.019343695},
	{12, 0.8727115716753022, 0.0, 0.019343695},
	{15.733333333333333, 0.8208981001727116, 0.0, 0.025906735},
	{19.466666666666665, 0.7317789291882556, 0.0, 0.027772885},
	{23.233333333333334, 0.7283246977547495, 0.0, 0.033506045},
	{27, 0.7096718480138169, 0.0, 0.04145078},
	{31.2, 0.7207253886010362, 0.0, 0.064248705},
};

vector<vector<double>> GMp_PR12_07_108_solid_circle = {
	{2.5, 1.0578583765112264, 0.0, 0.012435233},
	{3.2333333333333334, 1.053022452504318, 0.0, 0.015198618},
	{4, 1.0398963730569948, 0.0, 0.014853195},
	{5, 1.0281519861830744, 0.0, 0.015198618},
	{6, 1.0018998272884283, 0.0, 0.019343696},
	{7, 0.9728842832469775, 0.0, 0.022107081},
};

vector<vector<double>> GMp_PR12_07_108_hollow_circle ={
	{1.4666666666666666, 0.9715025906735751, 0.0, 0.11433506},
	{1.9666666666666666, 0.9791018998272885, 0.0, 0.074265976},
	{2.5, 1.010880829015544, 0.0, 0.034196891},
	{3.7333333333333334, 0.970811744386874, 0.0, 0.041450777},
};

vector<vector<double>> GMp_PR12_07_108_downward_triangle = {
	{0.13333333333333333, 0.9272884283246978, 0.0, 0.026943005},
	{0.23333333333333334, 0.9417962003454232, 0.0, 0.026943005},
	{0.16666666666666666, 0.9604490500863558, 0.0, 0.026943005},
	{0.26666666666666666, 0.968048359240069, 0.0, 0.026943005},
	{0.5333333333333333, 0.9694300518134715, 0.0, 0.026943005},
	{0.3, 0.9728842832469775, 0.0, 0.026943005},
	{0.5333333333333333, 0.9901554404145079, 0.0, 0.026943005},
	{0.2, 1.0025906735751295, 0.0, 0.024870466},
	{0.7, 1.0184801381692574, 0.0, 0.026943005},
	{0.6666666666666666, 1.0288428324697756, 0.0, 0.026943005},
	{0.7666666666666666, 1.0378238341968913, 0.0, 0.026943005},
	{0.8333333333333334, 1.086873920552677, 0.0, 0.01761658},
	{0.23333333333333334, 0.9355785837651123, 0.0, 0.026943005},
	{0.23333333333333334, 0.9390328151986184, 0.0, 0.026943005},
	{0.26666666666666666, 0.9424870466321243, 0.0, 0.026943005},
	{0.3, 0.9673575129533679, 0.0, 0.026943005},
	{0.5666666666666667, 0.9701208981001728, 0.0, 0.026943005},
	{0.3333333333333333, 0.9728842832469775, 0.0, 0.026943005},
	{0.6333333333333333, 0.990846286701209, 0.0, 0.026943005},
	{0.6666666666666666, 0.9853195164075994, 0.0, 0.026943005},
	{0.7, 1.0018998272884283, 0.0, 0.026943005},
	{0.7333333333333333, 1.017098445595855, 0.0, 0.026943005},
	{1.1666666666666667, 1.0184801381692574, 0.0, 0.026943005},
	{1.4, 1.0329879101899828, 0.0, 0.026943005},
	{0.3333333333333333, 0.9493955094991364, 0.0, 0.026943005},
	{0.5333333333333333, 0.9694300518134715, 0.0, 0.026943005},
	{0.3, 0.9728842832469775, 0.0, 0.026943005},
	{0.7, 0.9860103626943005, 0.0, 0.026943005},
	{0.4, 0.9804835924006909, 0.0, 0.026943005},
	{0.5, 0.9894645941278066, 0.0, 0.026943005},
	{0.9, 1.017098445595855, 0.0, 0.026943005},
	{1.1333333333333333, 1.0191709844559587, 0.0, 0.026943005},
	{0.43333333333333335, 0.9604490500863558, 0.0, 0.026943005},
};

vector<vector<double>> GMp_PR12_07_108_solid_square = {
	{0.7333333333333333, 0.9639032815198618, 0.0, 0.022452504},
	{1.1666666666666667, 1.0177892918825562, 0.0, 0.022452504},
	{1.5333333333333332, 1.030915371329879, 0.0, 0.022452504},
	{1.7666666666666666, 1.0516407599309154, 0.0, 0.022452504},
	{2.966666666666667, 1.0516407599309154, 0.0, 0.022452504},
	{0.9666666666666667, 1.0129533678756477, 0.0, 0.022452504},
	{0.3333333333333333, 0.9507772020725389, 0.0, 0.022452504},
	{0.4, 0.9625215889464594, 0.0, 0.022452504},
	{0.5333333333333333, 0.9887737478411054, 0.0, 0.022452504},
	{0.9666666666666667, 1.013644214162349, 0.0, 0.022452504},
	{1.1333333333333333, 1.017098445595855, 0.0, 0.022452504},
};

void plot_GMp_PR12_07_108(){
	TCanvas *c_GMp_PR12_07_108 = new TCanvas("c_GMp_PR12_07_108", "GMp_PR12_07_108", 600, 500);
	const int n_hollow_square = GMp_PR12_07_108_hollow_square.size();
	const int n_solid_circle = GMp_PR12_07_108_solid_circle.size();
	const int n_hollow_circle = GMp_PR12_07_108_hollow_circle.size();
	const int n_downward_triangle = GMp_PR12_07_108_downward_triangle.size();
	const int n_solid_square = GMp_PR12_07_108_solid_square.size();

	double Q2_hollow_square[n_hollow_square], Q2_hollow_square_err[n_hollow_square], GMn_hollow_square[n_hollow_square], GMn_hollow_square_err[n_hollow_square];
	double Q2_solid_circle[n_solid_circle], Q2_solid_circle_err[n_solid_circle], GMn_solid_circle[n_solid_circle], GMn_solid_circle_err[n_solid_circle];
	double Q2_hollow_circle[n_hollow_circle], Q2_hollow_circle_err[n_hollow_circle], GMn_hollow_circle[n_hollow_circle], GMn_hollow_circle_err[n_hollow_circle];
	double Q2_downward_triangle[n_downward_triangle], Q2_downward_triangle_err[n_downward_triangle], GMn_downward_triangle[n_downward_triangle], GMn_downward_triangle_err[n_downward_triangle];
	double Q2_solid_square[n_solid_square], Q2_solid_square_err[n_solid_square], GMn_solid_square[n_solid_square], GMn_solid_square_err[n_solid_square];

	for(int i = 0; i < n_hollow_square; i++ ){
		Q2_hollow_square[i] = GMp_PR12_07_108_hollow_square[i][0];
		GMn_hollow_square[i] = GMp_PR12_07_108_hollow_square[i][1];
		Q2_hollow_square_err[i] = GMp_PR12_07_108_hollow_square[i][2];
		GMn_hollow_square_err[i] = GMp_PR12_07_108_hollow_square[i][3];
	}

	for(int i = 0; i < n_solid_circle; i++ ){
		Q2_solid_circle[i] = GMp_PR12_07_108_solid_circle[i][0];
		GMn_solid_circle[i] = GMp_PR12_07_108_solid_circle[i][1];
		Q2_solid_circle_err[i] = GMp_PR12_07_108_solid_circle[i][2];
		GMn_solid_circle_err[i] = GMp_PR12_07_108_solid_circle[i][3];
	}

	for(int i = 0; i < n_hollow_circle; i++ ){
		Q2_hollow_circle[i] = GMp_PR12_07_108_hollow_circle[i][0];
		GMn_hollow_circle[i] = GMp_PR12_07_108_hollow_circle[i][1];
		Q2_hollow_circle_err[i] = GMp_PR12_07_108_hollow_circle[i][2];
		GMn_hollow_circle_err[i] = GMp_PR12_07_108_hollow_circle[i][3];
	}

	for(int i = 0; i < n_downward_triangle; i++ ){
		Q2_downward_triangle[i] = GMp_PR12_07_108_downward_triangle[i][0];
		GMn_downward_triangle[i] = GMp_PR12_07_108_downward_triangle[i][1];
		Q2_downward_triangle_err[i] = GMp_PR12_07_108_downward_triangle[i][2];
		GMn_downward_triangle_err[i] = GMp_PR12_07_108_downward_triangle[i][3];
	}

	for(int i = 0; i < n_solid_square; i++ ){
		Q2_solid_square[i] = GMp_PR12_07_108_solid_square[i][0];
		GMn_solid_square[i] = GMp_PR12_07_108_solid_square[i][1];
		Q2_solid_square_err[i] = GMp_PR12_07_108_solid_square[i][2];
		GMn_solid_square_err[i] = GMp_PR12_07_108_solid_square[i][3];
	}

	TGraphErrors *tg_GMp_PR12_07_108_hollow_square = new TGraphErrors(n_hollow_square, Q2_hollow_square, GMn_hollow_square, Q2_hollow_square_err, GMn_hollow_square_err);
	// tg_GMp_PR12_07_108_hollow_square->SetMarkerColor(kGreen);
	// tg_GMp_PR12_07_108_hollow_square->SetLineColor(kGreen);
	tg_GMp_PR12_07_108_hollow_square->SetMarkerStyle(54);
	tg_GMp_PR12_07_108_hollow_square->SetMarkerSize(1);
	tg_GMp_PR12_07_108_hollow_square->SetMinimum(0.0);
	tg_GMp_PR12_07_108_hollow_square->SetTitle("Rosenbluth Separation GMp, Normalized G_{M}^{p}");
	tg_GMp_PR12_07_108_hollow_square->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMp_PR12_07_108_hollow_square->GetYaxis()->SetTitle("G_{M}^{p}/#mu_{p}G_{D}");

	for(int i = 0; i < n_hollow_square; i++ ){
		tg_GMp_PR12_07_108_hollow_square->SetPoint(i, Q2_hollow_square[i], GMn_hollow_square[i]);
		tg_GMp_PR12_07_108_hollow_square->SetPointError(i, Q2_hollow_square_err[i], GMn_hollow_square_err[i]);
	}

	TGraphErrors *tg_GMp_PR12_07_108_solid_circle = new TGraphErrors(n_solid_circle, Q2_solid_circle, GMn_solid_circle, Q2_solid_circle_err, GMn_solid_circle_err);
	// tg_GMp_PR12_07_108_solid_circle->SetMarkerColor(kMagenta);
	// tg_GMp_PR12_07_108_solid_circle->SetLineColor(kMagenta);
	tg_GMp_PR12_07_108_solid_circle->SetMarkerStyle(8);
	tg_GMp_PR12_07_108_solid_circle->SetMarkerSize(1);
	tg_GMp_PR12_07_108_solid_circle->SetMinimum(0.0);
	tg_GMp_PR12_07_108_solid_circle->SetTitle("Rosenbluth Separation GMp, Normalized G_{M}^{p}");
	tg_GMp_PR12_07_108_solid_circle->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMp_PR12_07_108_solid_circle->GetYaxis()->SetTitle("G_{M}^{p}/#mu_{p}G_{D}");

	for(int i = 0; i < n_solid_circle; i++ ){
		tg_GMp_PR12_07_108_solid_circle->SetPoint(i, Q2_solid_circle[i], GMn_solid_circle[i]);
		tg_GMp_PR12_07_108_solid_circle->SetPointError(i, Q2_solid_circle_err[i], GMn_solid_circle_err[i]);
	}

	TGraphErrors *tg_GMp_PR12_07_108_hollow_circle = new TGraphErrors(n_hollow_circle, Q2_hollow_circle, GMn_hollow_circle, Q2_hollow_circle_err, GMn_hollow_circle_err);
	// tg_GMp_PR12_07_108_hollow_circle->SetMarkerColor(kYellow);
	// tg_GMp_PR12_07_108_hollow_circle->SetLineColor(kYellow);
	tg_GMp_PR12_07_108_hollow_circle->SetMarkerStyle(53);
	tg_GMp_PR12_07_108_hollow_circle->SetMarkerSize(1);
	tg_GMp_PR12_07_108_hollow_circle->SetMinimum(0.0);
	tg_GMp_PR12_07_108_hollow_circle->SetTitle("Rosenbluth Separation GMp, Normalized G_{M}^{p}");
	tg_GMp_PR12_07_108_hollow_circle->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMp_PR12_07_108_hollow_circle->GetYaxis()->SetTitle("G_{M}^{p}/#mu_{p}G_{D}");

	for(int i = 0; i < n_hollow_circle; i++ ){
		tg_GMp_PR12_07_108_hollow_circle->SetPoint(i, Q2_hollow_circle[i], GMn_hollow_circle[i]);
		tg_GMp_PR12_07_108_hollow_circle->SetPointError(i, Q2_hollow_circle_err[i], GMn_hollow_circle_err[i]);
	}

	TGraphErrors *tg_GMp_PR12_07_108_downward_triangle = new TGraphErrors(n_downward_triangle, Q2_downward_triangle, GMn_downward_triangle, Q2_downward_triangle_err, GMn_downward_triangle_err);
	// tg_GMp_PR12_07_108_downward_triangle->SetMarkerColor(kBlack);
	// tg_GMp_PR12_07_108_downward_triangle->SetLineColor(kBlack);
	tg_GMp_PR12_07_108_downward_triangle->SetMarkerStyle(23);
	tg_GMp_PR12_07_108_downward_triangle->SetMarkerSize(1);
	tg_GMp_PR12_07_108_downward_triangle->SetMinimum(0.0);
	tg_GMp_PR12_07_108_downward_triangle->SetTitle("Rosenbluth Separation GMp, Normalized G_{M}^{p}");
	tg_GMp_PR12_07_108_downward_triangle->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMp_PR12_07_108_downward_triangle->GetYaxis()->SetTitle("G_{M}^{p}/#mu_{p}G_{D}");

	for(int i = 0; i < n_downward_triangle; i++ ){
		tg_GMp_PR12_07_108_downward_triangle->SetPoint(i, Q2_downward_triangle[i], GMn_downward_triangle[i]);
		tg_GMp_PR12_07_108_downward_triangle->SetPointError(i, Q2_downward_triangle_err[i], GMn_downward_triangle_err[i]);
	}

	TGraphErrors *tg_GMp_PR12_07_108_solid_square = new TGraphErrors(n_solid_square, Q2_solid_square, GMn_solid_square, Q2_solid_square_err, GMn_solid_square_err);
	// tg_GMp_PR12_07_108_solid_square->SetMarkerColor(4);
	// tg_GMp_PR12_07_108_solid_square->SetLineColor(4);
	tg_GMp_PR12_07_108_solid_square->SetMarkerStyle(21);
	tg_GMp_PR12_07_108_solid_square->SetMarkerSize(1);
	tg_GMp_PR12_07_108_solid_square->SetMinimum(0.0);
	tg_GMp_PR12_07_108_solid_square->SetTitle("Rosenbluth Separation GMp, Normalized G_{M}^{p}");
	tg_GMp_PR12_07_108_solid_square->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	tg_GMp_PR12_07_108_solid_square->GetYaxis()->SetTitle("G_{M}^{p}/#mu_{p}G_{D}");

	for(int i = 0; i < n_solid_square; i++ ){
		tg_GMp_PR12_07_108_solid_square->SetPoint(i, Q2_solid_square[i], GMn_solid_square[i]);
		tg_GMp_PR12_07_108_solid_square->SetPointError(i, Q2_solid_square_err[i], GMn_solid_square_err[i]);
	}

	auto *mg_GMp_PR12_07_108 = new TMultiGraph();
	mg_GMp_PR12_07_108->SetMinimum(0.6);
	mg_GMp_PR12_07_108->SetMaximum(1.15);
	mg_GMp_PR12_07_108->SetTitle("Rosenbluth Separation GMp, Normalized G_{M}^{p}; Q^{2} (GeV/c)^{2}; G_{M}^{p}/#mu_{p}G_{D}");
	mg_GMp_PR12_07_108->Add(tg_GMp_PR12_07_108_downward_triangle);
	mg_GMp_PR12_07_108->Add(tg_GMp_PR12_07_108_hollow_square);
	mg_GMp_PR12_07_108->Add(tg_GMp_PR12_07_108_solid_circle);
	mg_GMp_PR12_07_108->Add(tg_GMp_PR12_07_108_hollow_circle);
	mg_GMp_PR12_07_108->Add(tg_GMp_PR12_07_108_solid_square);
	mg_GMp_PR12_07_108->Draw("AP");
	mg_GMp_PR12_07_108->GetXaxis()->SetLimits(-1.5, 32.5);
	c_GMp_PR12_07_108->Modified();

	TLine *tl_GMp_PR12_07_108_horiz_1 = new TLine(-1.5, 1.0, 32.5, 1.0);
	tl_GMp_PR12_07_108_horiz_1->Draw("same");

}
// {epsilon, tau(G_M)^2 + eps(G_E)^2}
vector<vector<double>> nTPE_RS = {
	{0.798389481, 0.00385}, //SBS8
	{0.513874804, 0.00375}, //SBS9

};

void plot_nTPE_RS(){
	int num_points = 2;
	double eps[2] = {0.798389481, 0.513874804};
	double RS[2] = {0.00385, 0.00375};

	TGraph *tg_nTPE_RS = new TGraph(num_points, eps, RS);
	tg_nTPE_RS->SetTitle("nTPE - Reduced C.S., #sigma_{R} vs Virtual Photon Polarization, #epsilon");
	tg_nTPE_RS->GetXaxis()->SetLimits(0.0, 1.0);
	tg_nTPE_RS->GetXaxis()->SetTitle("#epsilon");

	tg_nTPE_RS->GetYaxis()->SetTitle("#tauG_{M}^{2} + #epsilonG_{E}^{2}");
	tg_nTPE_RS->GetYaxis()->SetRangeUser(3.0e-3, 4.3e-3);
	tg_nTPE_RS->GetYaxis()->SetMaxDigits(2);

	TCanvas *c_nTPE_RS = new TCanvas("c_nTPE_RS", "c_nTPE_RS", 600, 500);
	tg_nTPE_RS->Draw("AP*");

	TF1 *tf_nTPE_RS = new TF1("tf_nTPE_RS", "pol1", 0.0, 1.0);
	tf_nTPE_RS->SetLineColor(1);
	tg_nTPE_RS->Fit("tf_nTPE_RS", "R+");

	double y_intercept = tf_nTPE_RS->GetParameter(0);
	double slope = tf_nTPE_RS->GetParameter(1);

	TLine *tl_eps_1 = new TLine(0.513874804, 3e-3, 0.513874804, 4.3e-3);
	tl_eps_1->SetLineStyle(6);
	tl_eps_1->Draw("same");

	TLine *tl_eps_2 = new TLine(0.798389481, 3e-3, 0.798389481, 4.3e-3);
	tl_eps_2->SetLineStyle(6);
	tl_eps_2->Draw("same");

	TText *txt_sbs8 = new TText(0.798389481 - 0.0005, 3.1e-3, "SBS8");
	txt_sbs8->SetTextAngle(90);
	txt_sbs8->SetTextSize(0.035f);
	txt_sbs8->Draw("same");

	TText *txt_sbs9 = new TText(0.513874804 - 0.0005, 3.1e-3, "SBS9");
	txt_sbs9->SetTextAngle(90);
	txt_sbs9->SetTextSize(0.035f);
	txt_sbs9->Draw("same");

	TArrow *ta_tauGM2 = new TArrow(0.025, 3e-3, 0.025, y_intercept, 0.02, "<|>");
	ta_tauGM2->SetFillColor(1);
	ta_tauGM2->Draw();

	TLatex txt_tauGMN2;
	txt_tauGMN2.SetTextSize(0.05f);
	txt_tauGMN2.SetTextAlign(13);
	txt_tauGMN2.DrawLatex(0.03, (y_intercept + 3e-3)/2.0, "#tauG_{M}^{2}");

	TArrow *ta_GE2 = new TArrow(0.52, 0.0037, 0.79, 0.0038, 0.02, "<|>");
	ta_GE2->SetFillColor(1);
	ta_GE2->Draw();

	TLatex txt_GE2;
	txt_GE2.SetTextSize(0.05f);
	txt_GE2.SetTextAlign(13);
	txt_GE2.DrawLatex(0.64, 0.00375, "G_{E}^{2}");

	cout << "----------------------------------------" << endl;
	cout << Form("Intercept ( #tau (G_{M})^{2} ) = %E", tf_nTPE_RS->GetParameter(0)) << endl;
	cout << Form("Slope ( (G_{E})^{2} ) = %E", tf_nTPE_RS->GetParameter(1) ) << endl;
}

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

void plot_nTPE_proposal_fig2(bool plot_jboyd = false, bool overlay_LT_PT = false ){

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

auto mg_nTPE_LT_PT_overlay = new TMultiGraph();
TGraphErrors *tg_nTPE_proposal_fig2_LT;
TGraphErrors *tg_nTPE_proposal_fig2_PT;

void plot_nTPE_proposal_fig2_LT_PT_overlay(bool plot_jboyd = false, bool plot_exp_info = false ){

	TCanvas *c_ntpe_LT_PT_overlay = new TCanvas("c_ntpe_LT_PT_overlay", "c_ntpe_LT_PT_overlay", 700, 500);
	const int n_LT = nTPE_proposal_fig2_green.size();
	const int n_PT = nTPE_proposal_fig2_green_V2.size();

	double eps_LT[n_LT], eps_LT_err[n_LT], ratio_LT[n_LT], ratio_LT_err[n_LT];
	double eps_PT[n_PT], eps_PT_err[n_PT], ratio_PT[n_PT], ratio_PT_err[n_PT];

	for(int i = 0; i < n_LT; i++ ){
		eps_LT[i] = nTPE_proposal_fig2_green[i][0];
		ratio_LT[i] = nTPE_proposal_fig2_green[i][1];
		eps_LT_err[i] = nTPE_proposal_fig2_green[i][2];
		ratio_LT_err[i] = nTPE_proposal_fig2_green[i][3];
	}

	for(int i = 0; i < n_PT; i++ ){
		eps_PT[i] = nTPE_proposal_fig2_green_V2[i][0];
		ratio_PT[i] = nTPE_proposal_fig2_green_V2[i][1];
		eps_PT_err[i] = nTPE_proposal_fig2_green_V2[i][2];
		ratio_PT_err[i] = nTPE_proposal_fig2_green_V2[i][3];
	}

	tg_nTPE_proposal_fig2_LT = new TGraphErrors(n_LT, eps_LT, ratio_LT, eps_LT_err, ratio_LT_err);
	tg_nTPE_proposal_fig2_LT->SetMarkerColor(kBlue);
	tg_nTPE_proposal_fig2_LT->SetLineColor(kBlue);
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
	tg_nTPE_proposal_fig2_PT->SetMarkerColor(kViolet);
	tg_nTPE_proposal_fig2_PT->SetLineColor(kViolet);
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
	mg_nTPE_LT_PT_overlay->SetTitle("LT & PT Data - Neutron FF Ratio vs Momentum Transfer, Q^{2}; Q^{2} (GeV/c)^{2}; #propto_{n} G_{E}^{n}/G_{M}^{n}");
	mg_nTPE_LT_PT_overlay->Add(tg_nTPE_proposal_fig2_LT);
	mg_nTPE_LT_PT_overlay->Add(tg_nTPE_proposal_fig2_PT);

	mg_nTPE_LT_PT_overlay->Draw("APL");

	TLegend *tl_nTPE_LT_PT_overlay = new TLegend(0.12, 0.70, 0.52, 0.90);
	tl_nTPE_LT_PT_overlay->AddEntry(tg_nTPE_proposal_fig2_LT, "LT Separation Data");
	tl_nTPE_LT_PT_overlay->AddEntry(tg_nTPE_proposal_fig2_PT, "Polarization Transfer Data");
	tl_nTPE_LT_PT_overlay->Draw("same");

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
		txt_overlay_LT.SetTextColor(kBlue);
		txt_overlay_LT.SetTextSize(0.035f);
		txt_overlay_LT.DrawLatex(4.54, 0.35, Form("LT Separation --> %0.3f", tf_nTPE_LT_simple_linear->Eval(4.5)));
		txt_overlay_LT.Draw("same");

		TLatex txt_overlay_PT;
		txt_overlay_PT.SetTextColor(kViolet);
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
	tf_nTPE_proposal_fig2_LT->SetLineColor(kBlue);

	TF1 *tf_nTPE_proposal_fig2_PT = new TF1("tf_nTPE_proposal_fig2_PT", "pol2", 0, 6.5);
	tg_nTPE_proposal_fig2_PT->Fit("tf_nTPE_proposal_fig2_PT", "R+0");
	tf_nTPE_proposal_fig2_PT->SetLineStyle(3);
	tf_nTPE_proposal_fig2_PT->SetLineColor(kViolet);	

	mg_nTPE_LT_PT_overlay->Draw("AP");
	tf_nTPE_proposal_fig2_LT->Draw("same");
	tf_nTPE_proposal_fig2_PT->Draw("same");

	tl_nTPE_LT_PT_overlay->Draw("same");

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

TF1 *ye_arrington_GMn_fit;
TGraph *tg_ye_arrington;

vector<double> ye_GMn_vec;
vector<double> ye_Q2_vec;

void plot_ye_arrington( int FF = 4, bool plot_jboyd = false ){

	TCanvas *c_ye_arrington = new TCanvas("c_ye_arrington", "c_ye_arrington", 600, 500);

	double exp_all_Q2[33] = {
		2.5, 4.0, 6.0, 8.0, 10.0,
		1.75, 2.5, 3.25, 4.00,
		0.200, 0.235, 0.503, 0.652, 0.784,
		0.07100, 0.125, 0.359, 0.894,
		0.10, 0.193, 0.30, 0.40, 0.50, 0.60,
		1.13637, 1.59, 2.04547, 2.5, 2.95557, 3.40910, 3.86363, 4.31820, 4.69695
	};

	double exp_all_Q2_err[33];

	double exp_all_Gmn[33] = {
		0.94308, 0.90717, 0.90717, 0.90717, 0.61150,
		1.05200, 1.01400, 0.96700, 0.92300,
		0.93100, 0.96200, 1.03200, 1.03700,1.04300,
		0.99000, 0.967, 0.989, 1.062,
		0.94810, 0.95110, 0.95770, 0.9694, 0.9689, 0.99390, 
		1.10294, 1.01627, 1.00364, 1.01620, 0.9940, 0.98781, 0.99029, 0.92990, 1.10180
	};

	double exp_all_GMn_err[33] = {
		0.06714, 0.06605, 0.06685, 0.10785, 0.15681,
		0.052, 0.0440, 0.061, 0.081,
		0.05725, 0.011, 0.019, 0.016, 0.016,
		0.016, 0.016, 0.019, 0.024,
		0.02663, 0.02664, 0.05536, 0.03058, 0.02700, 0.02864,
		0.0231, 0.01971, 0.02069, 0.02268, 0.02219, 0.02165, 0.02222, 0.02442, 0.07999
	};

	TGraphErrors *tg_GMn_exp_all = new TGraphErrors(33, exp_all_Q2, exp_all_Gmn, exp_all_Q2_err, exp_all_GMn_err);
	tg_GMn_exp_all->SetLineColor(kBlue);
	tg_GMn_exp_all->SetMarkerColor(kBlue);
	tg_GMn_exp_all->SetMarkerStyle(8);

	auto *mg_ye_arrington = new TMultiGraph();
	mg_ye_arrington->Add(tg_GMn_exp_all);

	TLegend *tleg_GMn_jboyd;
	if( plot_jboyd ){

		cout << "Plotting JBoyd's GMn points...." << endl;
		double jboyd_SBS8 = 0.87007995; //0.8783 //0.844, //0.859
		double jboyd_SBS8_err = 0.07056551; 
		double Q2_SBS8 = 4.5038;

		double jboyd_SBS9 = 0.92007;
		double jboyd_SBS9_err = 0.075;
		double Q2_SBS9 = 4.5047;

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

		mg_ye_arrington->Add(tg_GMn_jboyd_SBS8);
		mg_ye_arrington->Add(tg_GMn_jboyd_SBS9);
	}

	c_ye_arrington->SetLogx();
	mg_ye_arrington->GetXaxis()->SetLimits(0.01, 11);
	mg_ye_arrington->GetYaxis()->SetLimits(0.2, 1.2);

//Parameterization
	double GMn_a0 = 0.257758326959;
	double GMn_a1 = -1.079540642058;
	double GMn_a2 = 1.182183812195;
	double GMn_a3 = 0.711015085833;
	double GMn_a4 = -1.348080936796;
	double GMn_a5 = -1.662444025208;
	double GMn_a6 = 2.624354426029;
	double GMn_a7 = 1.751234494568;
	double GMn_a8 = -4.922300878888;
	double GMn_a9 = 3.197892727312;
	double GMn_a10 = -0.712072389946;
	double GMn_a11 = 0.0;
	double GMn_a12 = 0.0;

	vector<double> GMn_aN_vec = {GMn_a0, GMn_a1, GMn_a2, GMn_a3, GMn_a4, GMn_a5, GMn_a6, GMn_a7, GMn_a8, GMn_a9, GMn_a10};

	double tcut = 4*pow(0.13957, 2);
	double t0 = -0.7;

	// double z_SBS8 = (sqrt(tcut + pow(Q2_SBS8_arr[0], 2)) - sqrt(tcut - t0) )/(sqrt(tcut + pow(Q2_SBS8_arr[0], 2)) + sqrt(tcut - t0) );


	ye_arrington_GMn_fit = new TF1("ye_arrington_fit", "pol10", 0, 11);
	ye_arrington_GMn_fit->SetParameter(0, GMn_a0);
	ye_arrington_GMn_fit->SetParameter(1, GMn_a1);
	ye_arrington_GMn_fit->SetParameter(2, GMn_a2);
	ye_arrington_GMn_fit->SetParameter(3, GMn_a3);
	ye_arrington_GMn_fit->SetParameter(4, GMn_a4);
	ye_arrington_GMn_fit->SetParameter(5, GMn_a5);
	ye_arrington_GMn_fit->SetParameter(6, GMn_a6);
	ye_arrington_GMn_fit->SetParameter(7, GMn_a7);
	ye_arrington_GMn_fit->SetParameter(8, GMn_a8);
	ye_arrington_GMn_fit->SetParameter(9, GMn_a9);
	ye_arrington_GMn_fit->SetParameter(10, GMn_a10);

	double GMn_fit_max = 11;
	double GMn_fit_step_size = 0.05;
	int GMn_fit_counts = 0;

	double GMn_calc_val = 0;
	double Q2_calc_val = 0;
	double delta = 0.71;

	for( double q2 = 0.01; q2 < GMn_fit_max; q2 += GMn_fit_step_size ){

		double gmn_calc = 0.0;
		double z = 0.0;
		double G_D = 0.0;
		G_D = 1/(pow( ( 1 + (pow(q2, 2)/delta)), 2));
			
		z = (sqrt(tcut + pow(q2, 2)) - sqrt(tcut - t0) )/(sqrt(tcut + pow(q2, 2)) + sqrt(tcut - t0) );

		for( size_t k = 0; k <= 10; k++ ){
			gmn_calc += GMn_aN_vec[k]*pow( z, k );
		}

		ye_GMn_vec.push_back(gmn_calc);
		ye_Q2_vec.push_back(q2);

	}

	tg_ye_arrington = new TGraph(ye_GMn_vec.size(), &ye_GMn_vec[0], &ye_Q2_vec[0]);
	mg_ye_arrington->Add(tg_ye_arrington);

	mg_ye_arrington->Draw("AP");

}
#endif

