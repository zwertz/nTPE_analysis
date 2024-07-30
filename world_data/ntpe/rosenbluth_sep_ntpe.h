#ifndef ROSENBLUTH_SEP_NTPE_H
#define ROSENBLUTH_SEP_NTPE_H

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

#endif