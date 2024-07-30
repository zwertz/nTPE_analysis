#ifndef CHRISTY_NTPE_GMP_H
#define CHRISTY_NTPE_GMP_H


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

#endif