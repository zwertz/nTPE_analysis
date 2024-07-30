#ifndef NTPE_RS_H
#define NTPE_RS_H

// {epsilon, tau(G_M)^2 + eps(G_E)^2}
vector<vector<double>> nTPE_RS = {
	{0.798389481, 0.00385}, //SBS8
	{0.513874804, 0.00375}, //SBS9

};

void plot_nTPE_RS(){
	int num_points = 2;
	double eps[2] = {0.798389481, 0.513874804};
	// double RS[2] = {0.00385, 0.00375};
	double RS[2] = {0.001613, 0.001558};

	double SBS8_eps = 0.798389481;
	double SBS9_eps = 0.513874804;

	TGraph *tg_nTPE_RS = new TGraph(num_points, eps, RS);
	tg_nTPE_RS->SetTitle("Neutron Reduced C.S. From Parameterization, #sigma_{R} vs Virtual Photon Polarization, #epsilon");
	tg_nTPE_RS->GetXaxis()->SetLimits(0.0, 1.0);
	tg_nTPE_RS->GetXaxis()->SetTitle("#epsilon");

	tg_nTPE_RS->GetYaxis()->SetTitle("#tauG_{M}^{2} + #epsilonG_{E}^{2}");
	double y_min = 1.25e-3;
	double y_max = 1.7e-3;
	tg_nTPE_RS->GetYaxis()->SetRangeUser(y_min, y_max);
	tg_nTPE_RS->GetYaxis()->SetMaxDigits(2);

	TCanvas *c_nTPE_RS = new TCanvas("c_nTPE_RS", "c_nTPE_RS", 600, 500);
	tg_nTPE_RS->Draw("AP*");

	TF1 *tf_nTPE_RS = new TF1("tf_nTPE_RS", "pol1", 0.0, 1.0);
	tf_nTPE_RS->SetLineColor(1);
	tg_nTPE_RS->Fit("tf_nTPE_RS", "R+");

	double y_intercept = tf_nTPE_RS->GetParameter(0);
	double slope = tf_nTPE_RS->GetParameter(1);

	TLine *tl_eps_1 = new TLine(0.513874804, y_min, 0.513874804, y_max);
	tl_eps_1->SetLineStyle(6);
	tl_eps_1->Draw("same");

	TLine *tl_eps_2 = new TLine(0.798389481, y_min, 0.798389481, y_max);
	tl_eps_2->SetLineStyle(6);
	tl_eps_2->Draw("same");

	TText *txt_sbs8 = new TText(0.798389481 - 0.0005, 1.03*y_min, "SBS8");
	txt_sbs8->SetTextAngle(90);
	txt_sbs8->SetTextSize(0.035f);
	txt_sbs8->Draw("same");

	TText *txt_sbs9 = new TText(0.513874804 - 0.0005, 1.03*y_min, "SBS9");
	txt_sbs9->SetTextAngle(90);
	txt_sbs9->SetTextSize(0.035f);
	txt_sbs9->Draw("same");

	TArrow *ta_tauGM2 = new TArrow(0.025, y_min, 0.025, y_intercept, 0.02, "<|>");
	ta_tauGM2->SetFillColor(1);
	ta_tauGM2->Draw();

	TLatex txt_tauGMN2;
	txt_tauGMN2.SetTextSize(0.05f);
	txt_tauGMN2.SetTextAlign(13);
	txt_tauGMN2.DrawLatex(0.03, (y_intercept + y_min)/2.0, "#tauG_{M}^{2}");

	TArrow *ta_GE2 = new TArrow(0.82, RS[1], 0.82, RS[0], 0.01, "<|>");
	ta_GE2->SetFillColor(1);
	ta_GE2->Draw();

	double RS_diff = fabs(RS[0]-RS[1]);

	Double_t angle_rad = atan(1000*RS_diff/(SBS8_eps - SBS9_eps));
    if (angle_rad < 0) {
        angle_rad += TMath::Pi();
    }
    Double_t angle_deg = 2*angle_rad*(180/TMath::Pi());

	TEllipse *GEn_arc = new TEllipse( SBS9_eps, RS[1], 0.15, 0.00007, 0, angle_deg );
	GEn_arc->SetFillStyle(3003);
	GEn_arc->SetFillColorAlpha(16, 0.5);
	// GEn_arc->SetLineColor(16);
	// GEn_arc->SetLineWidth(1);
	GEn_arc->Draw("only");

	TLine *tl_G2_horiz = new TLine(0.513874804, RS[1], 0.798389481, RS[1]);
	tl_G2_horiz->SetLineStyle(3);
	tl_G2_horiz->SetLineColor(14);
	tl_G2_horiz->Draw("same");

	TLatex txt_GE2;
	txt_GE2.SetTextSize(0.05f);
	txt_GE2.SetTextAlign(13);
	txt_GE2.DrawLatex(0.835, RS[0], "G_{E}^{2}");

	cout << "----------------------------------------" << endl;
	cout << Form("Intercept ( #tau (G_{M})^{2} ) = %E", tf_nTPE_RS->GetParameter(0)) << endl;
	cout << Form("Slope ( (G_{E})^{2} ) = %E", tf_nTPE_RS->GetParameter(1) ) << endl;
}

#endif