#ifndef CHRISTY_FFR_PARAM_H
#define CHRISTY_FFR_PARAM_H

double par_christy_RS[2] = {-0.46, 0.12};

Double_t christy_RS(Double_t *x, Double_t *par_christy_RS ){
	double Mp = 0.938272;

	double tau = x[0]/(4*Mp*Mp);

	double RS = 1 + (-0.46)*tau + (0.12)*tau*tau;

	return sqrt(RS);
}

Double_t calc_christy_RS( double q2 ){
	double Mp = 0.938272;

	double tau = q2/(4*Mp*Mp);

	double RS = 1 + (-0.46)*tau + (0.12)*tau*tau;

	return sqrt(RS);
}

Double_t calc_christy_RS_error( double q2 ){

	double Mp = 0.938272;
	double tau = q2/(4*Mp*Mp);

	double val1 = (-0.46)*tau;
	double c1_error = 0.12;

	double val2 = (0.12)*tau*tau;
	double c2_error = 0.10;

	double error = CalculateErrorAdditionSubtraction(val1, c1_error, val2, c2_error);

	return error;
}

Double_t christy_RS_error_high_fit( Double_t *x, Double_t *par_christy_RS ){
	return calc_christy_RS( x[0] ) + calc_christy_RS_error( x[0] );
}

Double_t christy_RS_error_low_fit( Double_t *x, Double_t *par_christy_RS ){
	return calc_christy_RS( x[0] ) - calc_christy_RS_error( x[0] );
}

TF1 *tf1_christy_RS, *tf1_christy_RS_error_high, *tf1_christy_RS_error_low;

void plot_christy_ffr_param( bool plot_jboyd = false){

	double christy_RS_x_min = 0;
	double christy_RS_x_max = 18;
	double christy_RS_y_min = -1.0;
	double christy_RS_y_max = 2.5;
	double christy_RS_fit_Npx = 5000;

	TCanvas *c_christy_RS = new TCanvas("c_christy_RS","c_christy_RS", 600, 500);
	
	TH1D *h_christy_RS = new TH1D("h_christy_RS", "Christy RS Parameterization", christy_RS_fit_Npx, christy_RS_x_min, christy_RS_x_max );
	h_christy_RS->SetLineColor(0);
	h_christy_RS->SetStats(0);
	h_christy_RS->GetYaxis()->SetRangeUser(christy_RS_y_min, christy_RS_y_max);
	h_christy_RS->GetYaxis()->SetTitle("#sqrt{RS}( = (#mu_{p}G_{E}/G_{M})^{2} in OPE)");
	h_christy_RS->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
	h_christy_RS->Draw();

	tf1_christy_RS = new TF1("tf1_christy_RS", christy_RS, christy_RS_x_min, christy_RS_x_max, 0);
	tf1_christy_RS->SetNpx(christy_RS_fit_Npx);
	tf1_christy_RS->SetLineColor(kBlack);
	tf1_christy_RS->SetLineWidth(1);

	tf1_christy_RS_error_high = new TF1("tf1_christy_RS_error_high", christy_RS_error_high_fit, christy_RS_x_min, christy_RS_x_max, 0 );

	// tf1_christy_RS->Draw();

}

#endif