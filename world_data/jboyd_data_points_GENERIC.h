#ifndef JBOYD_DATA_POINTS_H
#define JBOYD_DATA_POINTS_H

double SBS8_Q2_theory = 4.5;
double jboyd_Q2_SBS8 = 4.5; 
double jboyd_Q2_SBS8_error = 0.333;

double SBS9_Q2_theory = 4.5;
double jboyd_Q2_SBS9 = 4.5;
double jboyd_Q2_SBS9_error = 0.333;

double SBS4_Q2_theory = 3.0;
double jboyd_Q2_SBS4 = 3.0;
double jboyd_Q2_SBS4_error = 0.333;

vector<double> jboyd_Q2_vec = {jboyd_Q2_SBS8, jboyd_Q2_SBS9, jboyd_Q2_SBS4};
vector<double> jboyd_Q2_error_vec = {jboyd_Q2_SBS8_error, jboyd_Q2_SBS9_error, jboyd_Q2_SBS4_error};

//NO APPROXIMATION

double jboyd_GMn_SBS4 = 0.980; 
double jboyd_GMn_SBS4_stat_error = 0.0300;
double jboyd_GMn_SBS4_syst_error = 0.0150;

//-------------------------------
	// SBS8
//-------------------------------

double jboyd_GMn_SBS8 = 0.92000;
double jboyd_GMn_SBS8_syst_error = 0.0300;
double jboyd_GMn_SBS8_stat_error = 0.0150; 

	//----------------------------------

double SBS8_yield_ratio = 0.3333;
double SBS8_yield_ratio_syst_error = 0.03333;
double SBS8_yield_ratio_stat_error = 0.01500;

double SBS8_scale_factor_ratio = 1.0000; 
double SBS8_scale_factor_ratio_syst_error = 0.0300;
double SBS8_scale_factor_ratio_stat_error = 0.0150; 

//----------------------------------------------------------------
//-------------------------------
	// SBS9
//-------------------------------

double jboyd_GMn_SBS9 = 0.92;
double jboyd_GMn_SBS9_syst_error = 0.0300;
double jboyd_GMn_SBS9_stat_error = 0.0150;

	//----------------------------------

double SBS9_yield_ratio = 0.3333;
double SBS9_yield_ratio_syst_error = 0.0300;
double SBS9_yield_ratio_stat_error = 0.0150;
	//----------------------------------
double SBS9_scale_factor_ratio = 1.000;
double SBS9_scale_factor_ratio_syst_error = 0.0300;
double SBS9_scale_factor_ratio_stat_error = 0.0150;
//-------------------------------
//-------------------------------

//WITH APPROXIATION GEN = 0:
double jboyd_GMn_SBS8_approx = 0.9500;
double jboyd_GMn_SBS8_approx_error = 0.0450;

double jboyd_GMn_SBS9_approx = 0.9500;
double jboyd_GMn_SBS9_approx_error = 0.0450;
//-------------------------------
//-------------------------------

vector<double> jboyd_GMn_vec = {jboyd_GMn_SBS8, jboyd_GMn_SBS9, jboyd_GMn_SBS4};

//Addng the systematic error to the statistical error so that when we plot the statistical error first its error bars hopefully poke out from the back. 
vector<double> jboyd_GMn_stat_error_vec = {jboyd_GMn_SBS8_syst_error + jboyd_GMn_SBS8_stat_error, jboyd_GMn_SBS9_syst_error + jboyd_GMn_SBS9_stat_error, jboyd_GMn_SBS4_syst_error + jboyd_GMn_SBS4_stat_error};
vector<double> jboyd_GMn_syst_error_vec = {jboyd_GMn_SBS8_syst_error, jboyd_GMn_SBS9_syst_error, jboyd_GMn_SBS4_syst_error};

vector<double> jboyd_GMn_approx_vec = {jboyd_GMn_SBS8_approx, jboyd_GMn_SBS9_approx};
vector<double> jboyd_GMn_approx_error_vec = {jboyd_GMn_SBS8_approx_error, jboyd_GMn_SBS9_approx_error};

TGraph *tg_GMn_SBS8_jboyd, *tg_GMn_SBS9_jboyd, *tg_GMn_SBS4_jboyd; 

TGraphErrors *tge_GMn_jboyd;
TGraphErrors *tge_GMn_jboyd_stat;
TGraphErrors *tge_GMn_SBS8_jboyd;
TGraphErrors *tge_GMn_SBS8_jboyd_stat;
TGraphErrors *tge_GMn_SBS8_approx_jboyd;

TGraphErrors *tge_GMn_SBS9_jboyd;
TGraphErrors *tge_GMn_SBS9_jboyd_stat;
TGraphErrors *tge_GMn_SBS9_approx_jboyd;

TGraphErrors *tge_GMn_SBS4_jboyd;
TGraphErrors *tge_GMn_SBS4_jboyd_stat;
TGraphErrors *tge_GMn_SBS4_approx_jboyd;

//----------------------------------------------------------
//----------------------------------------------------------

//NO APPROXIMATION
double jboyd_SBS8_R = 0.322429046; //0.322429046; //0.322429046;
double jboyd_SBS8_R_error = 0.001397691; //0.001397691;// 0.001397691;

double jboyd_SBS8_Rcorr = 0.332163553; // 0.363192502; //0.327055282; //0.327055282;
double jboyd_SBS8_Rcorr_error = 0.015998118; //0.016118127; //0.016118127;

double jboyd_SBS8_ScaleR = 0.0;

//-----------



void plot_jboyd_data(bool bool_plot_jboyd_GMn4 = false){

	cout << "Plotting JBoyd's GMn points...." << endl;
	int linewidth = 2;
// We will plot the tge with statistical error first so that it is behind the systematic.
// We will add the stat error to the syst error so that it hopefully pops out further than the syst. alone

	tge_GMn_jboyd = new TGraphErrors(2, &jboyd_Q2_vec[0], &jboyd_GMn_vec[0], &jboyd_Q2_error_vec[0], &jboyd_GMn_syst_error_vec[0]);

	tg_GMn_SBS8_jboyd = new TGraph(1, &jboyd_Q2_vec[0], &jboyd_GMn_vec[0]);
	tg_GMn_SBS8_jboyd->SetPoint(0, jboyd_Q2_vec[0], jboyd_GMn_vec[0]);
	tg_GMn_SBS8_jboyd->SetMarkerStyle(kFullDotLarge);
	tg_GMn_SBS8_jboyd->SetMarkerColor(kViolet);
	tg_GMn_SBS8_jboyd->SetMarkerSize(1.00000f);

	tg_GMn_SBS9_jboyd = new TGraph(1, &jboyd_Q2_vec[1], &jboyd_GMn_vec[1]);
	tg_GMn_SBS9_jboyd->SetPoint(0, jboyd_Q2_vec[1], jboyd_GMn_vec[1]);
	tg_GMn_SBS9_jboyd->SetMarkerStyle(kFullDotLarge);
	tg_GMn_SBS9_jboyd->SetMarkerColor(4);
	tg_GMn_SBS9_jboyd->SetMarkerSize(1.00000f);

	tg_GMn_SBS4_jboyd = new TGraph(1, &jboyd_Q2_vec[2], &jboyd_GMn_vec[2]);
	tg_GMn_SBS4_jboyd->SetPoint(0, jboyd_Q2_vec[2], jboyd_GMn_vec[2]);
	tg_GMn_SBS4_jboyd->SetMarkerStyle(kFullDotLarge);
	tg_GMn_SBS4_jboyd->SetMarkerColor(8);
	tg_GMn_SBS4_jboyd->SetMarkerSize(1.00000f);

	tge_GMn_jboyd->SetPoint(0, jboyd_Q2_vec[0], jboyd_GMn_vec[0]);
	tge_GMn_jboyd->SetPoint(1, jboyd_Q2_vec[1], jboyd_GMn_vec[1]);
	
	if( bool_plot_jboyd_GMn4 ){
		tge_GMn_jboyd->SetPoint(2, jboyd_Q2_vec[2], jboyd_GMn_vec[2]);
	}

	int markerstyle = 20;
	int markersize = 1;

	tge_GMn_jboyd->SetMarkerStyle(markerstyle);
	tge_GMn_jboyd->SetMarkerSize(3);
	tge_GMn_jboyd->SetMarkerColor(kRed);
	tge_GMn_jboyd->SetLineColor(kRed);

	tge_GMn_SBS8_jboyd = new TGraphErrors(1, &jboyd_Q2_vec[0], &jboyd_GMn_vec[0], &jboyd_Q2_error_vec[0], &jboyd_GMn_syst_error_vec[0]);
	tge_GMn_SBS8_jboyd->SetPoint(0, jboyd_Q2_vec[0], jboyd_GMn_vec[0]);
	tge_GMn_SBS8_jboyd->SetMarkerStyle(markerstyle);
	tge_GMn_SBS8_jboyd->SetMarkerSize(markersize);
	tge_GMn_SBS8_jboyd->SetMarkerColor(kViolet);
	tge_GMn_SBS8_jboyd->SetLineColor(kViolet);
	tge_GMn_SBS8_jboyd->SetLineWidth(linewidth);

	tge_GMn_SBS9_jboyd = new TGraphErrors(1, &jboyd_Q2_vec[1], &jboyd_GMn_vec[1], &jboyd_Q2_error_vec[1], &jboyd_GMn_syst_error_vec[1]);
	tge_GMn_SBS9_jboyd->SetPoint(0, jboyd_Q2_vec[1], jboyd_GMn_vec[1]);
	tge_GMn_SBS9_jboyd->SetMarkerStyle(markerstyle);
	tge_GMn_SBS9_jboyd->SetMarkerSize(markersize);
	tge_GMn_SBS9_jboyd->SetMarkerColor(kBlue);
	tge_GMn_SBS9_jboyd->SetLineColor(kBlue);
	tge_GMn_SBS9_jboyd->SetLineWidth(linewidth);

	tge_GMn_SBS4_jboyd = new TGraphErrors(1, &jboyd_Q2_vec[2], &jboyd_GMn_vec[2], &jboyd_Q2_error_vec[2], &jboyd_GMn_syst_error_vec[2]);
	tge_GMn_SBS4_jboyd->SetPoint(0, jboyd_Q2_vec[2], jboyd_GMn_vec[2]);
	tge_GMn_SBS4_jboyd->SetMarkerStyle(markerstyle);
	tge_GMn_SBS4_jboyd->SetMarkerSize(markersize);
	tge_GMn_SBS4_jboyd->SetMarkerColor(8);
	tge_GMn_SBS4_jboyd->SetLineColor(8);
	tge_GMn_SBS4_jboyd->SetLineWidth(linewidth);

	tge_GMn_SBS8_approx_jboyd = new TGraphErrors(1, &jboyd_Q2_vec[0], &jboyd_GMn_approx_vec[0], &jboyd_Q2_error_vec[0], &jboyd_GMn_approx_error_vec[0]);
	tge_GMn_SBS8_approx_jboyd->SetPoint(0, jboyd_Q2_vec[0], jboyd_GMn_approx_vec[0]);
	tge_GMn_SBS8_approx_jboyd->SetMarkerStyle(markerstyle);
	tge_GMn_SBS8_approx_jboyd->SetMarkerSize(markersize);
	tge_GMn_SBS8_approx_jboyd->SetMarkerColor(kRed);
	tge_GMn_SBS8_approx_jboyd->SetLineColor(kRed);

	tge_GMn_SBS9_approx_jboyd = new TGraphErrors(1, &jboyd_Q2_vec[1], &jboyd_GMn_approx_vec[1], &jboyd_Q2_error_vec[1], &jboyd_GMn_approx_error_vec[1]);
	tge_GMn_SBS9_approx_jboyd->SetPoint(0, jboyd_Q2_vec[1], jboyd_GMn_approx_vec[1]);
	tge_GMn_SBS9_approx_jboyd->SetMarkerStyle(markerstyle);
	tge_GMn_SBS9_approx_jboyd->SetMarkerSize(markersize);
	tge_GMn_SBS9_approx_jboyd->SetMarkerColor(kBlue);
	tge_GMn_SBS9_approx_jboyd->SetLineColor(kBlue);

	/////----------------------------

	tge_GMn_jboyd_stat = new TGraphErrors(2, &jboyd_Q2_vec[0], &jboyd_GMn_vec[0], &jboyd_Q2_error_vec[0], &jboyd_GMn_stat_error_vec[0]);

	tge_GMn_jboyd_stat->SetPoint(0, jboyd_Q2_vec[0], jboyd_GMn_vec[0]);
	tge_GMn_jboyd_stat->SetPoint(1, jboyd_Q2_vec[1], jboyd_GMn_vec[1]);
	
	if( bool_plot_jboyd_GMn4 ){
		tge_GMn_jboyd_stat->SetPoint(2, jboyd_Q2_vec[2], jboyd_GMn_vec[2]);
	}

	tge_GMn_SBS8_jboyd_stat = new TGraphErrors(1, &jboyd_Q2_vec[0], &jboyd_GMn_vec[0], &jboyd_Q2_error_vec[0], &jboyd_GMn_stat_error_vec[0]);
	tge_GMn_SBS8_jboyd_stat->SetPoint(0, jboyd_Q2_vec[0], jboyd_GMn_vec[0]);
	tge_GMn_SBS8_jboyd_stat->SetMarkerStyle(markerstyle);
	tge_GMn_SBS8_jboyd_stat->SetMarkerSize(markersize);
	tge_GMn_SBS8_jboyd_stat->SetMarkerColor(kViolet);
	tge_GMn_SBS8_jboyd_stat->SetLineColor(kRed);
	tge_GMn_SBS8_jboyd_stat->SetLineWidth(linewidth);

	tge_GMn_SBS9_jboyd_stat = new TGraphErrors(1, &jboyd_Q2_vec[1], &jboyd_GMn_vec[1], &jboyd_Q2_error_vec[1], &jboyd_GMn_stat_error_vec[1]);
	tge_GMn_SBS9_jboyd_stat->SetPoint(0, jboyd_Q2_vec[1], jboyd_GMn_vec[1]);
	tge_GMn_SBS9_jboyd_stat->SetMarkerStyle(markerstyle);
	tge_GMn_SBS9_jboyd_stat->SetMarkerSize(markersize);
	tge_GMn_SBS9_jboyd_stat->SetMarkerColor(kBlue);
	tge_GMn_SBS9_jboyd_stat->SetLineColor(kRed);
	tge_GMn_SBS9_jboyd_stat->SetLineWidth(linewidth);

	tge_GMn_SBS4_jboyd_stat = new TGraphErrors(1, &jboyd_Q2_vec[2], &jboyd_GMn_vec[2], &jboyd_Q2_error_vec[2], &jboyd_GMn_stat_error_vec[2]);
	tge_GMn_SBS4_jboyd_stat->SetPoint(0, jboyd_Q2_vec[2], jboyd_GMn_vec[2]);
	tge_GMn_SBS4_jboyd_stat->SetMarkerStyle(markerstyle);
	tge_GMn_SBS4_jboyd_stat->SetMarkerSize(markersize);
	tge_GMn_SBS4_jboyd_stat->SetMarkerColor(8);
	tge_GMn_SBS4_jboyd_stat->SetLineColor(kRed);
	tge_GMn_SBS4_jboyd_stat->SetLineWidth(linewidth);
	
}

void pretty_points( TCanvas *canvas ){

	canvas->cd();

	tge_GMn_SBS4_jboyd_stat->Draw("*same");
	tge_GMn_SBS4_jboyd->Draw("*same");
	tge_GMn_SBS4_jboyd->SetMarkerStyle(kFullDotLarge);
	tge_GMn_SBS4_jboyd->Draw("same");
	tg_GMn_SBS4_jboyd->Draw("*same");
	tg_GMn_SBS4_jboyd->SetLineColor(1);
	tg_GMn_SBS4_jboyd->SetMarkerColor(1);
	tg_GMn_SBS4_jboyd->SetMarkerStyle(kOpenCircle);
	tg_GMn_SBS4_jboyd->Draw("same");

	tge_GMn_SBS8_jboyd_stat->Draw("*same");
	tge_GMn_SBS8_jboyd->Draw("*same");
	tge_GMn_SBS8_jboyd->SetMarkerStyle(kFullDotLarge);
	tge_GMn_SBS8_jboyd->Draw("same");
	tg_GMn_SBS8_jboyd->Draw("*same");
	tg_GMn_SBS8_jboyd->SetLineColor(1);
	tg_GMn_SBS8_jboyd->SetMarkerColor(1);
	tg_GMn_SBS8_jboyd->SetMarkerStyle(kOpenCircle);
	tg_GMn_SBS8_jboyd->Draw("same");

	tge_GMn_SBS9_jboyd_stat->Draw("*same");
	tge_GMn_SBS9_jboyd->Draw("*same");
	tge_GMn_SBS9_jboyd->SetMarkerStyle(kFullDotLarge);
	tge_GMn_SBS9_jboyd->Draw("same");
	tg_GMn_SBS9_jboyd->Draw("*same");
	tg_GMn_SBS9_jboyd->SetLineColor(1);
	tg_GMn_SBS9_jboyd->SetMarkerColor(1);
	tg_GMn_SBS9_jboyd->SetMarkerStyle(kOpenCircle);
	tg_GMn_SBS9_jboyd->Draw("same");
}

#endif