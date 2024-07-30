#ifndef JBOYD_DATA_POINTS_H
#define JBOYD_DATA_POINTS_H

double SBS8_Q2_theory = 4.5;
double jboyd_Q2_SBS8 = 4.4071; //4.368;
double jboyd_Q2_SBS8_error = 0.310;

double SBS9_Q2_theory = 4.5;
double jboyd_Q2_SBS9 = 4.3927; //4.392;
double jboyd_Q2_SBS9_error = 0.1857;

double SBS4_Q2_theory = 3.0;
double jboyd_Q2_SBS4 = 2.9373;
double jboyd_Q2_SBS4_error = 0.165;

vector<double> jboyd_Q2_vec = {jboyd_Q2_SBS8, jboyd_Q2_SBS9, jboyd_Q2_SBS4};
vector<double> jboyd_Q2_error_vec = {jboyd_Q2_SBS8_error, jboyd_Q2_SBS9_error, jboyd_Q2_SBS4_error};

//NO APPROXIMATION

double jboyd_GMn_SBS4 = 0.978484; //1.00108; //0.9882; //0.963793; //0.9594; //using NonWeighted //weighted: .9473; //0.986481; //0.88548169; //0.870183447; //0.875815108; //0.856942989; //0.870183447; //0.868644662; //0.87007995;
double jboyd_GMn_SBS4_stat_error = 0.013068300; //using NonWeighted //weighted:0.000717; //0.0653; //0.00915483;//0.001182874; //0.010842118; //0.050072968; //0.000447308; //0.041910973; //0.044067905; //0.05056551;
double jboyd_GMn_SBS4_syst_error = 0.018860000; //using NonWeighted //weighted:0.00422; //0.0653; //0.00915483;//0.001182874; //0.010842118; //0.050072968; //0.000447308; //0.041910973; //0.044067905; //0.05056551;

//-------------------------------
	// SBS8
//-------------------------------

double jboyd_GMn_SBS8 = 0.900931; //0.9294; //using NonWeighted //weighted:0.8992; //0.8659; //0.887; //0.908872683; //0.869313056; //0.870183447; //0.875815108; //0.856942989; //0.870183447; //0.868644662; //0.87007995;
double jboyd_GMn_SBS8_syst_error = 0.018900000; //using NonWeighted //weighted:0.0156; //0.00915483; //0.00618213; //0.010842118; //0.050072968; //0.000447308; //0.041910973; //0.044067905; //0.05056551;
double jboyd_GMn_SBS8_stat_error = 0.012095000; //using NonWeighted //weighted:0.0079; //0.064; //0.00915483; //0.00618213; //0.010842118; //0.050072968; //0.000447308; //0.041910973; //0.044067905; //0.05056551;

	//----------------------------------

double SBS8_yield_ratio = 0.3256; //0.3310; //0.3322;
double SBS8_yield_ratio_syst_error = 0.0125; //0.0121;
double SBS8_yield_ratio_stat_error = 0.0040; //0.0339; //0.0038;

double SBS8_scale_factor_ratio = 0.970; //1.0193; //1.0297778;
double SBS8_scale_factor_ratio_syst_error = 0.0384; // 0.037583;
double SBS8_scale_factor_ratio_stat_error = 0.0070; //0.0589; //0.0065;

//----------------------------------------------------------------
//-------------------------------
	// SBS9
//-------------------------------

double jboyd_GMn_SBS9 = 1.00921; //0.9687; // 0.9587; //using NonWeighted //weighted:0.9516; //0.968790674; //0.867767465; //0.88564948; //0.89487207; //0.867767465; //0.88564948; //0.893231519; //0.925263924; //0.91079894;
double jboyd_GMn_SBS9_syst_error = 0.012285000; //2*0.0081; //using NonWeighted //weighted:0.00686; //0.0084 //0.0091544; //0.00855204; //0.011112426; //0.050095563; //0.000574477; //0.042572282; //0.043045523; //0.043729448; //0.05056551;
double jboyd_GMn_SBS9_stat_error = 0.0032800000; // 0.002969; //using NonWeighted //weighted:0.004762; //0.0929; //0.0091544; //0.00855204; //0.011112426; //0.050095563; //0.000574477; //0.042572282; //0.043045523; //0.043729448; //0.05056551;

	//----------------------------------

double SBS9_yield_ratio = 0.3759; //0.4157; //0.3721;
double SBS9_yield_ratio_syst_error = 0.0091; //0.0061;
double SBS9_yield_ratio_stat_error = 0.0010; //0.0060; //0.0008;
	//----------------------------------
double SBS9_scale_factor_ratio = 1.248; //1.2269; //1.1284444;
double SBS9_scale_factor_ratio_syst_error = 0.0267; //0.018281;
double SBS9_scale_factor_ratio_stat_error = 0.0012; //0.0078; //0.0013;
//-------------------------------
//-------------------------------

//WITH APPROXIATION GEN = 0:
double jboyd_GMn_SBS8_approx = 0.885700051; //0.868644662; //0.87007995;
double jboyd_GMn_SBS8_approx_error = 0.060361527; //0.044067905; //0.05056551;

double jboyd_GMn_SBS9_approx = 0.895494512; //0.893231519; //0.925263924; //0.91079894;
double jboyd_GMn_SBS9_approx_error = 0.061255666; //0.043045523; //0.043729448; //0.05056551;
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