//cutvar.C
//Author: Ezekiel Wertz
//Companion implementation.A class to hold important quanties related to cut stability and systematic studies. Cutvar should only handle functions/processes on a case-by-case basis for data, mc p, mc n, and not all at once. A companion class to this is stability_analysis.h and stability_analysis.C. The goal of the cutvar class is to organize delta x and the given cut variable of interest and create histograms within the slices of the cut variable of interest. Ideally it will create a delta x distribution for the given range of the cut variable and the function will return a vector of these histograms. This classes supports mode o which is narrow slices and mode 1 which handles very wide slices. This array of histograms is then sent to the stability_analysis class which handles fitting the delta x distributions and determining Rsf values

#include "../include/cutvar.h"

//public constructor implementation. Constructor, will initialize a cutvar object. Which will be centrally used to handle a lot of the necessary info for the stability study. All private variables will be based on the cutvar input through conditionals.

cutvar::cutvar(TString myVar,TString datCut,TString daFlag,int daMode, TString Left_Right, TString daKine,int daField,double dx_low, double dx_high, double dx_bin, double W2_low, double W2_high, double W2_bin,TChain* C ){

	kine = daKine;
	sbs_field = daField;

	//Store the cut variable independantly
	CutVar = myVar;
	CutString = datCut;
	DatavMC_flag = daFlag;
	
	//Initialize other private variables
	dx_hist_low = dx_low;
	dx_hist_high = dx_high;
	dx_hist_bin = dx_bin;
	W2_hist_low = W2_low;
	W2_hist_high = W2_high;
	W2_hist_bin = W2_bin;

	//Initialize the cutvar slicing mode value
	cutvar_mode = daMode;
	left_right = Left_Right;

	//Call function to initialize all of the potential xMinxMax ranges
	cutvar::initialize_xMin_xMax();

	//For each parameter we want to evaluate, make a conditional. Then initialize all private variables. This will need to be expanded as cut variables are investigated.
	if(myVar == "dy"){
	dxHistogramName = "hcal_dx__hcal_dy";
	W2HistogramName = "W2__hcal_dy";
	AxisTitle = "hcal_dy";
	xMin_xMax_range = xMin_xMax_dy; 
	cut_hist_low = -2;
	cut_hist_high = 2;
	cut_hist_bin =  100;
	}else if(myVar == "BBtr_vz"){
        dxHistogramName = "hcal_dx__tr_vz";
        W2HistogramName = "W2__tr_vz";
	AxisTitle = "track vz";
        xMin_xMax_range = xMin_xMax_vz;
        cut_hist_low = -0.1;
        cut_hist_high = 0.1;
        cut_hist_bin = 200;
	}else if(myVar == "BBps_e"){
        dxHistogramName = "hcal_dx__ps_e";
        W2HistogramName = "W2__ps_e";
	AxisTitle = "preshower energy";
        xMin_xMax_range = xMin_xMax_ps;
        cut_hist_low = 0;
        cut_hist_high = 2.6;
        cut_hist_bin = 250;
	}else if(myVar == "BBgem_nhits"){
	dxHistogramName = "hcal_dx__BBgem_nhits";
	W2HistogramName = "W2__BBgem_nhits";
	AxisTitle ="BBgem_nhits";
	xMin_xMax_range = xMin_xMax_BBgem_nhits;
	cut_hist_low = 0.0;
	cut_hist_high = 6.0;
	cut_hist_bin = 6;
	}else if(myVar == "BBgem_chi2ndf"){
	dxHistogramName = "hcal_dx__BBgem_chi2ndf";
	W2HistogramName = "W2__BBgem_chi2ndf";
	AxisTitle = "BBgem_chi2ndf";
	xMin_xMax_range = xMin_xMax_BBgem_chi2ndf;
	cut_hist_low = 0;
	cut_hist_high = 100;
	cut_hist_bin = 250;
	}else if(myVar == "ehcal"){
        dxHistogramName = "hcal_dx__hcal_e";
        W2HistogramName = "W2__hcal_e";
	AxisTitle = "hcal energy";
        xMin_xMax_range = xMin_xMax_hcal_e;
        cut_hist_low = 0;
        cut_hist_high = 1;
        cut_hist_bin = 200;
	}else if(myVar == "W2"){
        dxHistogramName = "hcal_dx__W2";
        W2HistogramName = "W2__W2";
	AxisTitle = "W^{2}";
        xMin_xMax_range = xMin_xMax_W2;
        cut_hist_low = W2_low;
        cut_hist_high = W2_high;
        cut_hist_bin = W2_bin;
	}else if(myVar == "xexp"){
        dxHistogramName = "hcal_dx__hcal_xexp";
        W2HistogramName = "W2__hcal_xexp";
	AxisTitle = "x expect";
        xMin_xMax_range = xMin_xMax_x_exp;
        cut_hist_low = -3;
        cut_hist_high = 3;
        cut_hist_bin = 600;
	}else if(myVar == "yexp"){
        dxHistogramName = "hcal_dx__yexp";
        W2HistogramName = "W2__yexp";
	AxisTitle = "y expect";
        xMin_xMax_range = xMin_xMax_y_exp;
        cut_hist_low = -2;
        cut_hist_high = 2;
        cut_hist_bin = 400;
	}else if(myVar == "BB_E_over_p"){
	dxHistogramName = "hcal_dx__e_over_p";
	W2HistogramName = "W2__e_over_p";
	AxisTitle = "e_over_p";
	xMin_xMax_range = xMin_xMax_e_over_p;
	cut_hist_low = 0.0;
	cut_hist_high = 1.5;
	cut_hist_bin = 200;
	}else if(myVar == "hcal_sh_atime_diff"){
	dxHistogramName = "hcal_dx__hcal_sh_atime_diff";
	W2HistogramName = "W2__hcal_sh_atime_diff";
	AxisTitle = "hcal_sh_atime_diff";
	xMin_xMax_range = xMin_xMax_coin_diff;
	cut_hist_low = -20;
	cut_hist_high = 20;
	cut_hist_bin = 400;
	}else if(myVar == "BBtr_r_x-BBtr_r_th*0.9"){
	dxHistogramName = "hcal_dx__BB_Optics_X";
        W2HistogramName = "W2__BB_Optics_X";
        AxisTitle = "BBtr_r_x-BBtr_r_th*0.9";
        xMin_xMax_range = xMin_xMax_Optics_x;
        cut_hist_low = -0.5;
        cut_hist_high = 0.5;
        cut_hist_bin = 200;
	}else if(myVar == "BBtr_r_y-0.9*BBtr_r_ph"){
	dxHistogramName = "hcal_dx__BB_Optics_Y";
        W2HistogramName = "W2__BB_Optics_Y";
        AxisTitle = "BBtr_r_y-0.9*BBtr_r_ph";
        xMin_xMax_range = xMin_xMax_Optics_y;
        cut_hist_low = -0.13;
        cut_hist_high = 0.13;
        cut_hist_bin = 200;
	}else{
	cout << "Error in cut variable string " << myVar << "." << endl;
	return;
	}

	//Now create the cut histograms we care about
	dx_hist = cutvar::make2DdxCutHisto(C);

	W2_hist = cutvar::make2DW2CutHisto(C);


}//end constructor

//Copy constructor that takes a reference as input
cutvar::cutvar(cutvar &myCut){

	CutVar = myCut.getCutVar();
        CutString = myCut.getCutString();
        DatavMC_flag = myCut.getDataMCFlag();
        dx_hist_low = myCut.get_dx_hist_low();
        dx_hist_high = myCut.get_dx_hist_high();
        dx_hist_bin = myCut.get_dx_hist_bin();
        W2_hist_low = myCut.get_W2_hist_low();
        W2_hist_high = myCut.get_W2_hist_high();
        W2_hist_bin = myCut.get_W2_hist_bin();
	dxHistogramName = myCut.getdxHistoName();
        W2HistogramName = myCut.getW2HistoName();
        AxisTitle = myCut.getAxisTitle();
        xMin_xMax_range = myCut.getXMinXMaxRange();
        cut_hist_low = myCut.get_cut_hist_low();
        cut_hist_high = myCut.get_cut_hist_high();
        cut_hist_bin =  myCut.get_cut_hist_bin();
	dx_hist = myCut.get2DdxCutHisto();
	W2_hist = myCut.get2DW2CutHisto();
	kine = myCut.getKine();
	sbs_field = myCut.get_sbs_field();
	cutvar_mode = myCut.get_cutvar_mode();

}

//Destructor
//Might have dynamically allocated memory or pointers. May have to worry about that.
cutvar::~cutvar(){
}

//Function that will create a histogram of dx vs the cut variable. Relies on initialization of histo info from constructor.
TH2D* cutvar::make2DdxCutHisto(TChain* C){
//Draw the histogram with the given cut conditions
	if(DatavMC_flag == "data"){
	C->Draw(Form("dx:%s>>%s(%f,%f,%f,%f,%f,%f)",CutVar.Data(),dxHistogramName.Data(),cut_hist_bin,cut_hist_low,cut_hist_high,dx_hist_bin,dx_hist_low,dx_hist_high),CutString.Data(),"COLZ");
	}else if(DatavMC_flag == "mc_p" || DatavMC_flag == "mc_n"){
	C->Draw(Form("dx:%s>>%s(%f,%f,%f,%f,%f,%f)",CutVar.Data(),dxHistogramName.Data(),cut_hist_bin,cut_hist_low,cut_hist_high,dx_hist_bin,dx_hist_low,dx_hist_high),Form("Corr_MC_weight*(%s)",CutString.Data()),"COLZ");
	}else {
	cout << "Error: Given a data mc type flag that I can't handle: " << DatavMC_flag << endl;
	}

//Retrieve and customize the histogram
TH2D* mydxHist = (TH2D*) gDirectory->Get(dxHistogramName.Data());
mydxHist->SetXTitle(AxisTitle.Data());
mydxHist->SetYTitle("hcal_dx");

//Modify the histogram names for data or mc depending on the input parameter.
TString newdxHistoName = DatavMC_flag +"_"+ dxHistogramName;

mydxHist->SetName(newdxHistoName.Data());
return mydxHist;
}//end make dx Histo

TH2D* cutvar::make2DW2CutHisto(TChain* C){
//Draw the histogram with the given cut conditions
        if(DatavMC_flag == "data"){
        C->Draw(Form("W2:%s>>%s(%f,%f,%f,%f,%f,%f)",CutVar.Data(),W2HistogramName.Data(),cut_hist_bin,cut_hist_low,cut_hist_high,400.0*5.0,W2_hist_low,5.0),CutString.Data(),"COLZ");
        }else if(DatavMC_flag == "mc_p" || DatavMC_flag == "mc_n"){
        C->Draw(Form("W2:%s>>%s(%f,%f,%f,%f,%f,%f)",CutVar.Data(),W2HistogramName.Data(),cut_hist_bin,cut_hist_low,cut_hist_high,400.0*5.0,W2_hist_low,5.0),Form("Corr_MC_weight*(%s)",CutString.Data()),"COLZ");
        }else {
        cout << "Error: Given a data mc type flag that I can't handle: " << DatavMC_flag << endl;
        }

//Retrieve and customize the histogram
TH2D* myW2Hist = (TH2D*) gDirectory->Get(W2HistogramName.Data());
myW2Hist->SetXTitle(AxisTitle.Data());
myW2Hist->SetYTitle("W2");

//Modify the histogram names for data or mc depending on the input parameter.
TString newW2HistoName = DatavMC_flag +"_"+ W2HistogramName;

myW2Hist->SetName(newW2HistoName.Data());
return myW2Hist;
}//end make W2 Histo


//Function that will create a histogram of W2 vs the cut variable. Relies on initialization of histo info from constructor.
//TH2D* cutvar::make2DW2CutHisto(TChain* C);

TString cutvar::getdxHistoName(){ return dxHistogramName; }

TString cutvar::getW2HistoName(){ return W2HistogramName; }

TString cutvar::getAxisTitle(){return AxisTitle; }

TString cutvar::getCutVar(){return CutVar; }

TString cutvar::getDataMCFlag(){return DatavMC_flag;}

TString cutvar::getCutString(){return CutString;}

TString cutvar::getKine(){return kine;}

double cutvar::get_dx_hist_low(){return dx_hist_low;}

double cutvar::get_dx_hist_high(){return dx_hist_high;}

double cutvar::get_dx_hist_bin(){return dx_hist_bin;}

double cutvar::get_W2_hist_low(){return W2_hist_low;}

double cutvar::get_W2_hist_high(){return W2_hist_high;}

double cutvar::get_W2_hist_bin(){return W2_hist_bin;}

double cutvar::get_cut_hist_low(){return cut_hist_low;}

double cutvar::get_cut_hist_high(){return cut_hist_high;}

double cutvar::get_cut_hist_bin(){return cut_hist_bin;}

int cutvar::get_cutvar_mode(){return cutvar_mode;}

int cutvar::get_sbs_field(){return sbs_field;}

vector<pair<double,double>> cutvar::getXMinXMaxRange(){return xMin_xMax_range; }

TH2D* cutvar::get2DdxCutHisto(){return dx_hist;}

TH2D* cutvar::get2DW2CutHisto(){return W2_hist;}

//A function that projects the vertical axis of the TH2D into a given xMin, xMax range. In our case the vertical axis is delta x and the horizontal axis is a given cut variable. So this function is what actually handles the generation of the vector of 1-dimensional delta x distributions for a given cut variable value range.
vector<TH1D*> cutvar::sliceAndProjectHisto_xMinxMax(TH2D* histo2D,TString xAxisName, TString yAxisName ){
	//create the vector to hold all of our slice histograms
	vector<TH1D*> myHistVec;
	//loop over the slice vector associated with the class variable
	for(int i=0; i<xMin_xMax_range.size(); i++){
	//Get xMin and the associated bin
	double xMin = xMin_xMax_range[i].first;
	int binMin = histo2D->GetXaxis()->FindBin(xMin);

	//Get xMax and the associated bin
	double xMax = xMin_xMax_range[i].second;
	int binMax = histo2D->GetXaxis()->FindBin(xMax);

	//Make a specific name for the histogram
	TString histName = yAxisName +"__"+ xAxisName +"_"+ utility::doubToTString(xMin,4) +"_to_"+ utility::doubToTString(xMax,4) +"_"+DatavMC_flag; 
	//Make the 1 dimensional histogram from a projection of the 2D histogram
	TH1D *projY = histo2D->ProjectionY(histName.Data(),binMin,binMax);
	myHistVec.push_back(projY);
	}
	//After the loop return the vector of histos
	return myHistVec;	
}//end slice and project

//Function that handles the initialization of the xMinXMax vectors dependent on the slice mode, the kinematic, and the magnetic field setting
void cutvar::initialize_xMin_xMax(){

  //First conditional for cutvar slice mode
  if(cutvar_mode == 0){
  //choosing mode 0 to handle the narrow slices, basically like differentials
  //primarily used for determing/optimizing cut selections

	//Second set of conditional for kinematic and magnetic field setting
	//Necessary to implement: SBS4 50%, SBS8 50%, 70%, 100%, SBS9 70%
	//All other kinematics and magnetic field settings not supported at this time
	if(kine == "SBS4" && sbs_field == 30){
	
	//BB preshower energy range
        xMin_xMax_ps = {{0.10,0.15},{0.15,0.20},{0.20,0.25},{0.25,0.30},{0.30,0.35},{0.35,0.4},{0.40,0.45},{0.45,0.50},{0.50,0.55},{0.55,0.60},{0.60,0.65},{0.65,0.70},{0.70,0.75},{0.75,0.80},{0.80,0.85},{0.85,0.90},{0.90,0.95},{0.95,1.00},{1.00,1.05},{1.05,1.10},{1.10,1.15},{1.15,1.20},{1.20,1.25},{1.25,1.30},{1.30,1.35},{1.35,1.40},{1.40,1.45},{1.45,1.50}};
	//number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
	//BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,2.5},{2.5,5},{5,7.5},{7.5,10},{10,12.5},{12.5,15},{15,17.5},{17.5,20},{20,22.5},{22.5,25},{25,27.5},{27.5,30},{30,32.5},{32.5,35},{35,37.5},{37.5,40},{40,42.5},{42.5,45}};
	//target vertex range
        xMin_xMax_vz = {{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09}};
	//W2 range
        xMin_xMax_W2 = {{0.35,0.4},{0.4,0.45},{0.45,0.5},{0.5,0.55},{0.55,0.6},{0.6,0.65},{0.65,0.7},{0.7,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3},{1.3,1.35},{1.35,1.4},{1.4,1.45},{1.45,1.5}};
	//E over p range
        xMin_xMax_e_over_p = {{0.65,0.70},{0.70,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2}};
	//HCal Energy range
        xMin_xMax_hcal_e = {{0.01,0.015},{0.015,0.02},{0.02,0.025},{0.025,0.03},{0.03,0.035},{0.035,0.04},{0.04,0.045},{0.045,0.05},{0.05,0.055},{0.055,0.06},{0.06,0.065},{0.065,0.07},{0.07,0.075},{0.075,0.08},{0.08,0.085},{0.085,0.09},{0.09,0.095},{0.095,0.1}};
	//xMin_xMax_hcal_e = {{0.01,0.0125},{0.0125,0.015},{0.015,0.0175},{0.0175,0.02},{0.02,0.0225},{0.0225,0.025},{0.025,0.0275},{0.0275,0.03},{0.03,0.0325},{0.0325,0.035},{0.035,0.0375},{0.0375,0.04},{0.04,0.0425},{0.0425,0.045},{0.045,0.0475},{0.0475,0.05},{0.05,0.0525},{0.0525,0.055},{0.055,0.0575},{0.0575,0.06},{0.06,0.0625},{0.0625,0.065},{0.065,0.0675},{0.0675,0.07},{0.07,0.0725},{0.0725,0.075},{0.075,0.0775},{0.0775,0.08},{0.08,0.0825},{0.0825,0.085},{0.085,0.0875},{0.0875,0.09},{0.09,0.0925},{0.0925,0.095},{0.095,0.0975},{0.0975,0.1}};
	//xMin_xMax_hcal_e = {{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1},{0.1,0.11},{0.11,0.12},{0.12,0.13},{0.13,0.14},{0.14,0.15},{0.15,0.16},{0.16,0.17},{0.17,0.18},{0.18,0.19},{0.19,0.2},{0.2,0.23},{0.23,0.26},{0.26,0.29},{0.29,0.32},{0.32,0.35},{0.35,0.38},{0.38,0.41}};
        //xMin_xMax_hcal_e = {{0.010,0.015},{0.015,0.020},{0.020,0.025},{0.025,0.030},{0.030,0.035},{0.035,0.040},{0.040,0.045},{0.045,0.050},{0.050,0.055},{0.055,0.060}};
	//HCal - Shower analog time diff range
        xMin_xMax_coin_diff = {{-7.0,-6.0},{-6.0,-5.0},{-5.0,-4.0},{-4.0,-3.0},{-3.0,-2.0},{-2.0,-1.0},{-1.0,0.0},{0.0,1.0},{1.0,2.0},{2.0,3.0}};
	//dy range
	xMin_xMax_dy = {{-0.9,-0.8},{-0.8,-0.7},{-0.7,-0.6},{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8},{0.8,0.9}};
	//Fid cut for x expect and y expect
	//Fid X, x expect range
	xMin_xMax_x_exp = {{-1.4,-1.3},{-1.3,-1.2},{-1.2,-1.1},{-1.1,-1.0},{-1.0,-0.9},{-0.9,-0.8},{-0.8,-0.7},{-0.7,-0.6},{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8}};
	//Fid Y, y expect range
	xMin_xMax_y_exp = {{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6}};

	xMin_xMax_Optics_x = {{-0.15,-0.1},{-0.1,-0.05},{-0.05,0.0},{0.0,0.05},{0.05,0.1},{0.1,0.15},{0.15,0.2},{0.2,0.25},{0.25,0.3},{0.3,0.35},{0.35,0.4}};
	xMin_xMax_Optics_y ={{-0.10,-0.09},{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1}};
	
	}else if(kine == "SBS4" && sbs_field == 50){
	//BB preshower energy range
        xMin_xMax_ps = {{0.10,0.15},{0.15,0.20},{0.20,0.25},{0.25,0.30},{0.30,0.35},{0.35,0.4},{0.40,0.45},{0.45,0.50},{0.50,0.55},{0.55,0.60},{0.60,0.65},{0.65,0.70},{0.70,0.75},{0.75,0.80},{0.80,0.85},{0.85,0.90},{0.90,0.95},{0.95,1.00},{1.00,1.05},{1.05,1.10},{1.10,1.15},{1.15,1.20},{1.20,1.25},{1.25,1.30},{1.30,1.35},{1.35,1.40},{1.40,1.45},{1.45,1.50}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,5},{5,10},{10,15},{15,20},{20,25},{25,30},{30,35},{35,40},{40,45}};
        //target vertex range
        xMin_xMax_vz = {{-0.09,-0.07},{-0.07,-0.05},{-0.05,-0.03},{-0.03,-0.01},{-0.01,0.01},{0.01,0.03},{0.03,0.05},{0.05,0.07},{0.07,0.09}};
        //W2 range
        xMin_xMax_W2 = {{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8},{0.8,0.9},{0.9,1.0},{1.0,1.1},{1.1,1.2},{1.2,1.3},{1.3,1.4},{1.4,1.5}};
        //E over p range
        xMin_xMax_e_over_p = {{0.55,0.65},{0.65,0.75},{0.75,0.85},{0.85,0.95},{0.95,1.05},{1.05,1.15},{1.15,1.25}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.01,0.03},{0.03,0.06},{0.06,0.09},{0.09,0.12},{0.12,0.15},{0.15,0.18},{0.18,0.21},{0.21,0.24},{0.24,0.27},{0.27,0.30},{0.30,0.33},{0.33,0.36},{0.36,0.39}};
        //xMin_xMax_hcal_e = {{0.010,0.015},{0.015,0.020},{0.020,0.025},{0.025,0.030},{0.030,0.035},{0.035,0.040},{0.040,0.045},{0.045,0.050},{0.050,0.055},{0.055,0.060}};
        //HCal - Shower analog time diff range
        xMin_xMax_coin_diff = {{-8.0,-6.0},{-6.0,-4.0},{-4.0,-2.0},{-2.0,0.0},{0.0,2.0},{2.0,4.0},{4.0,6.0}};
        //dy range
        xMin_xMax_dy = {{-1.0,-0.8},{-0.8,-0.6},{-0.6,-0.4},{-0.4,-0.2},{-0.2,0.0},{0.0,0.2},{0.2,0.4},{0.4,0.6},{0.6,0.8},{0.8,1.0}};
        //Fid cut for x expect and y expect
        //Fid X, x expect range
        xMin_xMax_x_exp = {{-1.6,-1.4},{-1.4,-1.2},{-1.2,-1.0},{-1.0,-0.8},{-0.8,-0.6},{-0.6,-0.4},{-0.4,-0.2},{-0.2,0.0},{0.0,0.2},{0.2,0.4},{0.4,0.6},{0.6,0.8},{0.8,1.0},{1.0,1.2}};
        //Fid Y, y expect range
        xMin_xMax_y_exp = {{-0.6,-0.4},{-0.4,-0.2},{-0.2,0.0},{0.0,0.2},{0.2,0.4},{0.4,0.6},{0.6,0.8}};
	}else if(kine == "SBS8" && sbs_field == 50){
	xMin_xMax_ps = {{0.10,0.15},{0.15,0.20},{0.20,0.25},{0.25,0.30},{0.30,0.35},{0.35,0.4},{0.40,0.45},{0.45,0.50},{0.50,0.55},{0.55,0.60},{0.60,0.65},{0.65,0.70},{0.70,0.75},{0.75,0.80},{0.80,0.85},{0.85,0.90},{0.90,0.95},{0.95,1.00},{1.00,1.05},{1.05,1.10},{1.10,1.15},{1.15,1.20},{1.20,1.25},{1.25,1.30},{1.30,1.35},{1.35,1.40},{1.40,1.45},{1.45,1.50},{1.5,1.55},{1.55,1.6},{1.6,1.65},{1.65,1.7},{1.7,1.75},{1.75,1.8},{1.8,1.85},{1.85,1.9},{1.9,1.95},{1.95,2.0},{2.0,2.05},{2.05,2.1},{2.1,2.15},{2.15,2.2}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,2.5},{2.5,5},{5,7.5},{7.5,10},{10,12.5},{12.5,15},{15,17.5},{17.5,20},{20,22.5},{22.5,25},{25,27.5},{27.5,30}};
        //target vertex range
        xMin_xMax_vz = {{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09}};
        //W2 range
        xMin_xMax_W2 = {{0.45,0.5},{0.5,0.55},{0.55,0.6},{0.6,0.65},{0.65,0.7},{0.7,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3},{1.3,1.35}};
        //E over p range
        xMin_xMax_e_over_p = {{0.65,0.70},{0.70,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.01,0.015},{0.015,0.02},{0.02,0.025},{0.025,0.03},{0.03,0.035},{0.035,0.04},{0.04,0.045},{0.045,0.05},{0.05,0.055},{0.055,0.06},{0.06,0.065},{0.065,0.07},{0.07,0.075},{0.075,0.08},{0.08,0.085},{0.085,0.09},{0.09,0.095},{0.095,0.1}};
        //xMin_xMax_hcal_e = {{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1},{0.1,0.11},{0.11,0.12},{0.12,0.13},{0.13,0.14},{0.14,0.15},{0.15,0.16},{0.16,0.17},{0.17,0.18},{0.18,0.19},{0.19,0.2},{0.2,0.23},{0.23,0.26},{0.26,0.29},{0.29,0.32},{0.32,0.35},{0.35,0.38},{0.38,0.41}};
        //HCal - Shower analog time diff range
        xMin_xMax_coin_diff = {{-8.0,-7.0},{-7.0,-6.0},{-6.0,-5.0},{-5.0,-4.0},{-4.0,-3.0},{-3.0,-2.0},{-2.0,-1.0},{-1.0,0.0},{0.0,1.0},{1.0,2.0},{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //dy range
        xMin_xMax_dy = {{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6}};
        //Fid cut for x expect and y expect
        //Fid X, x expect range
        xMin_xMax_x_exp = {{-1.8,-1.7},{-1.7,-1.6},{-1.6,-1.5},{-1.5,-1.4},{-1.4,-1.3},{-1.3,-1.2},{-1.2,-1.1},{-1.1,-1.0},{-1.0,-0.9},{-0.9,-0.8},{-0.8,-0.7},{-0.7,-0.6},{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8},{0.8,0.9},{0.9,1.0},{1.0,1.1}};
	//Fid Y, y expect range
        xMin_xMax_y_exp = {{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6}};
        xMin_xMax_Optics_x = {{-0.25,-0.2},{-0.2,-0.15},{-0.15,-0.1},{-0.1,-0.05},{-0.05,0.0},{0.0,0.05},{0.05,0.1},{0.1,0.15},{0.15,0.2},{0.2,0.25},{0.25,0.3},{0.3,0.35}};
        xMin_xMax_Optics_y ={{-0.10,-0.09},{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1}};
	}else if(kine == "SBS8" && sbs_field == 70){
	 //BB preshower energy range
        xMin_xMax_ps = {{0.10,0.15},{0.15,0.20},{0.20,0.25},{0.25,0.30},{0.30,0.35},{0.35,0.4},{0.40,0.45},{0.45,0.50},{0.50,0.55},{0.55,0.60},{0.60,0.65},{0.65,0.70},{0.70,0.75},{0.75,0.80},{0.80,0.85},{0.85,0.90},{0.90,0.95},{0.95,1.00},{1.00,1.05},{1.05,1.10},{1.10,1.15},{1.15,1.20},{1.20,1.25},{1.25,1.30},{1.30,1.35},{1.35,1.40},{1.40,1.45},{1.45,1.50},{1.5,1.55},{1.55,1.6},{1.6,1.65},{1.65,1.7},{1.7,1.75},{1.75,1.8},{1.8,1.85},{1.85,1.9},{1.9,1.95},{1.95,2.0},{2.0,2.05},{2.05,2.1},{2.1,2.15},{2.15,2.2},{2.2,2.25},{2.25,2.3}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,2.5},{2.5,5},{5,7.5},{7.5,10},{10,12.5},{12.5,15},{15,17.5},{17.5,20},{20,22.5},{22.5,25},{25,27.5},{27.5,30}};
        //target vertex range
        xMin_xMax_vz = {{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09}};
        //W2 range
        xMin_xMax_W2 = {{0.45,0.5},{0.5,0.55},{0.55,0.6},{0.6,0.65},{0.65,0.7},{0.7,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3},{1.3,1.35}};
        //E over p range
        xMin_xMax_e_over_p = {{0.65,0.70},{0.70,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.01,0.015},{0.015,0.02},{0.02,0.025},{0.025,0.03},{0.03,0.035},{0.035,0.04},{0.04,0.045},{0.045,0.05},{0.05,0.055},{0.055,0.06},{0.06,0.065},{0.065,0.07},{0.07,0.075},{0.075,0.08},{0.08,0.085},{0.085,0.09},{0.09,0.095},{0.095,0.1}};
        //xMin_xMax_hcal_e = {{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1},{0.1,0.11},{0.11,0.12},{0.12,0.13},{0.13,0.14},{0.14,0.15},{0.15,0.16},{0.16,0.17},{0.17,0.18},{0.18,0.19},{0.19,0.2},{0.2,0.23},{0.23,0.26},{0.26,0.29},{0.29,0.32},{0.32,0.35},{0.35,0.38},{0.38,0.41}};
	//HCal - Shower analog time diff range
        xMin_xMax_coin_diff = {{-8.0,-7.0},{-7.0,-6.0},{-6.0,-5.0},{-5.0,-4.0},{-4.0,-3.0},{-3.0,-2.0},{-2.0,-1.0},{-1.0,0.0},{0.0,1.0},{1.0,2.0},{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //dy range
        xMin_xMax_dy = {{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6}};
        //Fid cut for x expect and y expect
        //Fid X, x expect range
        xMin_xMax_x_exp = {{-1.5,-1.4},{-1.4,-1.3},{-1.3,-1.2},{-1.2,-1.1},{-1.1,-1.0},{-1.0,-0.9},{-0.9,-0.8},{-0.8,-0.7},{-0.7,-0.6},{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8}};
        //Fid Y, y expect range
        xMin_xMax_y_exp = {{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6}};
	xMin_xMax_Optics_x = {{-0.25,-0.2},{-0.2,-0.15},{-0.15,-0.1},{-0.1,-0.05},{-0.05,0.0},{0.0,0.05},{0.05,0.1},{0.1,0.15},{0.15,0.2},{0.2,0.25},{0.25,0.3}};
        xMin_xMax_Optics_y ={{-0.10,-0.09},{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1}};

	}else if(kine == "SBS8" && sbs_field == 100){
	xMin_xMax_ps = {{0.10,0.15},{0.15,0.20},{0.20,0.25},{0.25,0.30},{0.30,0.35},{0.35,0.4},{0.40,0.45},{0.45,0.50},{0.50,0.55},{0.55,0.60},{0.60,0.65},{0.65,0.70},{0.70,0.75},{0.75,0.80},{0.80,0.85},{0.85,0.90},{0.90,0.95},{0.95,1.00},{1.00,1.05},{1.05,1.10},{1.10,1.15},{1.15,1.20},{1.20,1.25},{1.25,1.30},{1.30,1.35},{1.35,1.40},{1.40,1.45},{1.45,1.50},{1.5,1.55},{1.55,1.6},{1.6,1.65},{1.65,1.7},{1.7,1.75},{1.75,1.8},{1.8,1.85},{1.85,1.9},{1.9,1.95},{1.95,2.0},{2.0,2.05},{2.05,2.1},{2.1,2.15},{2.15,2.2}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,2.5},{2.5,5},{5,7.5},{7.5,10},{10,12.5},{12.5,15},{15,17.5},{17.5,20},{20,22.5},{22.5,25},{25,27.5},{27.5,30}};
        //target vertex range
        xMin_xMax_vz = {{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09}};
        //W2 range
        xMin_xMax_W2 = {{0.45,0.5},{0.5,0.55},{0.55,0.6},{0.6,0.65},{0.65,0.7},{0.7,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3},{1.3,1.35}};
        //E over p range
        xMin_xMax_e_over_p = {{0.65,0.70},{0.70,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.01,0.015},{0.015,0.02},{0.02,0.025},{0.025,0.03},{0.03,0.035},{0.035,0.04},{0.04,0.045},{0.045,0.05},{0.05,0.055},{0.055,0.06},{0.06,0.065},{0.065,0.07},{0.07,0.075},{0.075,0.08},{0.08,0.085},{0.085,0.09},{0.09,0.095},{0.095,0.1}};
        //xMin_xMax_hcal_e = {{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1},{0.1,0.11},{0.11,0.12},{0.12,0.13},{0.13,0.14},{0.14,0.15},{0.15,0.16},{0.16,0.17},{0.17,0.18},{0.18,0.19},{0.19,0.2},{0.2,0.23},{0.23,0.26},{0.26,0.29},{0.29,0.32},{0.32,0.35},{0.35,0.38},{0.38,0.41}};
        //HCal - Shower analog time diff range
        xMin_xMax_coin_diff = {{-8.0,-7.0},{-7.0,-6.0},{-6.0,-5.0},{-5.0,-4.0},{-4.0,-3.0},{-3.0,-2.0},{-2.0,-1.0},{-1.0,0.0},{0.0,1.0},{1.0,2.0},{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //dy range
        xMin_xMax_dy = {{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6}};
        //Fid cut for x expect and y expect
        //Fid X, x expect range
        xMin_xMax_x_exp = {{-1.4,-1.3},{-1.3,-1.2},{-1.2,-1.1},{-1.1,-1.0},{-1.0,-0.9},{-0.9,-0.8},{-0.8,-0.7},{-0.7,-0.6},{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8},{0.8,0.9},{0.9,1.0},{1.0,1.1},{1.1,1.2},{1.2,1.3}};
        //Fid Y, y expect range
        xMin_xMax_y_exp = {{-0.8,-0.7},{-0.7,-0.6},{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7}};
	xMin_xMax_Optics_x = {{-0.3,-0.25},{-0.25,-0.2},{-0.2,-0.15},{-0.15,-0.1},{-0.1,-0.05},{-0.05,0.0},{0.0,0.05},{0.05,0.1},{0.1,0.15},{0.15,0.2},{0.2,0.25},{0.25,0.3}};
        xMin_xMax_Optics_y ={{-0.10,-0.09},{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1}};
	}else if(kine == "SBS9" && sbs_field == 70){
	//BB preshower energy range
        xMin_xMax_ps = {{0.15,0.20},{0.20,0.25},{0.25,0.30},{0.30,0.35},{0.35,0.4},{0.40,0.45},{0.45,0.50},{0.50,0.55},{0.55,0.60},{0.60,0.65},{0.65,0.70},{0.70,0.75},{0.75,0.80},{0.80,0.85},{0.85,0.90},{0.90,0.95},{0.95,1.00},{1.00,1.05},{1.05,1.10},{1.10,1.15},{1.15,1.20},{1.20,1.25},{1.25,1.30},{1.30,1.35},{1.35,1.40}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{2.0,3.0},{3.0,4.0},{4.0,5.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,2.5},{2.5,5},{5,7.5},{7.5,10},{10,12.5},{12.5,15},{15,17.5},{17.5,20},{20,22.5},{22.5,25},{25,27.5},{27.5,30}};
        //target vertex range
        xMin_xMax_vz = {{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09}};
        //W2 range
        xMin_xMax_W2 = {{0.45,0.5},{0.5,0.55},{0.55,0.6},{0.6,0.65},{0.65,0.7},{0.7,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3},{1.3,1.35}};
        //E over p range
        xMin_xMax_e_over_p = {{0.65,0.70},{0.70,0.75},{0.75,0.8},{0.8,0.85},{0.85,0.9},{0.9,0.95},{0.95,1.0},{1.0,1.05},{1.05,1.1},{1.1,1.15},{1.15,1.2},{1.2,1.25},{1.25,1.3}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.01,0.015},{0.015,0.02},{0.02,0.025},{0.025,0.03},{0.03,0.035},{0.035,0.04},{0.04,0.045},{0.045,0.05},{0.05,0.055},{0.055,0.06},{0.06,0.065},{0.065,0.07},{0.07,0.075},{0.075,0.08},{0.08,0.085},{0.085,0.09},{0.09,0.095},{0.095,0.1}};
        //xMin_xMax_hcal_e = {{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1},{0.1,0.11},{0.11,0.12},{0.12,0.13},{0.13,0.14},{0.14,0.15},{0.15,0.16},{0.16,0.17},{0.17,0.18},{0.18,0.19},{0.19,0.2},{0.2,0.23},{0.23,0.26},{0.26,0.29},{0.29,0.32},{0.32,0.35},{0.35,0.38},{0.38,0.41}};
        //HCal - Shower analog time diff range
        xMin_xMax_coin_diff = {{-8.0,-7.0},{-7.0,-6.0},{-6.0,-5.0},{-5.0,-4.0},{-4.0,-3.0},{-3.0,-2.0},{-2.0,-1.0},{-1.0,0.0},{0.0,1.0},{1.0,2.0},{2.0,3.0},{3.0,4.0},{4.0,5.0}};
        //dy range
        xMin_xMax_dy = {{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5}};
        //Fid cut for x expect and y expect
        //Fid X, x expect range
        xMin_xMax_x_exp = {{-1.2,-1.1},{-1.1,-1.0},{-1.0,-0.9},{-0.9,-0.8},{-0.8,-0.7},{-0.7,-0.6},{-0.6,-0.5},{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4},{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8},{0.8,0.9},{0.9,1.0},{1.0,1.1}};
        //Fid Y, y expect range
        xMin_xMax_y_exp = {{-0.5,-0.4},{-0.4,-0.3},{-0.3,-0.2},{-0.2,-0.1},{-0.1,0.0},{0.0,0.1},{0.1,0.2},{0.2,0.3},{0.3,0.4}};
        xMin_xMax_Optics_x = {{-0.3,-0.25},{-0.25,-0.2},{-0.2,-0.15},{-0.15,-0.1},{-0.1,-0.05},{-0.05,0.0},{0.0,0.05},{0.05,0.1},{0.1,0.15},{0.15,0.2},{0.2,0.25},{0.25,0.3},{0.3,0.35}};
        xMin_xMax_Optics_y ={{-0.10,-0.09},{-0.09,-0.08},{-0.08,-0.07},{-0.07,-0.06},{-0.06,-0.05},{-0.05,-0.04},{-0.04,-0.03},{-0.03,-0.02},{-0.02,-0.01},{-0.01,0.0},{0.0,0.01},{0.01,0.02},{0.02,0.03},{0.03,0.04},{0.04,0.05},{0.05,0.06},{0.06,0.07},{0.07,0.08},{0.08,0.09},{0.09,0.1}};

	}else{
	cout << "The kinematic: " << kine << " and magentic field setting: " << sbs_field << "is not implemented for the xMin_xMax initialization mode: " << cutvar_mode <<". Find the problem or make implementation!" << endl;
	}


  }else if(cutvar_mode == 1){
  //choosing mode 1 to handle the large slices, basically like integrals
  //primarily used for varying cut edges and quantifying any systematic change/effect

	//Second set of conditional for kinematic and magnetic field setting
        //Necessary to implement: SBS4 50%, SBS8 50%, 70%, 100%, SBS9 70%
        //All other kinematics and magnetic field settings not supported at this time
        if(kine == "SBS4" && sbs_field == 30){
        
	//BB preshower energy range
        xMin_xMax_ps = {{0.10,2.0},{0.12,2.0},{0.14,2.0},{0.16,2.0},{0.18,2.0},{0.20,2.0},{0.22,2.0},{0.24,2.0},{0.26,2.0},{0.28,2.0},{0.3,2.0}};	
	//number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{3.0,6.0},{4.0,6.0},{5.0,6.0}};		
	//BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,10.0},{0,11.0},{0,12.0},{0,13.0},{0,14.0},{0,15.0},{0,16.0},{0,17.0},{0,18.0},{0,19.0},{0,20.0}};
	//Target Vertex
	xMin_xMax_vz = {{-0.06,0.06},{-0.065,0.065},{-0.07,0.07},{-0.075,0.075},{-0.08,0.08}}; 
	//HCal Energy range
	xMin_xMax_hcal_e = {{0.035,1.0},{0.0375,1.0},{0.04,1.0},{0.0425,1.0},{0.045,1.0}};
	//xMin_xMax_hcal_e = {{0.01,1.0},{0.015,1.0},{0.02,1.0},{0.025,1.0},{0.03,1.0},{0.035,1.0},{0.04,1.0},{0.045,1.0},{0.05,1.0},{0.055,1.0},{0.06,1.0},{0.065,1.0},{0.07,1.0},{0.075,1.0},{0.08,1.0},{0.085,1.0},{0.09,1.0},{0.095,1.0},{0.1,1.0}};
	//HCal -Shower analog time diff range
	//xMin_xMax_coin_diff = {{-10.0,10.0},{-9.0,9.0},{-8.0,8.0},{-7.0,7.0},{-6.0,6.0}};
	//dy range
	xMin_xMax_dy = {{-0.36,0.36},{-0.34,0.34},{-0.32,0.32},{-0.3,0.3},{-0.28,0.28},{-0.26,0.26},{-0.24,0.24}};
	xMin_xMax_Optics_y = {{-0.1,0.1},{-0.095,0.095},{-0.09,0.09},{-0.085,0.085},{-0.08,0.08}};
		//Not everything requires a left right, if not just don't include
		if(left_right == "left"){
		//Left vary
		//Target Vertex
        	//xMin_xMax_vz = {{-0.095,0.075},{-0.085,0.075},{-0.075,0.075},{-0.065,0.075},{-0.055,0.075},{-0.045,0.075}};
		//W2 range
        	xMin_xMax_W2 = {{0.65,1.1},{0.7,1.1},{0.75,1.1},{0.8,1.1}};
		//E over p range
        	xMin_xMax_e_over_p = {{0.74,1.2},{0.76,1.2},{0.78,1.2},{0.8,1.2},{0.82,1.2},{0.84,1.2},{0.86,1.2}};
		//HCal - Shower analog time diff range
        	xMin_xMax_coin_diff = {{-12.0,10.0},{-11.0,10.0},{-10.0,10.0},{-9.0,10.0},{-8.0,10.0}};
		//dy range
        	//xMin_xMax_dy = {{-0.38,0.28},{-0.36,0.28},{-0.34,0.28},{-0.32,0.28},{-0.30,0.28},{-0.28,0.28},{-0.26,0.28},{-0.24,0.28},{-0.22,0.28},{-0.20,0.28},{-0.18,0.28}};
        	//Fid X, x expect range
        	//Bottom vary
       		xMin_xMax_x_exp = {{-1.6,0.79},{-1.55,0.79},{-1.5,0.79},{-1.45,0.79},{-1.4,0.79},{-1.35,0.79}};
		//Fid Y, y expect range
        	xMin_xMax_y_exp = {{-0.53,0.49},{-0.51,0.49},{-0.49,0.49},{-0.47,0.49},{-0.45,0.49}};
		xMin_xMax_Optics_x = {{-0.11,0.3},{-0.105,0.3},{-0.1,0.3},{-0.095,0.3},{-0.09,0.3}};
		}else if(left_right == "right"){
		//Right vary
		//Target Vertex
        	//xMin_xMax_vz = {{-0.075,0.045},{-0.075,0.055},{-0.075,0.065},{-0.075,0.075},{-0.075,0.085},{-0.075,0.095}};
		//W2 range
        	xMin_xMax_W2 = {{0.7,1.0},{0.7,1.05},{0.7,1.1},{0.7,1.15}};
		//E over p range
        	xMin_xMax_e_over_p = {{0.8,1.14},{0.8,1.16},{0.8,1.18},{0.8,1.20},{0.8,1.22},{0.8,1.24},{0.8,1.26}};
		 //HCal - Shower analog time diff range
        	xMin_xMax_coin_diff = {{-10.0,8.0},{-10.0,9.0},{-10.0,10.0},{-10.0,11.0},{-10.0,12.0}};
		//dy range
        	//xMin_xMax_dy = {{-0.28,0.18},{-0.28,0.20},{-0.28,0.22},{-0.28,0.24},{-0.28,0.26},{-0.28,0.28},{-0.28,0.30},{-0.28,0.32},{-0.28,0.34},{-0.28,0.36},{-0.28,0.38}};
        	//Fid X, x expect range
        	//Top vary
        	xMin_xMax_x_exp = {{-1.6,0.73},{-1.6,0.75},{-1.6,0.77},{-1.6,0.79},{-1.6,0.81},{-1.6,0.83},{-1.6,0.85}};
		//Fid Y, y expect range
        	xMin_xMax_y_exp ={{-0.49,0.45},{-0.49,0.47},{-0.49,0.49},{-0.49,0.51},{-0.49,0.53}};
		xMin_xMax_Optics_x = {{-0.1,0.29},{-0.1,0.295},{-0.1,0.3},{-0.1,0.305},{-0.1,0.31}};
		
		}else{
		cout << "The variable left_right: " << left_right << "was not properly implemented. Sort that out!" << endl;
		}

	}else if(kine == "SBS4" && sbs_field == 50){
	//BB preshower energy range
        xMin_xMax_ps = {{0.1,2.0},{0.12,2.0},{0.14,2.0},{0.16,2.0},{0.18,2.0},{0.20,2.0},{0.22,2.0},{0.24,2.0},{0.26,2.0},{0.28,2.0},{0.3,2.0}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{2.0,6.0},{3.0,6.0},{4.0,6.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,5},{0,10},{0,15},{0,16},{0,17},{0,18},{0,19},{0,20},{0,21},{0,22},{0,23},{0,24},{0,25}};
         //HCal Energy range
        xMin_xMax_hcal_e = {{0.01,1.0},{0.015,1.0},{0.02,1.0},{0.025,1.0},{0.03,1.0},{0.035,1.0},{0.04,1.0},{0.045,1.0},{0.05,1.0},{0.055,1.0}};

                //Not everything requires a left right, if not just don't include
                if(left_right == "left"){
                //Left vary
                //Target Vertex
                xMin_xMax_vz = {{-0.095,0.075},{-0.085,0.075},{-0.075,0.075},{-0.065,0.075},{-0.055,0.075},{-0.045,0.075}};
                //W2 range
                xMin_xMax_W2 = {{0.4,1.25},{0.45,1.25},{0.5,1.25},{0.55,1.25},{0.6,1.25},{0.65,1.25},{0.7,1.25},{0.75,1.25},{0.8,1.25}};
                //E over p range
                xMin_xMax_e_over_p = {{0.68,1.18},{0.70,1.18},{0.72,1.18},{0.74,1.18},{0.76,1.18},{0.78,1.18},{0.8,1.18},{0.82,1.18},{0.84,1.18},{0.86,1.18},{0.88,1.18}};
                //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-10.0,10.0},{-9.0,10.0},{-8.0,10.0},{-7.0,10.0},{-6.0,10.0},{-5.0,10.0},{-4.0,10.0}};
                //dy range
                xMin_xMax_dy = {{-0.38,0.28},{-0.36,0.28},{-0.34,0.28},{-0.32,0.28},{-0.30,0.28},{-0.28,0.28},{-0.26,0.28},{-0.24,0.28},{-0.22,0.28},{-0.20,0.28},{-0.18,0.28}};
                //Fid X, x expect range
                //Bottom vary
                xMin_xMax_x_exp = {{-1.6,0.79},{-1.55,0.79},{-1.5,0.79},{-1.45,0.79},{-1.4,0.79},{-1.35,0.79}};
                //Fid Y, y expect range
                xMin_xMax_y_exp = {{-0.59,0.49},{-0.57,0.49},{-0.55,0.49},{-0.53,0.49},{-0.51,0.49},{-0.49,0.49},{-0.47,0.49},{-0.45,0.49},{-0.43,0.49},{-0.41,0.49},{-0.39,0.49}};
		}else if(left_right == "right"){
                //Right vary
                //Target Vertex
                xMin_xMax_vz = {{-0.075,0.045},{-0.075,0.055},{-0.075,0.065},{-0.075,0.075},{-0.075,0.085},{-0.075,0.095}};
                //W2 range
                xMin_xMax_W2 = {{0.6,1.05},{0.6,1.10},{0.6,1.15},{0.6,1.2},{0.6,1.25},{0.6,1.3},{0.6,1.35},{0.6,1.4},{0.6,1.45}};
                //E over p range
                xMin_xMax_e_over_p = {{0.78,1.08},{0.78,1.10},{0.78,1.12},{0.78,1.14},{0.78,1.16},{0.78,1.18},{0.78,1.20},{0.78,1.22},{0.78,1.24},{0.78,1.26},{0.78,1.28}};
                 //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-10.0,4.0},{-10.0,5.0},{-10.0,6.0},{-10.0,7.0},{-10.0,8.0},{-10.0,9.0},{-10.0,10.0}};
                //dy range
                xMin_xMax_dy = {{-0.28,0.18},{-0.28,0.20},{-0.28,0.22},{-0.28,0.24},{-0.28,0.26},{-0.28,0.28},{-0.28,0.30},{-0.28,0.32},{-0.28,0.34},{-0.28,0.36},{-0.28,0.38}};
                //Fid X, x expect range
                //Top vary
                xMin_xMax_x_exp = {{-1.6,0.69},{-1.6,0.71},{-1.6,0.73},{-1.6,0.75},{-1.6,0.77},{-1.6,0.79},{-1.6,0.81},{-1.6,0.83},{-1.6,0.85},{-1.6,0.87},{-1.6,0.89}};
                //Fid Y, y expect range
                xMin_xMax_y_exp ={{-0.49,0.39},{-0.49,0.41},{-0.49,0.43},{-0.49,0.45},{-0.49,0.47},{-0.49,0.49},{-0.49,0.51},{-0.49,0.53},{-0.49,0.55},{-0.49,0.57},{-0.49,0.59}};
                }else{
                cout << "The variable left_right: " << left_right << "was not properly implemented. Sort that out!" << endl;
                }
	
	}else if(kine == "SBS8" && sbs_field == 50){
	//BB preshower energy range
        xMin_xMax_ps = {{0.10,2.6},{0.12,2.6},{0.14,2.6},{0.16,2.6},{0.18,2.6},{0.20,2.6},{0.22,2.6},{0.24,2.6},{0.26,2.6},{0.28,2.6},{0.3,2.6}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{3.0,6.0},{4.0,6.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,10.0},{0,11.0},{0,12.0},{0,13.0},{0,14.0},{0,15.0},{0,16.0},{0,17.0},{0,18.0},{0,19.0},{0,20.0}};
        //Target Vertex
        xMin_xMax_vz = {{-0.06,0.06},{-0.065,0.065},{-0.07,0.07},{-0.075,0.075},{-0.08,0.08}};
	//HCal Energy range
        xMin_xMax_hcal_e = {{0.05,1.0},{0.0525,1.0},{0.055,1.0},{0.0575,1.0},{0.06,1.0}};
	//dy range
        xMin_xMax_dy = {{-0.34,0.34},{-0.32,0.32},{-0.3,0.3},{-0.28,0.28},{-0.26,0.26}};

		//Not everything requires a left right, if not just don't include
                if(left_right == "left"){
                //Left vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.095,0.075},{-0.085,0.075},{-0.075,0.075},{-0.065,0.075},{-0.055,0.075},{-0.045,0.075}};
                //W2 range
                xMin_xMax_W2 = {{0.55,1.1},{0.6,1.1},{0.65,1.1},{0.7,1.1}};
                //E over p range
                xMin_xMax_e_over_p = {{0.76,1.2},{0.78,1.2},{0.8,1.2},{0.82,1.2},{0.84,1.2}};
                //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-12.0,10.0},{-11.0,10.0},{-10.0,10.0},{-9.0,10.0}};
                //dy range
                //xMin_xMax_dy = {{-0.38,0.28},{-0.36,0.28},{-0.34,0.28},{-0.32,0.28},{-0.30,0.28},{-0.28,0.28},{-0.26,0.28},{-0.24,0.28},{-0.22,0.28},{-0.20,0.28},{-0.18,0.28}};
                //Fid X, x expect range
                //Bottom vary
                xMin_xMax_x_exp = {{-2.36,0.83},{-2.34,0.83},{-2.32,0.83},{-2.30,0.83},{-2.28,0.83}};
                //Fid Y, y expect range
                xMin_xMax_y_exp = {{-0.55,0.51},{-0.53,0.51},{-0.51,0.51},{-0.49,0.51},{-0.47,0.51}};
                xMin_xMax_Optics_x = {{-0.16,0.3},{-0.155,0.3},{-0.15,0.3},{-0.145,0.3},{-0.14,0.3}};
		xMin_xMax_Optics_y = {{-0.1,0.09},{-0.095,0.09},{-0.09,0.09},{-0.085,0.09},{-0.08,0.09}};
		}else if(left_right == "right"){
                //Right vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.075,0.045},{-0.075,0.055},{-0.075,0.065},{-0.075,0.075},{-0.075,0.085},{-0.075,0.095}};
                //W2 range
                xMin_xMax_W2 = {{0.6,1.0},{0.6,1.05},{0.6,1.1},{0.6,1.15}};
                //E over p range
                xMin_xMax_e_over_p = {{0.8,1.16},{0.8,1.18},{0.8,1.20},{0.8,1.22},{0.8,1.24}};
                 //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-10.0,9.0},{-10.0,10.0},{-10.0,11.0},{-10.0,12.0}};
                //dy range
                //xMin_xMax_dy = {{-0.28,0.18},{-0.28,0.20},{-0.28,0.22},{-0.28,0.24},{-0.28,0.26},{-0.28,0.28},{-0.28,0.30},{-0.28,0.32},{-0.28,0.34},{-0.28,0.36},{-0.28,0.38}};
                //Fid X, x expect range
                //Top vary
                xMin_xMax_x_exp = {{-2.32,0.79},{-2.32,0.81},{-2.32,0.83},{-2.32,0.85},{-2.32,0.87}};
                //Fid Y, y expect range
                xMin_xMax_y_exp ={{-0.51,0.47},{-0.51,0.49},{-0.51,0.51},{-0.51,0.53},{-0.51,0.55}};
                xMin_xMax_Optics_x = {{-0.15,0.29},{-0.15,0.295},{-0.15,0.3},{-0.15,0.305},{-0.15,0.31}};
                xMin_xMax_Optics_y = {{-0.09,0.08},{-0.09,0.085},{-0.09,0.09},{-0.09,0.095},{-0.09,0.1}};
		}else{
                cout << "The variable left_right: " << left_right << "was not properly implemented. Sort that out!" << endl;
                }


	}else if(kine == "SBS8" && sbs_field == 70){

	//BB preshower energy range
        xMin_xMax_ps = {{0.10,2.6},{0.12,2.6},{0.14,2.6},{0.16,2.6},{0.18,2.6},{0.20,2.6},{0.22,2.6},{0.24,2.6},{0.26,2.6},{0.28,2.6},{0.3,2.6}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{3.0,6.0},{4.0,6.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,10.0},{0,11.0},{0,12.0},{0,13.0},{0,14.0},{0,15.0},{0,16.0},{0,17.0},{0,18.0},{0,19.0},{0,20.0}};
        //Target Vertex
        xMin_xMax_vz = {{-0.06,0.06},{-0.065,0.065},{-0.07,0.07},{-0.075,0.075},{-0.08,0.08}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.05,1.0},{0.0525,1.0},{0.055,1.0},{0.0575,1.0},{0.06,1.0}};
        //HCal -Shower analog time diff range
        //xMin_xMax_coin_diff = {{-10.0,10.0},{-9.0,9.0},{-8.0,8.0},{-7.0,7.0},{-6.0,6.0}};
        //dy range
        xMin_xMax_dy = {{-0.34,0.34},{-0.32,0.32},{-0.3,0.3},{-0.28,0.28},{-0.26,0.26}};

                //Not everything requires a left right, if not just don't include
                if(left_right == "left"){
                //Left vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.095,0.075},{-0.085,0.075},{-0.075,0.075},{-0.065,0.075},{-0.055,0.075},{-0.045,0.075}};
                //W2 range
                xMin_xMax_W2 = {{0.55,1.1},{0.6,1.1},{0.65,1.1},{0.7,1.1}};
                //E over p range
                xMin_xMax_e_over_p = {{0.76,1.2},{0.78,1.2},{0.8,1.2},{0.82,1.2},{0.84,1.2}};
                //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-12.0,10.0},{-11.0,10.0},{-10.0,10.0},{-9.0,10.0}};
                //dy range
                //xMin_xMax_dy = {{-0.38,0.28},{-0.36,0.28},{-0.34,0.28},{-0.32,0.28},{-0.30,0.28},{-0.28,0.28},{-0.26,0.28},{-0.24,0.28},{-0.22,0.28},{-0.20,0.28},{-0.18,0.28}};
                //Fid X, x expect range
                //Bottom vary
                xMin_xMax_x_exp = {{-2.36,0.83},{-2.34,0.83},{-2.32,0.83},{-2.3,0.83},{-2.28,0.83}};
                //Fid Y, y expect range
                xMin_xMax_y_exp = {{-0.55,0.51},{-0.53,0.51},{-0.51,0.51},{-0.49,0.51},{-0.47,0.51}};
		xMin_xMax_Optics_x = {{-0.16,0.25},{-0.155,0.25},{-0.15,0.25},{-0.145,0.25},{-0.14,0.25}};
		xMin_xMax_Optics_y = {{-0.1,0.09},{-0.095,0.09},{-0.09,0.09},{-0.085,0.09},{-0.08,0.09}};
		}else if(left_right == "right"){
                //Right vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.075,0.045},{-0.075,0.055},{-0.075,0.065},{-0.075,0.075},{-0.075,0.085},{-0.075,0.095}};
                //W2 range
                xMin_xMax_W2 = {{0.6,1.0},{0.6,1.05},{0.6,1.1},{0.6,1.15}};
                //E over p range
                xMin_xMax_e_over_p = {{0.8,1.16},{0.8,1.18},{0.8,1.20},{0.8,1.22},{0.8,1.24}};
                 //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-10.0,9.0},{-10.0,10.0},{-10.0,11.0},{-10.0,12.0}};
                //dy range
                //xMin_xMax_dy = {{-0.28,0.18},{-0.28,0.20},{-0.28,0.22},{-0.28,0.24},{-0.28,0.26},{-0.28,0.28},{-0.28,0.30},{-0.28,0.32},{-0.28,0.34},{-0.28,0.36},{-0.28,0.38}};
                //Fid X, x expect range
                //Top vary
                xMin_xMax_x_exp = {{-2.32,0.79},{-2.32,0.81},{-2.32,0.83},{-2.32,0.85},{-2.32,0.87}};
                //Fid Y, y expect range
                xMin_xMax_y_exp ={{-0.51,0.47},{-0.51,0.49},{-0.51,0.51},{-0.51,0.53},{-0.51,0.55}};
                xMin_xMax_Optics_x = {{-0.15,0.29},{-0.15,0.295},{-0.15,0.3},{-0.15,0.305},{-0.15,0.31}};
		xMin_xMax_Optics_y = {{-0.09,0.08},{-0.09,0.085},{-0.09,0.09},{-0.09,0.095},{-0.09,0.1}};
		}else{
                cout << "The variable left_right: " << left_right << "was not properly implemented. Sort that out!" << endl;
                }
	}else if(kine == "SBS8" && sbs_field == 100){
	//BB preshower energy range
        xMin_xMax_ps = {{0.10,2.6},{0.12,2.6},{0.14,2.6},{0.16,2.6},{0.18,2.6},{0.20,2.6},{0.22,2.6},{0.24,2.6},{0.26,2.6},{0.28,2.6},{0.3,2.6}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{3.0,6.0},{4.0,6.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,10.0},{0,11.0},{0,12.0},{0,13.0},{0,14.0},{0,15.0},{0,16.0},{0,17.0},{0,18.0},{0,19.0},{0,20.0}};
        //Target Vertex
        xMin_xMax_vz = {{-0.06,0.06},{-0.065,0.065},{-0.07,0.07},{-0.075,0.075},{-0.08,0.08}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.05,1.0},{0.0525,1.0},{0.055,1.0},{0.0575,1.0},{0.06,1.0}};
        //dy range
        xMin_xMax_dy = {{-0.34,0.34},{-0.32,0.32},{-0.3,0.3},{-0.28,0.28},{-0.26,0.26}};

                //Not everything requires a left right, if not just don't include
                if(left_right == "left"){
                //Left vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.095,0.075},{-0.085,0.075},{-0.075,0.075},{-0.065,0.075},{-0.055,0.075},{-0.045,0.075}};
                //W2 range
                xMin_xMax_W2 = {{0.55,1.1},{0.6,1.1},{0.65,1.1},{0.7,1.1}};
                //E over p range
                xMin_xMax_e_over_p = {{0.76,1.2},{0.78,1.2},{0.8,1.2},{0.82,1.2},{0.84,1.2}};
                //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-12.0,10.0},{-11.0,10.0},{-10.0,10.0},{-9.0,10.0}};
                //dy range
                //xMin_xMax_dy = {{-0.38,0.28},{-0.36,0.28},{-0.34,0.28},{-0.32,0.28},{-0.30,0.28},{-0.28,0.28},{-0.26,0.28},{-0.24,0.28},{-0.22,0.28},{-0.20,0.28},{-0.18,0.28}};
                //Fid X, x expect range
                //Bottom vary
                xMin_xMax_x_exp = {{-2.36,0.83},{-2.34,0.83},{-2.32,0.83},{-2.30,0.83},{-2.28,0.83}};
                //Fid Y, y expect range
                xMin_xMax_y_exp = {{-0.55,0.51},{-0.53,0.51},{-0.51,0.51},{-0.49,0.51},{-0.47,0.51}};
                xMin_xMax_Optics_x = {{-0.16,0.2},{-0.155,0.2},{-0.15,0.2},{-0.145,0.2},{-0.14,0.2}};
		xMin_xMax_Optics_y = {{-0.085,0.08},{-0.08,0.08},{-0.075,0.08},{-0.07,0.08}};
		}else if(left_right == "right"){
                //Right vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.075,0.045},{-0.075,0.055},{-0.075,0.065},{-0.075,0.075},{-0.075,0.085},{-0.075,0.095}};
                //W2 range
                xMin_xMax_W2 = {{0.6,1.0},{0.6,1.05},{0.6,1.1},{0.6,1.15}};
                //E over p range
                xMin_xMax_e_over_p = {{0.8,1.16},{0.8,1.18},{0.8,1.20},{0.8,1.22},{0.8,1.24}};
                 //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-10.0,9.0},{-10.0,10.0},{-10.0,11.0},{-10.0,12.0}};
                //dy range
                //xMin_xMax_dy = {{-0.28,0.18},{-0.28,0.20},{-0.28,0.22},{-0.28,0.24},{-0.28,0.26},{-0.28,0.28},{-0.28,0.30},{-0.28,0.32},{-0.28,0.34},{-0.28,0.36},{-0.28,0.38}};
                //Fid X, x expect range
                //Top vary
                xMin_xMax_x_exp = {{-2.32,0.79},{-2.32,0.81},{-2.32,0.83},{-2.32,0.85},{-2.32,0.87}};
                //Fid Y, y expect range
                xMin_xMax_y_exp ={{-0.51,0.47},{-0.51,0.49},{-0.51,0.51},{-0.51,0.53},{-0.51,0.55}};
                xMin_xMax_Optics_x = {{-0.15,0.19},{-0.15,0.195},{-0.15,0.2},{-0.15,0.205},{-0.15,0.21}};
		xMin_xMax_Optics_y = {{-0.08,0.07},{-0.08,0.075},{-0.08,0.08},{-0.08,0.085}};
		}else{

                cout << "The variable left_right: " << left_right << "was not properly implemented. Sort that out!" << endl;
                }
	}else if(kine == "SBS9" && sbs_field == 70){
		//BB preshower energy range
        xMin_xMax_ps = {{0.14,2.6},{0.16,2.6},{0.18,2.6},{0.20,2.6},{0.22,2.6},{0.24,2.6},{0.26,2.6}};
        //number of hits on BBGEM tracks range
        xMin_xMax_BBgem_nhits = {{3.0,6.0},{4.0,6.0},{5.0,6.0}};
        //BBGEM Track Chi2/ndf range
        xMin_xMax_BBgem_chi2ndf = {{0,10.0},{0,11.0},{0,12.0},{0,13.0},{0,14.0},{0,15.0},{0,16.0},{0,17.0},{0,18.0},{0,19.0},{0,20.0}};
        //Target Vertex
        xMin_xMax_vz = {{-0.06,0.06},{-0.065,0.065},{-0.07,0.07},{-0.075,0.075},{-0.08,0.08}};
        //HCal Energy range
        xMin_xMax_hcal_e = {{0.045,1.0},{0.0475,1.0},{0.05,1.0},{0.0525,1.0},{0.055,1.0}};
        //dy range
        xMin_xMax_dy = {{-0.34,0.34},{-0.32,0.32},{-0.3,0.3},{-0.28,0.28},{-0.26,0.26}};

                //Not everything requires a left right, if not just don't include
                if(left_right == "left"){
                //Left vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.095,0.075},{-0.085,0.075},{-0.075,0.075},{-0.065,0.075},{-0.055,0.075},{-0.045,0.075}};
                //W2 range
                xMin_xMax_W2 = {{0.65,1.1},{0.7,1.1},{0.75,1.1}};
                //E over p range
                xMin_xMax_e_over_p = {{0.76,1.2},{0.78,1.2},{0.8,1.2},{0.82,1.2},{0.84,1.2}};
                //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-12.0,10.0},{-11.0,10.0},{-10.0,10.0},{-9.0,10.0}};
                //dy range
                //xMin_xMax_dy = {{-0.38,0.28},{-0.36,0.28},{-0.34,0.28},{-0.32,0.28},{-0.30,0.28},{-0.28,0.28},{-0.26,0.28},{-0.24,0.28},{-0.22,0.28},{-0.20,0.28},{-0.18,0.28}};
                //Fid X, x expect range
                //Bottom vary
                xMin_xMax_x_exp = {{-2.36,0.84},{-2.34,0.84},{-2.32,0.84},{-2.30,0.84},{-2.28,0.84}};
                //Fid Y, y expect range
                xMin_xMax_y_exp = {{-0.54,0.5},{-0.52,0.5},{-0.5,0.5},{-0.48,0.5},{-0.46,0.5}};
                xMin_xMax_Optics_x = {{-0.16,0.3},{-0.155,0.3},{-0.15,0.3},{-0.145,0.3},{-0.14,0.3}};
                xMin_xMax_Optics_y = {{-0.095,0.09},{-0.09,0.09},{-0.085,0.09},{-0.08,0.09}};
		}else if(left_right == "right"){
                //Right vary
                //Target Vertex
                //xMin_xMax_vz = {{-0.075,0.045},{-0.075,0.055},{-0.075,0.065},{-0.075,0.075},{-0.075,0.085},{-0.075,0.095}};
                //W2 range
                xMin_xMax_W2 = {{0.65,1.0},{0.65,1.05},{0.65,1.1},{0.65,1.15}};
                //E over p range
                xMin_xMax_e_over_p = {{0.8,1.16},{0.8,1.18},{0.8,1.20},{0.8,1.22},{0.8,1.24}};
                 //HCal - Shower analog time diff range
                xMin_xMax_coin_diff = {{-10.0,9.0},{-10.0,10.0},{-10.0,11.0},{-10.0,12.0}};
                //dy range
                //xMin_xMax_dy = {{-0.28,0.18},{-0.28,0.20},{-0.28,0.22},{-0.28,0.24},{-0.28,0.26},{-0.28,0.28},{-0.28,0.30},{-0.28,0.32},{-0.28,0.34},{-0.28,0.36},{-0.28,0.38}};
                //Fid X, x expect range
                //Top vary
                xMin_xMax_x_exp = {{-2.32,0.80},{-2.32,0.82},{-2.32,0.84},{-2.32,0.86},{-2.32,0.88}};
                //Fid Y, y expect range
                xMin_xMax_y_exp ={{-0.5,0.46},{-0.5,0.48},{-0.5,0.5},{-0.5,0.52},{-0.5,0.54}};
                xMin_xMax_Optics_x = {{-0.15,0.29},{-0.15,0.295},{-0.15,0.3},{-0.15,0.305},{-0.15,0.31}};
                xMin_xMax_Optics_y = {{-0.09,0.08},{-0.09,0.085},{-0.09,0.09},{-0.09,0.095}};
                }else{

                cout << "The variable left_right: " << left_right << "was not properly implemented. Sort that out!" << endl;
                }
	}else{
        cout << "The kinematic: " << kine << " and magentic field setting: " << sbs_field << "is not implemented for the xMin_xMax initialization mode: " << cutvar_mode <<". Find the problem or make implementation!" << endl;
        }


  }else{
  //We only get here by error or if we are using a new slice mode which has not yet been implemented
  cout << "Error, the cutvar slice mode: " << cutvar_mode << " has not been implemented! Figure out why!" << endl;
  }


}//end initialize xMin_xMax
