//cutvar.C
//Author: Ezekiel Wertz
//Companion implementation. A class to hold important quanties related to cut stability and systematic studies. utvar should only handle functions/processes on a case-by-case basis for data, mc p, mc n, and not all at once. A companion class to this is stability_analysis.h and stability_analysis.C. Will be modified as stability and systematic studies progress.

#include "../include/cutvar.h"

//public constructor implementation. Constructor, will initialize a cutvar object. Which will be centrally used to handle a lot of the necessary info for the stability study. All private variables will be based on the cutvar input through conditionals.

cutvar::cutvar(TString myVar,TString datCut,TString daFlag, double dx_low, double dx_high, double dx_bin, double W2_low, double W2_high, double W2_bin,TChain* C ){

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

	//For each parameter we want to evaluate, make a conditional. Then initialize all private variables. This will need to be expanded as cut variables are investigated.
	if(myVar == "dy"){
	dxHistogramName = "hcal_dx__hcal_dy";
	W2HistogramName = "W2__hcal_dy";
	AxisTitle = "hcal_dy";
	xMin_xMax_range = xMin_xMax_dy; 
	cut_hist_low = -2;
	cut_hist_high = 2;
	cut_hist_bin =  100;
	}else if(myVar == "nsigdy"){
        dxHistogramName = "hcal_dx__hcal_nsigdy";
        W2HistogramName = "W2__hcal_nsigdy";
	AxisTitle = "hcal_nsigdy";
        xMin_xMax_range = xMin_xMax_nsigdy;
        cut_hist_low = -5;
        cut_hist_high = 5;
        cut_hist_bin = 100;
	}else if(myVar == "vz"){
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
        cut_hist_high = 2.2;
        cut_hist_bin = 250;
	}else if(myVar == "hcal_e"){
        dxHistogramName = "hcal_dx__hcal_e";
        W2HistogramName = "W2__hcal_e";
	AxisTitle = "hcal energy";
        xMin_xMax_range = xMin_xMax_hcal_e;
        cut_hist_low = 0;
        cut_hist_high = 2;
        cut_hist_bin = 2000;
	}else if(myVar == "w2"){
        dxHistogramName = "hcal_dx__W2";
        W2HistogramName = "W2__W2";
	AxisTitle = "W^{2}";
        xMin_xMax_range = xMin_xMax_w2;
        cut_hist_low = W2_low;
        cut_hist_high = W2_high;
        cut_hist_bin = W2_bin;
	}else if(myVar == "x_exp"){
        dxHistogramName = "hcal_dx__hcal_x_exp";
        W2HistogramName = "W2__hcal_x_exp";
	AxisTitle = "x expected";
        xMin_xMax_range = xMin_xMax_x_exp;
        cut_hist_low = -3;
        cut_hist_high = 3;
        cut_hist_bin = 600;
	}else if(myVar == "nsigx_fid"){
        dxHistogramName = "hcal_dx__nsigx_fid";
        W2HistogramName = "W2__nsigx_fid";
	AxisTitle = "nsigx_fid";
        xMin_xMax_range = xMin_xMax_nsigx_fid;
        cut_hist_low = -10;
        cut_hist_high = 10;
        cut_hist_bin = 200;
	}else if(myVar == "y_exp"){
        dxHistogramName = "hcal_dx__y_exp";
        W2HistogramName = "W2__y_exp";
	AxisTitle = "y expected";
        xMin_xMax_range = xMin_xMax_y_exp;
        cut_hist_low = -2;
        cut_hist_high = 2;
        cut_hist_bin = 400;
	}else if(myVar == "nsigy_fid"){
        dxHistogramName = "hcal_dx__nsigy_fid";
        W2HistogramName = "W2__nsigy_fid";
	AxisTitle = "nsigy_fid";
        xMin_xMax_range = xMin_xMax_nsigy_fid;
        cut_hist_low = -1;
        cut_hist_high = 4;
        cut_hist_bin = 50;
	}else if(myVar == "e_over_p"){
	dxHistogramName = "hcal_dx__e_over_p";
	W2HistogramName = "W2__e_over_p";
	AxisTitle = "e_over_p";
	xMin_xMax_range = xMin_xMax_e_over_p;
	cut_hist_low = -1;
	cut_hist_high = 3;
	cut_hist_bin = 400;
	}else if(myVar == "hcal_shower_atime_diff"){
	dxHistogramName = "hcal_dx__hcal_shower_atime_diff";
	W2HistogramName = "W2__hcal_shower_atime_diff";
	AxisTitle = "hcal_shower_atime_diff";
	xMin_xMax_range = xMin_xMax_coin_diff;
	cut_hist_low = -15;
	cut_hist_high = 15;
	cut_hist_bin = 300;
	}else{
	cout << "Error in cut variable string " << myVar << "." << endl;
	return;
	}

	//Now create the cut histograms we care about
	dx_hist = cutvar::make2DdxCutHisto(C);

	//Modify the histogram names for data or mc depending on the input parameter.
        TString newW2HistoName = DatavMC_flag +"_"+ W2HistogramName;


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
	dx_hist = get2DdxCutHisto();
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
	C->Draw(Form("dx:%s>>%s(%f,%f,%f,%f,%f,%f)",CutVar.Data(),dxHistogramName.Data(),cut_hist_bin,cut_hist_low,cut_hist_high,dx_hist_bin,dx_hist_low,dx_hist_high),Form("Final_MC_weight*(%s)",CutString.Data()),"COLZ");
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
}

//Function that will create a histogram of W2 vs the cut variable. Relies on initialization of histo info from constructor.
//TH2D* cutvar::make2DW2CutHisto(TChain* C);

TString cutvar::getdxHistoName(){ return dxHistogramName; }

TString cutvar::getW2HistoName(){ return W2HistogramName; }

TString cutvar::getAxisTitle(){return AxisTitle; }

TString cutvar::getCutVar(){return CutVar; }

TString cutvar::getDataMCFlag(){return DatavMC_flag;}

TString cutvar::getCutString(){return CutString;}

double cutvar::get_dx_hist_low(){return dx_hist_low;}

double cutvar::get_dx_hist_high(){return dx_hist_high;}

double cutvar::get_dx_hist_bin(){return dx_hist_bin;}

double cutvar::get_W2_hist_low(){return W2_hist_low;}

double cutvar::get_W2_hist_high(){return W2_hist_high;}

double cutvar::get_W2_hist_bin(){return W2_hist_bin;}

double cutvar::get_cut_hist_low(){return cut_hist_low;}

double cutvar::get_cut_hist_high(){return cut_hist_high;}

double cutvar::get_cut_hist_bin(){return cut_hist_bin;}

vector<pair<double,double>> cutvar::getXMinXMaxRange(){return xMin_xMax_range; }

TH2D* cutvar::get2DdxCutHisto(){return dx_hist;}

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
}
