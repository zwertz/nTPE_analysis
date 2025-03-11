//cutvar.C
//Author: Ezekiel Wertz
//Companion implementation. A class to hold important quanties related to cut stability and systematic studies. utvar should only handle functions/processes on a case-by-case basis for data, mc p, mc n, and not all at once. A companion class to this is stability_analysis.h and stability_analysis.C. Will be modified as stability and systematic studies progress.

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
        cut_hist_high = 2.2;
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
	cut_hist_low = -13;
	cut_hist_high = 13;
	cut_hist_bin = 200;
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
	dx_hist = myCut.get2DdxCutHisto();
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
}

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

	}else{
        cout << "The kinematic: " << kine << " and magentic field setting: " << sbs_field << "is not implemented for the xMin_xMax initialization mode: " << cutvar_mode <<". Find the problem or make implementation!" << endl;
        }


  }else{
  //We only get here by error or if we are using a new slice mode which has not yet been implemented
  cout << "Error, the cutvar slice mode: " << cutvar_mode << " has not been implemented! Figure out why!" << endl;
  }


}//end initialize xMin_xMax
