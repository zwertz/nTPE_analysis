//cutvar.C
//Author: Ezekiel Wertz
//Companion implementation. A class to hold important quanties related to cut stability and systematic studies. Will be modified as stability and systematic studies progress.

#include "../include/cutvar.h"
#include <iostream>

//public constructor implementation. Constructor, will initialize a cutvar object. Which will be centrally used to handle a lot of the necessary info for the stability study. All private variables will be based on the cutvar input through conditionals.

cutvar::cutvar(TString myVar){

	//Store the cut variable independantly
	CutVar = myVar;

	//For each parameter we want to evaluate, make a conditional. Then initialize all private variables. This will need to be expanded as cut variables are investigated.
	if(myVar == "dy"){
	HistogramName = "hcal_dx_hcal_dy";
	AxisTitle = "hcal_dy";
	xMin_xMax_range = xMin_xMax_dy; 
	}else if(myVar == "nsigdy"){
        HistogramName = "hcal_dx_hcal_nsigdy";
        AxisTitle = "hcal_nsigdy";
        xMin_xMax_range = xMin_xMax_nsigdy;
        }else if(myVar == "vz"){
        HistogramName = "hcal_dx_tr_vz";
        AxisTitle = "track vz";
        xMin_xMax_range = xMin_xMax_vz;
        }else if(myVar == "ps"){
        HistogramName = "hcal_dx_ps_e";
        AxisTitle = "preshower energy";
        xMin_xMax_range = xMin_xMax_ps;
        }else if(myVar == "hcal_e"){
        HistogramName = "hcal_dx_hcal_e";
        AxisTitle = "hcal energy";
        xMin_xMax_range = xMin_xMax_hcal_e;
        }else if(myVar == "w2"){
        HistogramName = "hcal_dx_W2";
        AxisTitle = "W^{2}";
        xMin_xMax_range = xMin_xMax_w2;
        }else if(myVar == "x_exp"){
        HistogramName = "hcal_dx_hcal_x_exp";
        AxisTitle = "x expected";
        xMin_xMax_range = xMin_xMax_x_exp;
        }else if(myVar == "nsigx_fid"){
        HistogramName = "hcal_dx_nsigx_fid";
        AxisTitle = "nsigx_fid";
        xMin_xMax_range = xMin_xMax_nsigx_fid;
        }else if(myVar == "y_exp"){
        HistogramName = "hcal_dx_y_exp";
        AxisTitle = "y expected";
        xMin_xMax_range = xMin_xMax_y_exp;
        }else if(myVar == "nsigy_fid"){
        HistogramName = "hcal_dx_nsigy_fid";
        AxisTitle = "nsigy_fid";
        xMin_xMax_range = xMin_xMax_nsigy_fid;
        }else{
	cout << "Error in cut variable string " << myVar << "." << endl;
	return;
	}

}//end constructor

TString cutvar::getHistoName(){ return HistogramName; }

TString cutvar::getAxisTitle(){return AxisTitle; }

TString cutvar::getCutVar(){return CutVar; }

vector<pair<double,double>> cutvar::getXMinXMaxRange(){return xMin_xMax_range; }
