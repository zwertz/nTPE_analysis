#ifndef CUTVAR_H
#define CUTVAR_H

//Author: Ezekiel Wertz
//A class to hold important quanties related to cut stability and systematic studies. Cutvar should only handle functions/processes on a case-by-case basis for data, mc p, mc n, and not all at once. A companion class to this is stability_analysis.h and stability_analysis.C. Will be modified as stability and systematic studies progress.

#include "TString.h"
#include <vector>
#include <utility>

class cutvar{
private:
TString CutVar,dxHistogramName,W2HistogramName,AxisTitle,CutString, DatavMC_flag;
vector<pair<double,double>> xMin_xMax_range;
double dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin, cut_hist_low, cut_hist_high, cut_hist_bin;
TH2D* dx_hist;
TH2D* W2_hist;

//// ranges for slices. These are what we want to change.
//// Started values taken from Maria, need to investigage each of these and figure out what I want
//In this case we are inherently assuming that the first in the pair is xMin and the second is XMax. A convention.
////****** dy ***************************************************************
vector<pair<double,double>> xMin_xMax_dy = {{-1,-0.5},{-0.5,0},{0,0.5},{0.5,1}};
////*****nsigdy***************************************************************
vector<pair<double,double>> xMin_xMax_nsigdy = {{-4.1,-2.1},{-3.1,-1.1},{-2,0.1},{-1,1},{0,2},{1.1,3},{2.1,4}};
///********vz**********************************************************************
vector<pair<double,double>> xMin_xMax_vz = {{-0.07,0.07},{-0.065,0.065},{-0.06,0.06},{-0.05,0.05},{-0.04,0.04},{-0.03,0.03}};
//// ******ps_e**********************************************************************
vector<pair<double,double>> xMin_xMax_ps = {{0.0,0.2},{0.2,0.4},{0.4,0.6},{0.6,0.8},{0.8,1.0},{1.2,1.4},{1.4,2.0}};
////******* hcal energy **************************************************************
vector<pair<double,double>> xMin_xMax_hcal_e = {{0.015,0.025},{0.020,0.030},{0.025,0.035},{0.030,0.040}};
//// fiducial x
//// fiducial y
//// ******w2****************************************************************
vector<pair<double,double>> xMin_xMax_w2 = {{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8},{0.8,0.9},{0.9,1.0},{1.0,1.1},{1.1,1.2},{1.2,1.3},{1.3,1.4},{1.4,1.5}};
// e_over_p
//Make educated guess when get here
vector<pair<double,double>> xMin_xMax_e_over_p;
//// *** x_expected*******************************************
vector<pair<double,double>> xMin_xMax_x_exp = {{-1.5,-1},{-1.25,-0.75},{-1,-0.5},{-0.75,-0.25},{-0.5,0},{-0.25,0.25},{0,0.5},{0.25,0.75},{0.5,1},{0.75,1.25}};
///****** nsigx_fid *******************************************
vector<pair<double,double>> xMin_xMax_nsigx_fid = {{0,1},{0.5,1.5},{1,2},{1.5,2.5},{2,3},{2.5,3.5},{3.5,4.5},{4,5},{5,6},{5.5,6.5},{6,7},{6.5,7.5}};
///******* y expected*************************************************************
vector<pair<double,double>> xMin_xMax_y_exp = {{-1,-0.5},{-0.75,-0.25},{-0.5,0},{-0.25,0.25},{0,0.5},{0.25,0.75},{0.5,1}};
//********* nsigy_fid ***************************************************************
vector<pair<double,double>> xMin_xMax_nsigy_fid = {{0,3},{0.5,3},{1,3},{1.5,3},{2,3}};
//********* coin_time ***************************************************************
//Make educated guess when get here
vector<pair<double,double>> xMin_xMax_coin_diff;
//********* spot_cuts ***************************************************************
//Maybe implement more ranges, to encompass all relevant cuts

//Function that will create a histogram of dx vs the cut variable. Relies on initialization of histo info from constructor.
TH2D* make2DdxCutHisto(TChain* C);
//Function that will create a histogram of W2 vs the cut variable. Relies on initialization of histo info from constructor.
TH2D* make2DW2CutHisto(TChain* C);

public:
//Constructor, will initialize a cutvar object. Which will be centrally used to handle a lot of the necessary info for the stability study. All private variables will be based on the cutvar input through conditionals.

cutvar(TString myVar,TString datCut, TString daFlag, double dx_low, double dx_high, double dx_bin, double W2_low, double W2_high, double W2_bin,TChain* C );

//Getter functions for each private variable. To provide controlle access. 

TString getdxHistoName();

TString getW2HistoName();

TString getAxisTitle();

TString getCutVar();

TString getDataMCFlag();

vector<pair<double,double>> getXMinXMaxRange();

TH2D* get2DdxCutHisto();

TH2D* get2DW2CutHisto();

vector<TH1D*> sliceAndProjectHisto_xMinxMax(TH2D* histo2D,TString xAxisName, TString yAxisName );
//Possibly more useful functions

};//end class
#endif
