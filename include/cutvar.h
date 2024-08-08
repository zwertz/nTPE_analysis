#ifndef CUTVAR_H
#define CUTVAR_H

//Author: Ezekiel Wertz
//A class to hold important quanties related to cut stability and systematic studies. Will be modified as stability and systematic studies progress.

#include "TString.h"
#include <vector>
#include <utility>

class cutvar{
private:
TString CutVar,HistogramName,AxisTitle;
vector<pair<double,double>> xMin_xMax_range;

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
vector<pair<double,double>> xMin_xMax_ps = {{0.1,0.2},{0.12,0.22},{0.14,0.24},{0.16,0.26},{0.18,0.28},{0.2,0.3},{0.22,0.32},{0.24,0.34},{0.26,0.36},{0.28,0.38}};
////******* hcal energy **************************************************************
vector<pair<double,double>> xMin_xMax_hcal_e = {{0.015,0.025},{0.020,0.030},{0.025,0.035},{0.030,0.040}};
//// fiducial x
//// fiducial y
//// ******w2****************************************************************
vector<pair<double,double>> xMin_xMax_w2 = {{0.4,0.5},{0.5,0.6},{0.6,0.7},{0.7,0.8},{0.8,0.9},{0.9,1.0},{1.0,1.1},{1.1,1.2},{1.2,1.3},{1.3,1.4},{1.4,1.5}};
// e_over_p
//// *** x_expected*******************************************
vector<pair<double,double>> xMin_xMax_x_exp = {{-1.5,-1},{-1.25,-0.75},{-1,-0.5},{-0.75,-0.25},{-0.5,0},{-0.25,0.25},{0,0.5},{0.25,0.75},{0.5,1},{0.75,1.25}};
///****** nsigx_fid *******************************************
vector<pair<double,double>> xMin_xMax_nsigx_fid = {{0,1},{0.5,1.5},{1,2},{1.5,2.5},{2,3},{2.5,3.5},{3.5,4.5},{4,5},{5,6},{5.5,6.5},{6,7},{6.5,7.5}};
///******* y expected*************************************************************
vector<pair<double,double>> xMin_xMax_y_exp = {{-1,-0.5},{-0.75,-0.25},{-0.5,0},{-0.25,0.25},{0,0.5},{0.25,0.75},{0.5,1}};
//********* nsigy_fid ***************************************************************
vector<pair<double,double>> xMin_xMax_nsigy_fid = {{0,3},{0.5,3},{1,3},{1.5,3},{2,3}};
//********* coin_time ***************************************************************

//********* spot_cuts ***************************************************************
//Maybe implement more ranges, to encompass all relevant cuts


public:
//Constructor, will initialize a cutvar object. Which will be centrally used to handle a lot of the necessary info for the stability study. All private variables will be based on the cutvar input through conditionals.

//Getter functions for each private variable. To provide controlle access. 

cutvar(TString myVar);

TString getHistoName();

TString getAxisTitle();

TString getCutVar();

vector<pair<double,double>> getXMinXMaxRange();

//Possibly more useful functions

};//end class
#endif
