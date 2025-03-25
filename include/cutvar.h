#ifndef CUTVAR_H
#define CUTVAR_H

//Author: Ezekiel Wertz
//A class to hold important quanties related to cut stability and systematic studies. Cutvar should only handle functions/processes on a case-by-case basis for data, mc p, mc n, and not all at once. A companion class to this is stability_analysis.h and stability_analysis.C. Will be modified as stability and systematic studies progress.

#include "TString.h"
#include <vector>
#include <utility>
#include <iostream>

class cutvar{
private:
TString CutVar,dxHistogramName,W2HistogramName,AxisTitle,CutString, DatavMC_flag,kine, left_right;
vector<pair<double,double>> xMin_xMax_range;
double dx_hist_low, dx_hist_high, dx_hist_bin, W2_hist_low, W2_hist_high, W2_hist_bin, cut_hist_low, cut_hist_high, cut_hist_bin;
TH2D* dx_hist;
TH2D* W2_hist;
int sbs_field;
int cutvar_mode; //If equals 0, then mode is small slices like differential. If equal 1, then mode is large slices like integral.


//// ranges for slices. These are what we want to change.
//In this case we are inherently assuming that the first in the pair is xMin and the second is XMax. A convention.
////****** dy ***************************************************************
vector<pair<double,double>> xMin_xMax_dy;
///********vz**********************************************************************
vector<pair<double,double>> xMin_xMax_vz ;
///********BBgem_nhits**********************************************************************
vector<pair<double,double>> xMin_xMax_BBgem_nhits;
//// ******ps_e**********************************************************************
vector<pair<double,double>> xMin_xMax_ps;
////******* hcal energy **************************************************************
vector<pair<double,double>> xMin_xMax_hcal_e;
//// ******w2****************************************************************
vector<pair<double,double>> xMin_xMax_W2;
// e_over_p
vector<pair<double,double>> xMin_xMax_e_over_p;
//Fid cut for x expect and y expect
//// *** x_expected*******************************************
vector<pair<double,double>> xMin_xMax_x_exp;
///******* y expected*************************************************************
vector<pair<double,double>> xMin_xMax_y_exp;
//********* coin_time ***************************************************************
vector<pair<double,double>> xMin_xMax_coin_diff;
/////******* track chi2ndf ***********************************************************
vector<pair<double,double>> xMin_xMax_BBgem_chi2ndf;
/////******* Optics X ***********************************************************
vector<pair<double,double>> xMin_xMax_Optics_x;
/////******* Optics Y ***********************************************************
vector<pair<double,double>> xMin_xMax_Optics_y;


//Function that will create a histogram of dx vs the cut variable. Relies on initialization of histo info from constructor.
TH2D* make2DdxCutHisto(TChain* C);
//Function that will create a histogram of W2 vs the cut variable. Relies on initialization of histo info from constructor.
TH2D* make2DW2CutHisto(TChain* C);

//Function that handles the initialization of the xMinXMax vectors dependent on the slice mode, the kinematic, and the magnetic field setting
void initialize_xMin_xMax();

public:
//Constructor, will initialize a cutvar object. Which will be centrally used to handle a lot of the necessary info for the stability study. All private variables will be based on the cutvar input through conditionals.

cutvar(TString myVar,TString datCut, TString daFlag,int daMode,TString Left_Right, TString daKine,int daField,double dx_low, double dx_high, double dx_bin, double W2_low, double W2_high, double W2_bin,TChain* C );

//Destructor
//Might have dynamically allocated memory or pointers. May have to worry about that.
~cutvar();

//Copy constructor that takes a reference as input
cutvar(cutvar &myCut);

//Getter functions for each private variable. To provide controlle access. 

TString getdxHistoName();

TString getW2HistoName();

TString getAxisTitle();

TString getCutVar();

TString getCutString();

TString getDataMCFlag();

TString getKine();

TString get_left_right();

double get_dx_hist_low();

double get_dx_hist_high();

double get_dx_hist_bin();

double get_W2_hist_low();

double get_W2_hist_high();

double get_W2_hist_bin();

double get_cut_hist_low();

double get_cut_hist_high();

double get_cut_hist_bin();

int get_cutvar_mode();

int get_sbs_field();

vector<pair<double,double>> getXMinXMaxRange();

TH2D* get2DdxCutHisto();

TH2D* get2DW2CutHisto();

vector<TH1D*> sliceAndProjectHisto_xMinxMax(TH2D* histo2D,TString xAxisName, TString yAxisName );
//Possibly more useful functions

};//end class
#endif
