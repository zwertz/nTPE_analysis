#ifndef CUTS_H
#define CUTS_H

#include <vector>
#include "TString.h"
#include <iostream>

namespace cuts{

//For HCal Active Area cut and data
vector<double> hcal_ActiveArea_data(int num_blk_x, int num_blk_y, TString pass);//number of blocks to exclude from top and bottom, number of blocks to exclude from left and right

//For HCal Active Area cut and MC
vector<double> hcal_ActiveArea_MC(int num_blk_x, int num_blk_y);//number of blocks to exclude from top and bottom, number of blocks to exclude from left and right

//Function that checks per event if the position is on HCal. Returns True if detected
bool hcalaa_ON (double xhcal, double yhcal, vector<double> hcalaa);


//Function that defines HCal Fiducial region
vector<double> hcalfid(double dxsig_p, double dxsig_n, double dysig, vector<double> hcalaa,double num_sig_x, double num_sig_y );

//Funtion that defines good W2 elastic cut
bool goodW2(double W2, double W2_low, double W2_high);

}//end namespace
#endif
