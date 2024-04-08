#ifndef CUTS_H
#define CUTS_H

#include <vector>
#include "TString.h"
#include <iostream>
#include "TTree.h"

namespace cuts{

//function to define HCal physical position boundaries
vector<double> hcal_Position_data(TString pass);
vector<double> hcal_Position_MC();

//For HCal Active Area cut and data
vector<double> hcal_ActiveArea_data(int num_blk_x, int num_blk_y, TString pass);//number of blocks to exclude from top and bottom, number of blocks to exclude from left and right

//For HCal Active Area cut and MC
vector<double> hcal_ActiveArea_MC(int num_blk_x, int num_blk_y);//number of blocks to exclude from top and bottom, number of blocks to exclude from left and right

//Function that checks per event if the position is on HCal. Returns True if detected
bool hcalaa_ON (double xhcal, double yhcal, vector<double> hcalaa);


//Function that defines HCal Fiducial region
vector<double> hcalfid(double dxsig_p, double dxsig_n, double dysig, vector<double> hcalaa,double num_sig_x, double num_sig_y );

//function determine if event is inside the fiducial region. That is equally count neutrons and protons
bool hcalfid_IN(double xhcal_expect, double yhcal_expect, double dx_pn, vector<double> hcalfid);

//Funtion that defines good W2 elastic cut
bool goodW2(double W2, double W2_low, double W2_high);

//Function used in defining globabl cut
bool failedGlobal(TTreeFormula *GlobalCut);

//Function used to define good coin cut
bool passCoin(double coin_bestclus,double coin_mean, double coin_sig_fac, double coin_sigma);

//Function to define a good dy cut
bool good_dy(double dy_bestclus,double dyO_p, double dysig_cut_fac, double dysig_p);

//Function to check if above min HCal E value
bool passHCalE(double hcal_e,double hcalemin);

//Function to check if number of cluster in HCal is above a min
bool passHCal_NClus(double nclus_hcal,int hcalnclusmin);



}//end namespace
#endif
