//cuts.C
//Author: Ezekiel Wertz
//The companion implementation for the header file
//A file to handle implementations for the cuts we need or care about for analysis

#include "../include/cuts.h"

namespace cuts {
//functions to handle defining HCal active area
vector<double> hcal_ActiveArea_data(int num_blk_x, int num_blk_y, TString pass){

vector<double> hcalaa;
double hcalaa_Xi; 
double hcalaa_Xf;
double hcalaa_Yi;
double hcalaa_Yf;

	if(pass == "pass0"){
	//slightly different calibration, need different parameters
        hcalaa_Xi = exp_constants::hcalposXi_p0 + num_blk_x * exp_constants::hcalblk_h_p0;    
	hcalaa_Xf = exp_constants::hcalposXf_p0 - num_blk_x * exp_constants::hcalblk_h_p0;
	hcalaa_Yi = exp_constants::hcalposYi_p0 + num_blk_y * exp_constants::hcalblk_w_p0;
	hcalaa_Yf = exp_constants::hcalposYf_p0 - num_blk_y * exp_constants::hcalblk_w_p0;
	//better calibration, differerent parameters
	}else if(pass == "pass1" || pass == "pass2"){
	hcalaa_Xi = exp_constants::hcalposXi_mc + num_blk_x * exp_constants::hcalblk_div_v;      
	hcalaa_Xf = exp_constants::hcalposXf_mc - num_blk_x * exp_constants::hcalblk_div_v;
	hcalaa_Yi = exp_constants::hcalposYi_mc + num_blk_y * exp_constants::hcalblk_div_h;
	hcalaa_Yf = exp_constants::hcalposYf_mc - num_blk_y * exp_constants::hcalblk_div_h;
	//catch all case in the event we get something not a pass we can handle
	}else{
	//put an error message here
	cout << "Error HCal Active Area can't handle the pass: " << pass << "! Do something to fix it." << endl; 	
	}// end conditional

//store the values in a vector to use the information later
hcalaa.push_back(hcalaa_Xi);
hcalaa.push_back(hcalaa_Xf);
hcalaa.push_back(hcalaa_Yi);
hcalaa.push_back(hcalaa_Yf);

return hcalaa;
}//end function

vector<double> hcal_ActiveArea_MC(int num_blk_x, int num_blk_y){
vector<double> hcalaa;
double hcalaa_Xi = exp_constants::hcalposXi_mc + num_blk_x * exp_constants::hcalblk_div_v; 
double hcalaa_Xf = exp_constants::hcalposXf_mc - num_blk_x * exp_constants::hcalblk_div_v;
double hcalaa_Yi = exp_constants::hcalposYi_mc + num_blk_y * exp_constants::hcalblk_div_h;
double hcalaa_Yf = exp_constants::hcalposYf_mc - num_blk_y * exp_constants::hcalblk_div_h;

//store the values in a vector to use the information later
hcalaa.push_back(hcalaa_Xi);
hcalaa.push_back(hcalaa_Xf);
hcalaa.push_back(hcalaa_Yi);
hcalaa.push_back(hcalaa_Yf);

return hcalaa;
}//end function

//check position per event against hcal active area. Is True if detected on active area
bool hcalaa_ON (double xhcal, double yhcal, vector<double> hcalaa){

double hcalx_top = hcalaa[0];
double hcalx_bot = hcalaa[1];
double hcaly_left = hcalaa[2];
double hcaly_right = hcalaa[3];

//actually make the boolean value. Could be either true or false
bool detect = (yhcal > hcaly_left) && (yhcal < hcaly_right) && (xhcal > hcalx_top) && (xhcal < hcalx_bot); 

return detect;
}//end function

//function to define HCal fiducial area
vector<double> hcalfid(double dxsig_p, double dxsig_n, double dysig, vector<double> hcalaa,double num_sig_x, double num_sig_y ){
vector<double> hcalfid;

double hcalx_top = hcalaa[0] + num_sig_x * dxsig_p; //-x margin (relevant for proton)
double hcalx_bot = hcalaa[1] - num_sig_x * dxsig_n; //+x margin (relevant for neutron)
double hcaly_left = hcalaa[2] + num_sig_y * dysig; // -y margin
double hcaly_right = hcalaa[3] - num_sig_y * dysig; // +y margin

hcalfid.push_back(hcalx_top);
hcalfid.push_back(hcalx_bot);
hcalfid.push_back(hcaly_left);
hcalfid.push_back(hcaly_right);

return hcalfid;
}//end function

//function determine if event is inside the fiducial region. That is equally count neutrons and protons
bool hcalfid_IN(double xhcal_expect, double yhcal_expect, double dx_pn, vector<double> hcalfid){

double hcalx_top = hcalfid[0]; //-x fid
double hcalx_bot = hcalfid[1]; //+x fid
double hcaly_left = hcalfid[2]; //-y fid
double hcaly_right = hcalfid[3]; //+y fid

//for proton hypothesis
double xhcal_exp_pro = xhcal_expect - dx_pn; //define the proton expect postion based on the difference between the proton and neutron peaks

bool inFid = (yhcal_expect > hcaly_left) && (yhcal_expect < hcaly_right) &&  //check that dy is the same for protons and neutrons
	     (xhcal_expect > hcalx_top) && (xhcal_expect < hcalx_bot) &&  //check that expected position if it were a neutron
	     (xhcal_exp_pro > hcalx_top) && (xhcal_exp_pro < hcalx_bot);  //check that expected position if it were a proton

return inFid;
}//end function

//Funtion that defines good W2 elastic cut
bool goodW2(double W2, double W2_low, double W2_high){
bool goodW2 = (W2 >= W2_low) && (W2 <= W2_high);
return goodW2;
}



}//end namespace

