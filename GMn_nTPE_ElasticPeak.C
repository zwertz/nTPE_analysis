#include <TMath.h>
#include <iostream>
#include <vector>
#include "analysis_utility_functions.h"

const Double_t M_p = 0.938272; // GeV
const Double_t M_n = 0.939565; // GeV
const Double_t M_e = 0.00051; // GeV
const Double_t ff = 0.05; // Added arbitrary smearing factor to account for beam energy fluctuations and fermi motion in downstream estimates
const Double_t hcal_width = 1.70434; // m
const Double_t hcal_sampfrac = 0.0795; // m
const Double_t hcal_threshconv = 6.914; // MeV/mV
const Double_t bbcal_threshconv = 7.2; // MeV/mV

void GMn_nTPE_ElasticPeak(TString Kinematic){
//Use the kinematic objectic class which keeps track of useful information. Load that information in via object member functions. This inherently looks in the config file for the information.
kinematic_obj myKin = New kinematic_obj("All_Kinematic.cfg",Kinematic);
double Ebeam = myKin.getBeamEnergy();
double BBAngle_rad = myKin.getBBAngle_Rad();
double BBDist = myKin.getBBDist();
double Q2 = myKin.getQ2();
double HCalAngle_rad = myKin.getHCalAngle_Rad();
double HCalDist = myKin.getHCalDist();

//determine min/max angle for HCal by trying to approximate the arclength
double HCalAngle_min = HCalAngle_rad - ((hcal_width/2)/HCalDist); //approx with arclength
double HCalAngle_max = HCalAngle_rad + ((hcal_width/2)/HCalDist); //approx with arclength

//dist to the front of BBCal
double sh_faceDist = 3.1 + BBDist; //1.2 to GEMs, another 1.9 to BBCal

double sh_yposition[7] = {-0.2565, -0.171, -.0855, 0.0, 0.0855, 0.171, 0.2565}; //Relative postion of shower columns 

//Initialize some arrays to some information
double effective_BBang[7] = {0};
double hcal_p_projang[7] = {0};
double hcal_n_projang[7] = {0};

double deltaBBAng = 0;

//determine the effective BBAngle based on trigonometry
for(int shcol=0;shcol<7;shcol++){
effective_BBang[shco]= (sh_yposition[shcol]/sh_faceDist)+BBAngle_rad;
}

double elasticPeak_e[7] = {0};
double elasticPeak_n[7] = {0};
double elasticPeak_p[7] = {0};

double epsilon_n[7]={0};
double epsilon_p[7]={0};

double nu_n[7] = {0};
double nu_p[7] = {0};

double p_e[7] = {0};
double p_n[7] = {0};
double p_p[7] = {0};

double hcal_minKE = 12345;
bool hcalON[7] = {false};

for(int shcol=0; shcol<7;shcol++){

double ePeak_ep = Ebeam/(1+(2*Ebeam/M_p)*(pow(sin(effective_BBang[shcol]/2),2)));
elasticPeak_e[shcol] = ePeak_ep;

//Calculate for proton
nu_p[shcol]=Ebeam -ePeak_ep;

double phel_p = sqrt(pow(nu_p[shcol],2)+2*M_p*nu_p[shcol]);

p_p[shcol]= phel_p;

epsilon_p[shcol] = pow(1+2*(1+Q2/4*M_p)*pow(tan(effective_BBang[shcol]/2),2),-1);
elasticPeak_p[shcol] = nu_p[shcol]+M_p;
hcal_p_projang[shcol] = acos((Ebeam - ePeak_ep*cos(effective_BBang[shcol]))/phel_p);

//Neutron
double ePeak_en = Ebeam/(1+(2*Ebeam/M_n)*(pow(sin(effective_BBang[shcol]/2),2)));
nu_n[shcol]=Ebeam -ePeak_en;

double phel_n = sqrt(pow(nu_n[shcol],2)+2*M_n*nu_n[shcol]);

p_n[shcol]= phel_n;

epsilon_n[shcol] = pow(1+2*(1+Q2/4*M_n)*pow(tan(effective_BBang[shcol]/2),2),-1);
elasticPeak_n[shcol] = nu_n[shcol]+M_n;
hcal_n_projang[shcol] = acos((Ebeam - ePeak_en*cos(effective_BBang[shcol]))/phel_n);

hcalON[shcol] = (hcal_p_projang[shcol]>hcal_minang) && (hcal_p_projang[shcol]<hcal_maxang) ;

//Check for lowest KE in HCal
if((nu_p[shcol] < hcal_minKE) && hcalON[shcol]){
hcal_minKE = nu_p[shcol];
}

if((nu_n[shcol] < hcal_minKE) && hcalON[shcol]){
hcal_minKE = nu_n[shcol];
}

}

cout << endl;

cout << "// BBSH col  |  e' angle  |  e' p   |  elastic proton KE  |  elastic neutron KE  |  had proj ang  | on HCal?  |  p epsilon  |  n epsilon  //" << endl; 

for(int i=0; i<7; i++){

std:cout << std::setprecision(5) << std::fixed;

cout << "//  "<< i << "    |  " << RadToDeg(effective_BBang[i]) << "  |  " << elasticPeak_e[i] << "  |  " << nu_p[i] << "    |  " << nu_n[i] << "    |   " << RadToDeg(hcal_p_projang[i]) << "   |  " << hcalON[i] << "   |  " << epsiolon_p[i] << "  |  " << epsilon_n[i] << "  //" << endl;

}

double hcal_ms = hcal_minKE * hcal_sampfrac*1000; //convert from GeV to MeV

cout << endl << endl << "Lowest energy sampled in HCal with estimated smearing (MeV): " << hcal_ms -hcal_ms*ff << "." << endl;
cout << endl << endl << "Suggested HCal threshold with estimated smearing: " (hcal_ms -hcal_ms*ff)/hcal_threshconv << "." << endl;

cout << endl << endl <<  "Lowest energy sampled in HCal (MeV): " << hcal_ms << "." << endl;
cout << endl << endl << "Suggest HCal threshold: " << hcal_ms/hcal_threshconv << "." << endl;

}
