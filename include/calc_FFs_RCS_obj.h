#ifndef CALC_FFS_RCS_OBJ.H
#define CALC_FFS_RCS_OBJ.H

#include "TString.h"
#include <vector>
#include "physics_constants.h"
#include "../../src/kinematic_obj.C"


//Author: Ezekiel Wertz
//A class to facilite calculating electromagnetic form factors via various supported parameterization. The most critical being what is used in SIMC. The class will also calculate the reduced cross-section. The FFs and the Reduced Cross-Section will be dependent on user-input which will allow for multiple supposed parameterizations. Initially this will include Kelly, Riordan, and Ye. The class object will hold these all as class variables. Public and private functions and the constructor will enable this class to be compatible with other scripts and classes in this analysis framework.

//Supported Parameterizations: Ye, Kelly, SIMC (Kelly and Riordan)

//Kelly Parameterization Citation
//(Phys. Rev. C 70, 068202)

//Riordan Parameterization Citation
//(Phys. Rev. Lett. 105, 262302)

//Ye Parameterization Citation
//( Physics Letters B 777, 8 (2018))

class calc_FFs_RCS_obj{
private:
//Will need to add class variables
TString param_type, kin, kinematic_file;
double reduced_cross_section_n, reduced_cross_section_p, reduced_cross_section_ratio, reduced_cross_section_n_err, reduced_cross_section_p_err, reduced_cross_section_ratio_err;
double Ebeam, Eprime, BB_angle_rad, tau_p, tau_n, Q2, epsilon_p, epsilon_n;
double GEp, GEp_Div_GD, GMp, GMp_Div_mup, GMp_Div_mup_GD, GEn, GEn_Div_GD, GMn, GMn_Div_mun, GMn_Div_mun_GD;
double GEp_err, GEp_Div_GD_err, GMp_err, GMp_Div_mup_err, GMp_Div_mup_GD_err, GEn_err, GEn_Div_GD_err, GMn_err, GMn_Div_mun_err, GMn_Div_mun_GD_err;
double GD_const = 0.71;
//For Kelly we assume the vector has the format a_1, b_1, b_2, b_3
//For GEp
vector<double> Kelly_GEp = {-0.24, 10.98, 12.82, 21.97};
vector<double> Kelly_GEp_err = {0.12, 0.19, 1.1, 6.8}; 
//For GMp/mu_p
vector<double> Kelly_GMp = {0.12, 10.97,1 8.86, 6.55};
vector<double> Kelly_GMp_err = {0.04, 0.11, 0.28, 1.2};
//For GMn/mu_n
vector<double> Kelly_GMn = {2.33, 14.72, 24.20, 84.10};
vector<double> Kelly_GMn_err = {1.4, 1.7, 9.8, 41 };
//GEn will be from Galster Parameterization
double Kelly_A = 1.70;
double Kelly_A_err = 0.04;
double Kelly_B = 3.30;
double Kelly_B_err = 0.32;
//Kelly has a dependent error analysis. Would need to implement a covariance matrix analysis. However, since this is for Rsim. We should not implement an uncertainty.


//For Ye we assume the vector has the format a_0, a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8, a_9, a_10, a_11, a_12
//For GEp
vector<double> GEp_Ye_coef = {0.239163298067, -1.109858574410, 1.444380813060, 0.479569465603, -2.286894741870, 1.126632984980, 1.250619843540, -3.631020471590, 4.082217023790, 0.504097346499, -5.085120460510, 3.967742543950, -0.981529071103};
vector<double> GEp_Ye_err_coef = {-1.97750308, -4.46566998*pow(10,-1), 2.94508717*pow(10,-1), 1.54467525, 9.05268347*pow(10,-1), -6.00008111*pow(10,-1), -1.10732394, -9.85982716*pow(10,-2), 4.63035988*pow(10,-1), 1.37729116*pow(10,-1), -7.82991627*pow(10,-2), -3.63056932*pow(10,-2), 2.64219326*pow(10,-3), 3.13261383*pow(10,-3), 3.89593858*pow(10,-4)};
//For GMp/mu_p
vector<double> GMp_Div_mup_Ye_coef = {0.264142994136, -1.095306122120, 1.218553781780, 0.661136493537, -1.405678925030, -1.356418438880, 1.447029155340, 4.235669735900, -5.334045653410, -2.916300520960, 8.707403067570, -5.706999943750, 1.280814375890};
vector<double> GMp_Div_mup_Ye_err_coef = {-1.76549673, 1.67218457*pow(10,-1), -1.20542733, -4.72244127*pow(10,-1), 1.41548871, 6.61320779*pow(10,-1), -8.16422909*pow(10,-1), -3.73804477*pow(10,-1), 2.62223992*pow(10,-1), 1.28886639*pow(10,-1), -3.90901510*pow(10,-2), -2.44995181*pow(10,-2), 8.34270064*pow(10,-4), 1.88226433*pow(10,-3), 2.43073327*pow(10,-4)};
//For GEn
vector<double> GEn_Ye_coef = {0.048919981379, -0.064525053912, -0.240825897382, 0.392108744873, 0.300445258602, -0.661888687179, -0.175639769687, 0.624691724461, -0.077684299367, -0.236003975259, 0.090401973470};
vector<double> GEn_Ye_err_coef = {-2.07194073, 1.13809127, 1.01431277, -3.13301380*pow(10,-1), -2.73293676*pow(10,-1), 2.57350595*pow(10,-1), -2.06042113*pow(10,-1), -1.68497332*pow(10,-1), 1.37784515*pow(10,-1), 7.57591964*pow(10,-2), -2.67511301*pow(10,-2), -1.72573088*pow(10,-2), 7.03581500*pow(10,-4), 1.47962095*pow(10,-3), 1.97375221*pow(10,-4)};
//For GMn/mu_n
vector<double> GMn_Div_mun_Ye_coef = {0.257758326959, -1.079540642058, 1.182183812195, 0.711015085833, -1.348080936796, -1.662444025208, 2.624354426029, 1.751234494568, -4.922300878888, 3.197892727312, -0.712072389946};
vector<double> GMn_Div_mun_Ye_err_coef = {-2.06920873, 6.43156400*pow(10,-2), -3.55593786*pow(10,-1), 4.14897660*pow(10,-1), 1.95746824, 2.70525700*pow(10,-1), -1.52685784, -4.43527359*pow(10,-1), 5.16884065*pow(10,-1), 2.07915837*pow(10,-1), -7.48665703*pow(10,-2), -4.25411431*pow(10,-2), 1.54965016*pow(10,-3), 3.25322279*pow(10,-3), .:20819518*pow(10,-4)};

double ye_tcut = 4*physcics_constants::M_pi*physcics_constants::M_pi;
double ye_tnot = -0.7; //GeV^2

//Helper function that calculates form factors based on the parameterization input and stores information in class variables. To be called by the constructor
void calculateFFs_witherr();

//Private function that calculates Q^2. This should be the central value of our data/simulation distribution
double calculate_Q2(double M);

//Private function that calculates the energy of the scattered electron
double calculate_Eprime(double M);

//Private function that calculates tau
double calculate_tau(double M);

//Private function that calculates epsilon
double calculate_epsilon(double tau);

//Private function that calculates the form factor according to Kelly parameterization and 4 input parameters. Handles GEp, GMp/mu_p, GMn/mu_n
double calcKelly(double a1, double b1, double b2, double b3);

//Private function that calculates GEp according to Kelly parameterization
double GEp_Kelly();

//Private function that calculates GMp/mu_p according to Kelly parameterization
double GMp_Kelly();

//Private function that calculates GMn/mu_n according to Kelly parameterization
double GMn_Kelly();

//Private function that calculates GEn according to Kelly
double GEn_Kelly();

//Private function that calculates GEn according to Riordan
double GEn_Riordan();

//Private function that calculates z value for the Ye z expansion
double calculate_YeZ();

//Private function that calculates x value for the Ye x expansion in error analysis
double calculate_YeX();

//Private function that  calculates Ye form factor from the coefs and the ye z.
// G = GEp, GMp/mu_p, GEn, or GMn/mu_n
double calculate_Ye(vector<double>& a_coef, double z);

//Private function that  calculates Ye form factor error  from the coefs and the ye x.
// G = GEp, GMp/mu_p, GEn, or GMn/mu_n
double calculate_Ye_err(vector<double>& b_coef, double x);

//Private function that calculates GEp according to Ye parameterization
double GEp_Ye();

//Private function that calculates GMp/mu_p according to Ye parameterization
double GMp_Ye();

//Private function that calculates GMn/mu_n according to Ye parameterization
double GMn_Ye();

//Private function that calculates GEn according to Ye parameterization
double GEn_Ye();

//Private function that calculates GEp error according to Ye parameterization
double GEp_Ye_err();

//Private function that calculates GMp/mu_p error according to Ye parameterization
double GMp_Ye_err();

//Private function that calculates GMn/mu_n error according to Ye parameterization
double GMn_Ye_err();

//Private function that calculates GEn error according to Ye parameterization
double GEn_Ye_err();

//Private function that calculates the reduced cross-section
double calculate_reduced_cross_section(double GM, double GE, double tau, double epsilon);

//Private function that calculates the error on the reduced cross-section
double calculate_reduced_cross_section_err(double GM, double GM_err, double GE, double GE_err, double tau, double epsilon);

public:
//constructor
calc_FFs_RCS_obj(TString my_param_type,TString kin_file, TString daKin);

//function to calculate the dipole form factor
double GD(double Q2);

//function to calculate the Galster Parameterization. Normally used for just GEn
double Galster(double tau, double Q2, double A, double B);


};//end class
#endif
