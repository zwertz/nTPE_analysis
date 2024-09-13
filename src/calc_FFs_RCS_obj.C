#include "../include/calc_FFs_RCS_obj.h"

//Author: Ezekiel Wertz
//A companion implementation. A class to facilite calculating electromagnetic form factors via various supported parameterization. The most critical being what is used in SIMC. The class will also calculate the reduced cross-section. The FFs and the Reduced Cross-Section will be dependent on user-input which will allow for multiple supposed parameterizations. Initially this will include Kelly, Riordan, and Ye. The class object will hold these all as class variables. Public and private functions and the constructor will enable this class to be compatible with other scripts and classes in this analysis framework.

//Kelly Parameterization Citation
//(Phys. Rev. C 70, 068202)

//Riordan Parameterization Citation
//(Phys. Rev. Lett. 105, 262302)

//Ye Parameterization Citation
//( Physics Letters B 777, 8 (2018))

//Helper function that calculates form factors based on the parameterization input and stores information in class variables. To be called by the constructor
void calc_FFs_RCS_obj::calculateFFs_witherr(){

	//Form factors are parameterization specific
	if(param_type == "Kelly"){
		GEp = GEp_Kelly();
		GMp_Div_mup = GMp_Kelly(); 
		GEn = GEn_Kelly();
		GMn_Div_mun = GMn_Kelly(); 
		//Use the standard physical constants
		GMp = GMp_Div_mup * physics_constants::mu_p;
		GMn = GMn_Div_mun * abs(physics_constants::mu_n);
		//Kelly does have a way to do error evaluation via a covariance matrix. However, this calculation is nontrival and not yet implemented. Defaulting to 0.
		GEp_err = 0;
		GMp_Div_mup_err = 0;
		GEn_err = 0;
		GMn_Div_mun_err = 0;
	}else if(param_type == "SIMC"){
		//Kelly for everything but GEn
		//GEn is according to Riordan
		GEp = GEp_Kelly();
                GMp_Div_mup = GMp_Kelly();
		GMn_Div_mun = GMn_Kelly();
		GEn = GEn_Riordan();
		//Hard-coded values for the proton and neutron magnetic moments so exactly the same as SIMC
		GMp = GMp_Div_mup * 2.79;
		GMn = GMn_Div_mun * abs(-1.913);
		//We should not assign an error value to Rsim and SIMC does not. Defaulting to 0.
		GEp_err = 0;
                GMp_Div_mup_err = 0;
                GEn_err = 0;
                GMn_Div_mun_err = 0;
	}else if(param_type == "Ye"){
		GEp = GEp_Ye(); 
                GMp_Div_mup = GMp_Ye();
                GEn = GEn_Ye();
                GMn_Div_mun = GMn_Ye();
		//Use the standard physical constants
                GMp = GMp_Div_mup * physics_constants::mu_p;
                GMn = GMn_Div_mun * abs(physics_constants::mu_n);
		//Ye has well developed and reproducible error analysis. 
		GEp_err = GEp_Ye_err();
                GMp_Div_mup_err = GMp_Ye_err();
                GEn_err = GEn_Ye_err();
                GMn_Div_mun_err = GMn_Ye_err();
	}else{
		cout << "The parameterization choice: " << param_type << " is not supported by this class. Either fix the issue or implement the parameterization!" << endl;
	}
	//Call the dipole for use here
	double GD = calculate_GD(Q2);
	//Modifications to form factors which can be done independently
	GEp_Div_GD = GEp/GD;
	GMp_Div_mup_GD = GMp_Div_mup/GD;
	GEn_Div_GD = GEn/GD;
	GMn_Div_mun_GD = GMn_Div_mun/GD;
	//Modification to the form factor errors, which can be done independently
	GEp_Div_GD_err = GEp_err/GD;
	GMp_err = GMp_Div_mup_err * physics_constants::mu_p;
	GMp_Div_mup_GD_err = GMp_Div_mup_err/GD;
	GEn_Div_GD_err = GEn_err/GD;
	GMn_err = GMn_Div_mun_err* abs(physics_constants::mu_n);
	GMn_Div_mun_GD_err = GMn_Div_mun_err/GD;
}

//Private function that calculates Q^2. This should be the central value of our data/simulation distribution
double calc_FFs_RCS_obj::calculate_Q2(double M){
double Q2 = 4 * Ebeam*Eprime*pow(sin(BB_angle_rad/2),2);
return Q2;
}


//Private function that calculates the energy of the scattered electron
double calc_FFs_RCS_obj::calculate_Eprime(double M){
double Eprime = Ebeam / (1+(2*Ebeam/M)*pow(sin(BB_angle_rad/2),2));
return Eprime;
}

//Private function that calculates tau
double calc_FFs_RCS_obj::calculate_tau(double M){
double tau = Q2 / (4 * M * M);;
return tau;
}


//Private function that calculates epsilon
double calc_FFs_RCS_obj::calculate_epsilon(double tau){
double epsilon = 1/(1+2*(1+tau)*pow(tan(BB_angle_rad/2),2));
return epsilon;
}

//Private function that calculates the form factor according to Kelly parameterization and 4 input parameters. Handles GEp, GMp/mu_p, GMn/mu_n
double calc_FFs_RCS_obj::calcKelly(double tau, double a1, double b1, double b2, double b3){
double numerator = 1 + a1*tau;
double denominator = 1 + b1*tau + b2*pow(tau,2) + b3*pow(tau,3);
double KellyFF =  numerator/denominator;
return KellyFF;
}

//Private function that calculates GEp according to Kelly parameterization
double calc_FFs_RCS_obj::GEp_Kelly(){
double a1 = Kelly_GEp[0];
double b1 = Kelly_GEp[1];
double b2 = Kelly_GEp[2];
double b3 = Kelly_GEp[3];

double GEp_Kelly = calcKelly(tau_p,a1,b1,b2,b3);
return GEp_Kelly;
}

//Private function that calculates GMp/mu_p according to Kelly parameterization
double calc_FFs_RCS_obj::GMp_Kelly(){
double a1 = Kelly_GMp[0];
double b1 = Kelly_GMp[1];
double b2 = Kelly_GMp[2];
double b3 = Kelly_GMp[3];

double GMp_Kelly = calcKelly(tau_p,a1,b1,b2,b3);
return GMp_Kelly;
}

//Private function that calculates GMn/mu_n according to Kelly parameterization
double calc_FFs_RCS_obj::GMn_Kelly(){
double a1 = Kelly_GMn[0];
double b1 = Kelly_GMn[1];
double b2 = Kelly_GMn[2];
double b3 = Kelly_GMn[3];

double GMn_Kelly = calcKelly(tau_p,a1,b1,b2,b3);
return GMn_Kelly;
}

//Private function that calculates GEn according to Kelly
double calc_FFs_RCS_obj::GEn_Kelly(){
double GEn_Kelly = calculate_Galster(tau_p, Q2, Kelly_A, Kelly_B);
return GEn_Kelly;
}

//Private function that calculates GEn according to Riordan
double calc_FFs_RCS_obj::GEn_Riordan(){
double GEn = calculate_GD(Q2)* (1.520*tau_p + 2.629*pow(tau_p,2) + 3.055*pow(tau_p,3))/(1.0+5.222*tau_p+0.040*pow(tau_p,2)+11.438*pow(tau_p,3));
return GEn;
}

//Private function that calculates z value for the Ye z expansion
double calc_FFs_RCS_obj::calculate_YeZ(){
double numerator = sqrt(ye_tcut + Q2) - sqrt(ye_tcut - ye_tnot);
double denominator = sqrt(ye_tcut + Q2) + sqrt(ye_tcut - ye_tnot);
double z = numerator/denominator;
return z;
}

//Private function that calculates x value for the Ye x expansion in error analysis
double calc_FFs_RCS_obj::calculate_YeX(){
double x;
	if (Q2 == 1){
	x = 0.00000001;
	}else{
      	x = log10(Q2);
	}
return x;
}

//Private function that  calculates Ye form factor from the coefs and the ye z.
// G = GEp, GMp/mu_p, GEn, or GMn/mu_n
double calc_FFs_RCS_obj::calculate_Ye(vector<double>& a_coef, double z){
double G = 0;
	for (int i = 0; i<a_coef.size();i++){
	G += a_coef[i]*pow(z, i);
	}
return G;
}

//Private function that  calculates Ye form factor error  from the coefs and the ye x.
// G = GEp, GMp/mu_p, GEn, or GMn/mu_n
double calc_FFs_RCS_obj::calculate_Ye_err(vector<double>& b_coef, double x){
//This is only implemented for a Q^2 range of 10^-3 <= Q^2 <= 10^2. It does not support the other 2 regions which are uncessary for this analysis
double exponent = 0;
   for (int i = 0; i<b_coef.size();i++){
      exponent += b_coef[i]*pow(x, i);
    }
   double G_err = pow(10,exponent) * calculate_GD(Q2);
   return G_err;
}

//Private function that calculates GEp according to Ye parameterization
double calc_FFs_RCS_obj::GEp_Ye(){
double z = calculate_YeZ();
double GEp_Ye = calculate_Ye(GEp_Ye_coef,z);
return GEp_Ye;
}

//Private function that calculates GMp/mu_p according to Ye parameterization
double calc_FFs_RCS_obj::GMp_Ye(){
double z = calculate_YeZ();
double GMp_Ye = calculate_Ye(GMp_Div_mup_Ye_coef,z);
return GMp_Ye;
}

//Private function that calculates GMn/mu_n according to Ye parameterization
double calc_FFs_RCS_obj::GMn_Ye(){
double z = calculate_YeZ();
double GMn_Ye = calculate_Ye(GMn_Div_mun_Ye_coef,z);
return GMn_Ye;
}

//Private function that calculates GEn according to Ye parameterization
double calc_FFs_RCS_obj::GEn_Ye(){
double z = calculate_YeZ();
double GEn_Ye = calculate_Ye(GEn_Ye_coef,z);
return GEn_Ye;
}

//Private function that calculates GEp error according to Ye parameterization
double calc_FFs_RCS_obj::GEp_Ye_err(){
double x = calculate_YeX();
double GEp_Ye_err = calculate_Ye_err(GEp_Ye_err_coef,x);
return GEp_Ye_err;
}

//Private function that calculates GMp/mu_p error according to Ye parameterization
double calc_FFs_RCS_obj::GMp_Ye_err(){
double x = calculate_YeX();
double GMp_Ye_err = calculate_Ye_err(GMp_Div_mup_Ye_err_coef,x);
return GMp_Ye_err;
}

//Private function that calculates GMn/mu_n error according to Ye parameterization
double calc_FFs_RCS_obj::GMn_Ye_err(){
double x = calculate_YeX();
double GMn_Ye_err = calculate_Ye_err(GMn_Div_mun_Ye_err_coef,x);
return GMn_Ye_err;
}

//Private function that calculates GEn error according to Ye parameterization
double calc_FFs_RCS_obj::GEn_Ye_err(){
double x = calculate_YeX();
double GEn_Ye_err = calculate_Ye_err(GEn_Ye_err_coef,x);
return GEn_Ye_err;
}

//Private function that calculates the reduced cross-section
double calc_FFs_RCS_obj::calculate_reduced_cross_section(double GM, double GE, double tau, double epsilon){
double RCS = GM*GM + (epsilon/tau)*GE*GE;
return RCS;
}

//Private function that calculates the error on the reduced cross-section
double calc_FFs_RCS_obj::calculate_reduced_cross_section_err(double GM, double GM_err, double GE, double GE_err, double tau, double epsilon){
double GM_partial = 2 * GM;
double GE_partial  = 2 * (epsilon/tau)* GE;

double variance = pow(GM_partial,2)* pow(GM_err,2) + pow(GE_partial,2)*pow(GE_err,2);
double stdev  =  sqrt(variance);
return stdev;
}

//Standard constructor implementation
calc_FFs_RCS_obj::calc_FFs_RCS_obj(TString my_param_type, TString kin_file, TString daKin){
//Initialize the logistical info
param_type = my_param_type;
kinematic_file = kin_file;
kin = daKin;

//First make a kinematic object to get and store the information we need about the kinematic in class variables
kinematic_obj myKin(kinematic_file,kin);
Ebeam = myKin.getBeamEnergy();
BB_angle_rad = myKin.getBBAngle_Rad();

//calculate average mass of proton and neutron
double M_avg = (physics_constants::M_p + physics_constants::M_n)/2;
//calculate the scattered electron energy using the average mass
Eprime = calculate_Eprime(M_avg);
//calculate the central Q^2 of the kinematics using the average mass
Q2 = calculate_Q2(M_avg);

//calculate the dimensionless Q2 quantities using proton and neutron masses
tau_p = calculate_tau(physics_constants::M_p);
tau_n = calculate_tau(physics_constants::M_n);

//calculate the polarization fo the virtual photon, epsilon, using proton and nuetron taus
epsilon_p = calculate_epsilon(tau_p);
epsilon_n = calculate_epsilon(tau_n);

//Calculate the form factors and the error according to the input parameterization type. Store them in the class variables
calculateFFs_witherr();

//Calculate the reduced cross-sections according the chosen parameterization and the stored class variables
reduced_cross_section_n = calculate_reduced_cross_section(GMn,GEn,tau_n,epsilon_n);
reduced_cross_section_p = calculate_reduced_cross_section(GMp,GEp,tau_p,epsilon_p);

//Calculate the errors on the reduced cross section
reduced_cross_section_n_err = calculate_reduced_cross_section_err(GMn,GMn_err,GEn,GEn_err,tau_n,epsilon_n);
reduced_cross_section_p_err = calculate_reduced_cross_section_err(GMp,GMp_err,GEp,GEp_err,tau_p,epsilon_p);

//Calculate the ratio of the n/p reduced cross-sections
reduced_cross_section_ratio = reduced_cross_section_n/reduced_cross_section_p;

//Calculate the error of the ratio of the n/p reduced cross-sections
reduced_cross_section_ratio_err = reduced_cross_section_ratio * sqrt(pow(reduced_cross_section_n_err/reduced_cross_section_n,2)+pow(reduced_cross_section_p_err/reduced_cross_section_p,2));

}
//destructor
//Since we don't have any dynamically allocated memory or pointers this should be ok
//This implicitly calls the destructors of the base class variables
calc_FFs_RCS_obj::~calc_FFs_RCS_obj(){
}

//function to calculate the dipole form factor
double calc_FFs_RCS_obj::calculate_GD(double Q2){
double GD = pow(1+Q2/GD_const,-2);
return GD;
}

//function to calculate the Galster Parameterization. Normally used for just GEn
double calc_FFs_RCS_obj::calculate_Galster(double tau, double Q2, double A, double B){
double Gal = ((A*tau)/(1+B*tau))*calculate_GD(Q2);
return Gal;
}

//Public getter functions for most if not all the private class variables

TString calc_FFs_RCS_obj::get_param_type(){return param_type;}

TString calc_FFs_RCS_obj::get_kin(){return kin;}

TString calc_FFs_RCS_obj::get_kinematic_file(){return kinematic_file;}

double calc_FFs_RCS_obj::get_reduced_cross_section_n(){return reduced_cross_section_n;}

double calc_FFs_RCS_obj::get_reduced_cross_section_n_err(){return reduced_cross_section_n_err;}

double calc_FFs_RCS_obj::get_reduced_cross_section_p(){return reduced_cross_section_p;}

double calc_FFs_RCS_obj::get_reduced_cross_section_p_err(){return reduced_cross_section_p_err;}

double calc_FFs_RCS_obj::get_reduced_cross_section_ratio(){return reduced_cross_section_ratio;}

double calc_FFs_RCS_obj::get_reduced_cross_section_ratio_err(){return reduced_cross_section_ratio_err;}

double calc_FFs_RCS_obj::get_Ebeam(){return Ebeam;}

double calc_FFs_RCS_obj::get_Eprime(){return Eprime;}

double calc_FFs_RCS_obj::get_BB_angle_rad(){return BB_angle_rad;}

double calc_FFs_RCS_obj::get_Q2(){return Q2;}

double calc_FFs_RCS_obj::get_tau_p(){return tau_p;}

double calc_FFs_RCS_obj::get_tau_n(){return tau_n;}

double calc_FFs_RCS_obj::get_epsilon_p(){return epsilon_p;}

double calc_FFs_RCS_obj::get_epsilon_n(){return epsilon_n;}

double calc_FFs_RCS_obj::get_GEp(){return GEp;}

double calc_FFs_RCS_obj::get_GEp_Div_GD(){return GEp_Div_GD;}

double calc_FFs_RCS_obj::get_GMp(){return GMp;}

double calc_FFs_RCS_obj::get_GMp_Div_mup(){return GMp_Div_mup;}

double calc_FFs_RCS_obj::get_GMp_Div_mup_GD(){return GMp_Div_mup_GD;}

double calc_FFs_RCS_obj::get_GEn(){return GEn;}

double calc_FFs_RCS_obj::get_GEn_Div_GD(){return GEn_Div_GD;}

double calc_FFs_RCS_obj::get_GMn(){return GMn;}

double calc_FFs_RCS_obj::get_GMn_Div_mun(){return GMn_Div_mun;}

double calc_FFs_RCS_obj::get_GMn_Div_mun_GD(){return GMn_Div_mun_GD;}

double calc_FFs_RCS_obj::get_GEp_err(){return GEp_err;}

double calc_FFs_RCS_obj::get_GEp_Div_GD_err(){return GEp_Div_GD_err;}

double calc_FFs_RCS_obj::get_GMp_err(){return GMp_err;}

double calc_FFs_RCS_obj::get_GMp_Div_mup_err(){return GMp_Div_mup_err;}

double calc_FFs_RCS_obj::get_GMp_Div_mup_GD_err(){return GMp_Div_mup_GD_err;}

double calc_FFs_RCS_obj::get_GEn_err(){return GEn_err;}

double calc_FFs_RCS_obj::get_GEn_Div_GD_err(){return GEn_Div_GD_err;}

double calc_FFs_RCS_obj::get_GMn_err(){return GMn_err;}

double calc_FFs_RCS_obj::get_GMn_Div_mun_err(){return GMn_Div_mun_err;}

double calc_FFs_RCS_obj::get_GMn_Div_mun_GD_err(){return GMn_Div_mun_GD_err;}

//Print all the relevant information. Mainly to compare to other calculations. QA check as well
void calc_FFs_RCS_obj::print(){
cout << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl
     << Form("Parameterization Type: %s",get_param_type().Data())												   << endl
     << Form("Kinematic: %s",get_kin().Data())															   << endl
     << Form("Beam Energy: %f", get_Ebeam())															   << endl
     << Form("Scattered Electron Energy: %f", get_Eprime())													   << endl 
     << Form("BigBite Angle (radians): %f", get_BB_angle_rad())													   << endl
     << Form("Q^{2}: %f",get_Q2())																   << endl
     << Form("Dimensionless Four Momentum Tau (for proton): %f", get_tau_p())											   << endl
     << Form("Dimensionless Four Momentum Tau (for neutron): %f", get_tau_n())                                                                                     << endl 
     << Form("Polarization of Virtual Photon Epsilon (for proton): %f", get_epsilon_p())									   << endl
     << Form("Polarization of Virtual Photon Epsilon (for neutron): %f", get_epsilon_n())                                                                          << endl     
     << Form("GD: %f",calculate_GD(Q2)) 		                                                                                                           << endl
     << Form("GEp: %f +/- %f",get_GEp(),get_GEp_err())														   << endl
     << Form("GEp/GD: %f +/- %f",get_GEp_Div_GD(),get_GEp_Div_GD_err())                                                                                            << endl 
     << Form("GMp: %f +/- %f",get_GMp(),get_GMp_err())                                                                                                             << endl
     << Form("GMp/mu_p: %f +/- %f",get_GMp_Div_mup(),get_GMp_Div_mup_err())											   << endl
     << Form("GMp/mu_pGD: %f +/- %f",get_GMp_Div_mup_GD(),get_GMp_Div_mup_GD_err())                                                                                << endl
     << Form("GEn: %f +/- %f",get_GEn(),get_GEn_err())                                                                                                             << endl
     << Form("GEn/GD: %f +/- %f",get_GEn_Div_GD(),get_GEn_Div_GD_err())                                                                                            << endl
     << Form("GMn: %f +/- %f",get_GMn(),get_GMn_err())                                                                                                             << endl
     << Form("GMn/mu_n: %f +/- %f",get_GMn_Div_mun(),get_GMn_Div_mun_err())                                                                                        << endl
     << Form("GMn/mu_nGD: %f +/- %f",get_GMn_Div_mun_GD(),get_GMn_Div_mun_GD_err())                                                                                << endl    
     << Form("MC Neutron reduced cross-section: %f +/- %f",get_reduced_cross_section_n(),get_reduced_cross_section_n_err())					   << endl
     << Form("MC Proton reduced cross-section: %f +/- %f",get_reduced_cross_section_p(),get_reduced_cross_section_p_err())                                         << endl
     << Form("MC n/p reduced cross-section ratio %f +/- %f",get_reduced_cross_section_ratio(), get_reduced_cross_section_ratio_err())				   << endl
     << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl;
}

//Print but for a report file
void calc_FFs_RCS_obj::print(std::ofstream& report){
report << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl
     << Form("Parameterization Type: %s",get_param_type().Data())                                                                                                  << endl
     << Form("Kinematic: %s",get_kin().Data())                                                                                                                     << endl
     << Form("Beam Energy: %f", get_Ebeam())                                                                                                                       << endl
     << Form("Scattered Electron Energy: %f", get_Eprime())                                                                                                        << endl
     << Form("BigBite Angle (radians): %f", get_BB_angle_rad())                                                                                                    << endl
     << Form("Q^{2}: %f",get_Q2())                                                                                                                                 << endl
     << Form("Dimensionless Four Momentum Tau (for proton): %f", get_tau_p())                                                                                      << endl
     << Form("Dimensionless Four Momentum Tau (for neutron): %f", get_tau_n())                                                                                     << endl
     << Form("Polarization of Virtual Photon Epsilon (for proton): %f", get_epsilon_p())                                                                           << endl
     << Form("Polarization of Virtual Photon Epsilon (for neutron): %f", get_epsilon_n())                                                                          << endl
     << Form("GD: %f",calculate_GD(Q2))                                                                                                                            << endl
     << Form("GEp: %f +/- %f",get_GEp(),get_GEp_err())                                                                                                             << endl
     << Form("GEp/GD: %f +/- %f",get_GEp_Div_GD(),get_GEp_Div_GD_err())                                                                                            << endl
     << Form("GMp: %f +/- %f",get_GMp(),get_GMp_err())                                                                                                             << endl
     << Form("GMp/mu_p: %f +/- %f",get_GMp_Div_mup(),get_GMp_Div_mup_err())                                                                                        << endl
     << Form("GMp/mu_pGD: %f +/- %f",get_GMp_Div_mup_GD(),get_GMp_Div_mup_GD_err())                                                                                << endl
     << Form("GEn: %f +/- %f",get_GEn(),get_GEn_err())                                                                                                             << endl
     << Form("GEn/GD: %f +/- %f",get_GEn_Div_GD(),get_GEn_Div_GD_err())                                                                                            << endl
     << Form("GMn: %f +/- %f",get_GMn(),get_GMn_err())                                                                                                             << endl
     << Form("GMn/mu_n: %f +/- %f",get_GMn_Div_mun(),get_GMn_Div_mun_err())                                                                                        << endl
     << Form("GMn/mu_nGD: %f +/- %f",get_GMn_Div_mun_GD(),get_GMn_Div_mun_GD_err())                                                                                << endl
     << Form("MC Neutron reduced cross-section: %f +/- %f",get_reduced_cross_section_n(),get_reduced_cross_section_n_err())                                        << endl
     << Form("MC Proton reduced cross-section: %f +/- %f",get_reduced_cross_section_p(),get_reduced_cross_section_p_err())                                         << endl
     << Form("MC n/p reduced cross-section ratio %f +/- %f",get_reduced_cross_section_ratio(), get_reduced_cross_section_ratio_err())                              << endl
     << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl;
}

