#include "../include/physics.h"


//Author: Ezekiel Wertz
//Companion implementation for the phyiscs functions

namespace physics{

 //Define HCal X-axis based on kinematic information
 TVector3 getHCal_xaxis(){
 //This is by convention
 TVector3 hcal_xaxis(0,-1,0);
 
 return hcal_xaxis;
 }

 //Define HCal Y-axis based on kinematic information
 TVector3 getHCal_yaxis(TVector3 hcal_xaxis, TVector3 hcal_zaxis){
 //Since we define hcal_xaxis and hcal_zaxis we require a perpendicular unit vector
 TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

 return hcal_yaxis;
 }

 //Define HCal Z-axis based on kinematic information
 TVector3 getHCal_zaxis(double hcaltheta){
 //Define the zaxis so it is perpendicular to the face of HCal?
 TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
 
 return hcal_zaxis;
 }
 
 //Define the origin for HCal, part of defining the HCal coordinate system
 TVector3 getHCal_origin(double hcaldist, double hcal_offset, TVector3 hcal_xaxis, TVector3 hcal_zaxis){
 TVector3 hcal_origin = hcaldist *hcal_zaxis +hcal_offset*hcal_xaxis;

 return hcal_origin;
 }

 //Try to calculate the effects of the SBS magnetic field
 double getBdL(int sbs_field){
 double BdL = (exp_constants::maxsbsfield)*(exp_constants::sbsdipolegap)*(sbs_field/100);

 return BdL;
 }

 //Calculate the outgoing energy based on losses in the target
 //Currently does not account for the aluminum target cell, requires modification
 double getEloss_outgoing(double bbtheta,TString target){
 double Eloss_outgoing;
        //LH2 
	if(target == "LH2"){
 	Eloss_outgoing = (exp_constants::celldiameter/2)/sin(bbtheta)*(exp_constants::lh2_rho_tgt)*(exp_constants::lh2_dEdx); //Approximately 1 MeV, one could correct more using the raster position
	//LD2
	}else if(target == "LD2"){
	Eloss_outgoing = (exp_constants::celldiameter/2)/sin(bbtheta)*(exp_constants::ld2_rho_tgt)*(exp_constants::ld2_dEdx);
	}else{
	//give an error statement for and don't calculate the Eloss_outgoing.
	cout << "Warning: Target " << target << " is not handle by energy loss function! Defaulting to zero." << endl;
	Eloss_outgoing = 0;
	}

 return Eloss_outgoing;
 }

 //Calculate energy loss in target using vertex position
 double getEloss(double vz, TString target){
 double Eloss;

 	//LH2
 	if(target == "LH2"){
	Eloss = (vz+(exp_constants::l_tgt/2))*(exp_constants::lh2_rho_tgt)*(exp_constants::lh2_dEdx) + (exp_constants::lh2_uwallthick)*(exp_constants::rho_Al)*(exp_constants::dEdx_Al);
	//LD2
	}else if(target == "LD2"){
	Eloss = (vz+(exp_constants::l_tgt/2))*(exp_constants::ld2_rho_tgt)*(exp_constants::ld2_dEdx) + (exp_constants::ld2_uwallthick)*(exp_constants::rho_Al)*(exp_constants::dEdx_Al);
	}else{
        //give an error statement for and don't calculate the Eloss.
         cout << "Warning: Target " << target << " is not handle by energy loss function! Defaulting to zero." << endl;
	Eloss = 0;
	}
 return Eloss;	
 }

 //Correct beam energy for energy loss in target
 double getEcorr(double Ebeam, double Eloss){
 double Ecorr = Ebeam - Eloss; 
 return Ecorr;
 }
 
 //Make a vector based on vertex information. Right now this is only in vz. Could change later and include vx, and vy. Instead of assuming 0
 TVector3 getVertex(double vz){
 TVector3 vertex( 0, 0, vz );

 return vertex;
 }

 //reconstructed momentum, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
 double getp_recon_corr(double tr_p, double Eloss_outgoing){
 double pcorr = tr_p - Eloss_outgoing;
 return pcorr;
 }

 //four momentum vector for electron beam with correted Energy value
 TLorentzVector getpBeam(double Ecorr){
 TLorentzVector Pbeam(0,0,Ecorr,Ecorr);
 return Pbeam;
 }
  
 //four momentum for scattered electron based on reconstruction
 TLorentzVector getp_eprime(double tr_px, double tr_py, double tr_pz, double tr_p, double pcorr){
 TLorentzVector p_eprime(pcorr*(tr_px/tr_p),(pcorr*(tr_py/tr_p),(pcorr*(tr_pz/tr_p),pcorr);
 return p_eprime;
 }


}//end namepspace
