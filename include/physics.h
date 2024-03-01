#ifndef PHYSICS_H
#define PHYSICS_H

#include "TLorentzVector.h"
#include "TVector3.h"
//Author Ezekiel Wertz
//A namespace to hold functions that calculate or update relevant physics analysis variables or quantities


namespace physics{

 //Define HCal X-axis based on kinematic information
 TVector3 getHCal_xaxis();

 //Define HCal Y-axis based on kinematic information
 TVector3 getHCal_yaxis(TVector3 hcal_xaxis, TVector3 hcal_yaxis);

 //Define HCal Z-axis based on kinematic information
 TVector3 getHCal_zaxis(double hcaltheta);

 //Define the origin for HCal, part of defining the HCal coordinate system
 TVector3 getHCal_origin(double hcaldist, double hcal_offset, TVector3 hcal_xaxis, TVector3 hcal_zxaxis);

 //Try to calculate the effects of the SBS magnetic field
 double getBdL(int sbs_field);

 //Calculate the outgoing energy based on losses in the target
 //Currently does not account for the aluminum target cell, requires modification
 double getEloss_outgoing(double bbtheta,TString target);

 //Calculate the energy loss in target using vertex position
 double getEloss(double vz, TString target);

 //Correct beam energy for energy loss in target
 double getEcorr(double Ebeam, double Eloss);
 
 //Make a vector based on vertex information. Right now this is only in vz. Could change later and include vx, and vy. Instead of assuming 0
 TVector3 getVertex(double vz); 

 //reconstructed momentum, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
 double getp_recon_corr(double tr_p, double Eloss_outgoing); 
 
 //four momentum vector for electron beam with correted Energy value
 TLorentzVector getpBeam(double Ecorr);

}//end namespace

#endif
