#ifndef KINEMATIC_OBJ_H
#define KINEMATIC_OBJ_H

//Author: Ezekiel Wertz
//A class that handles all the information necessary to describe an SBS kinematic for GMN/nTPE. Sort of a companion class to the data_obj.
//Don't see setter funcitons as necessary

class kinematic_obj{
private:
TString kinematic;
double Ebeam,bbtheta,bbdist,sbstheta,sbsdist,hcaltheta,hcaldist,Q2,electron_p,nucleon_p;

public:
//Constructor
kinematic_obj(const char *kinematic_file_name,TString Kin);

//desctructor
//no dyanmic memory and no points
~kinematic_obj();

//Necessary getter functions
TString getKinematic();

double getBeamEnergy();

double getBBAngle_Deg();

double getBBAngle_Rad();

double getBBDist();

double getSBSAngle_Deg();

double getSBSAngle_Rad();

double getSBSDist();

double getHCalAngle_Deg();

double getHCalAngle_Rad();

double getHCalDist();

double getQ2();

double getElectronP();

double getNucleonP();

void printKinInfo();

}; //end class

#endif
