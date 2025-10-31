#ifndef DATA_OBJECT_H
#define DATA_OBJECT_H

//Author: Ezekiel Wertz
//Class to handle all the information associated with a given data run. This class is essential to the data parser script.
//I have not written setter functions for this class as I have never found a need to.

#include "TString.h"

class data_object{
private:
TString pass,kinematic,target,input_file;
double Ebeam,bbtheta,bbdist,sbstheta,sbsdist,hcaltheta,hcaldist,Q2,electron_p,nucleon_p;
int run,sbs_field;

  //Function that forms the input file name for a given run based on the mass-replay pass of the data file
TString makeInputFileName();

public:
//constructor to search through the files we care about and store the information of use. Will handle some cases of unexpected behavior. Handle on a run by run basis. But one could store these in a vector and manipulate them
data_object(int runnum,const char *data_file_name,const char *kinematic_file_name,TString Kin, TString SBS_field, TString targ, TString daPass);

//destructor
//no dynamically allocated memory or pointers
~data_object();


  //getter functions to get back a single private class variable as needed.
int getRun();

TString getPass();

TString getKinematic();

TString getTarget();

int getSBSField();

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

TString getInputFile();

void printRunInfo();


};//end class
#endif
