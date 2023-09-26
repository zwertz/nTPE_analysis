#ifndef ANALYSIS_UTILITY_FUNCTIONS.H 

#define ANALYSIS_UTILITY_FUNCTIONS.H 

#include <iostream>
#include <fstream>
#include <sstream>
#include "TString.h"
#include "TMath.h"
#include <vector>

//////////////////////////////
//////Static Detector Parameters
//Trigger TDC
const Int_t maxTDCTrigChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
const Double_t tdiffwidecut = 50; // Set at 50ns from nominal 510ns (passed by user) from GMn
//HCal - Note that actual measurement of vertical length is 381.6cm, indicating that the MC figures are correct
const Int_t maxHCalChan = 288; // Total HCal channels
const Int_t maxHCalRows = 24; // Total HCal rows
const Int_t maxHCalCols = 12; // Total HCal cols
const Int_t maxClusters = 10; // Total HCal clusters with information saved
const Double_t HCALHeight = -0.2897; // Height of HCal above beamline in m
const Double_t HCalblk_l = 0.15; // Width and height of HCAL blocks in m
const Double_t HCalblk_l_h_MC = 0.15494; // Horizontal length of all HCAL blocks in m from MC database
const Double_t HCalblk_l_v_MC = 0.15875; // Vertical length of all HCAL blocks in m from MC database
const Double_t posHCalXi = -2.165; // Distance from beam center to top of HCal in m
const Double_t posHCalXf = 1.435; // Distance from beam center to bottom of HCal in m
const Double_t posHCalYi = -0.9; // Distance from beam center to opposite-beam side of HCal in m
const Double_t posHCalYf = 0.9; // Distance from beam center to beam side of HCal in m
const Double_t posHCalXi_MC = -2.3531; // Distance from beam center to top of HCal in m from MC database
const Double_t posHCalXf_MC = 1.45309; // Distance from beam center to bottom of HCal in m from MC database
const Double_t posHCalYi_MC = -0.93155; // Distance from beam center to opposite-beam side of HCal in m from MC database
const Double_t posHCalYf_MC = 0.93155; // Distance from beam center to beam side of HCal in m from MC database
const Double_t HCalSampFrac = 0.077;  //Re-evaluated with MC GEn settings using second to outermost shower column for kin2

//BBCal
const Int_t maxBBCalShChan = 189; // Total BBCal Shower Channels
const Int_t maxBBCalShRows = 27;
const Int_t maxBBCalShCols = 7;
const Int_t maxBBCalPSChan = 52; // Total BBCal Preshower Channels
const Int_t maxBBCalPSRows = 26;
const Int_t maxBBCalPSCols = 2;

//Beamline
const Int_t chargeConvert = 3318; // See D.Flay Doc DB sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf p.8
const Int_t clockActual = 103700; // Needed to convert the 104kHz clock to the actual counting rate

//SBS Magnet
const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
const Double_t maxSBSfield = 1.26; //Tesla
const Double_t SBSfield = 1.0; //fraction of max field. TODO: should be variable per run
const Double_t SBSdist = 2.25; //m
const Double_t dipGap = 1.22; //m
const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;

//GEMs
const Double_t GEMpitch = 10*TMath::DegToRad();

///////////////
/////Physics/Math
const Double_t PI = TMath::Pi();
const Double_t M_e = 0.00051;
const Double_t M_p = 0.938272;
const Double_t M_n = 0.939565;
const UInt_t us = 1000000; //For conversion to seconds used by reporting time delays
const Int_t proton = 2212; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
const Int_t neutron = 2112; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
const Int_t electron = 11; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
const Int_t photon = 22; //pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf

////////////////////////////
//////Static Target/Scattering Chamber Parameters
const Double_t l_tgt = 0.15; // Length of the target (m), keep in m and convert in equation
const Double_t rho_tgt = 0.0723; // Density of target (g/cc)
const Double_t rho_Al = 2.7; // Density of aluminum windows (g/cc)
const Double_t celldiameter = 1.6*2.54; //cm, this is to properly cancel units
const Double_t Ztgt = 1.0;
const Double_t Atgt = 1.0;
const Double_t Mmol_tgt = 1.008; //g/mol
const Double_t dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const Double_t dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const Double_t uwallthick_LH2 = 0.0145; //cm
const Double_t dwallthick_LH2 = 0.015; //cm
const Double_t cellthick_LH2 = 0.02; //cm, this is a guess;
const Double_t Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch 



//This belongs to the header file, not the class. So that way an output file can come from multiple data objects

 string out_dir_temp = "yields_output";
 TString output_directory = TString(out_dir_temp);
 string output_temp = "";
 TString output_file = TString(output_temp);

//Generically implement a couple types of fits
//4th-order polynomial
double poly4_fit(double *x, double *param){

double yint = param[0]; 
double p1 = param[1];
double p2 = param[2];
double p3 = param[3];
double p4 = param[4];

double func = p4*(pow(x[0],4))+p3*(pow(x[0],3))+ p2*(pow(x[0],2))+p1*(x[0])+yint;
return func;
}
//5th-order polynomial
double poly5_fit(double *x, double *param){

double yint = param[0];
double p1 = param[1];
double p2 = param[2];
double p3 = param[3];
double p4 = param[4];
double p5 = param[5];

double func = p5*(pow(x[0],5))+p4*(pow(x[0],4))+p3*(pow(x[0],3))+ p2*(pow(x[0],2))+p1*(x[0])+yint;
return func;
}
//6th-order polynomial
double poly6_fit(double *x, double *param){

double yint = param[0];
double p1 = param[1];
double p2 = param[2];
double p3 = param[3];
double p4 = param[4];
double p5 = param[5];
double p6 = param[6];

double func =p6*(pow(x[0],6))+p5*(pow(x[0],5))+p4*(pow(x[0],4))+p3*(pow(x[0],3))+ p2*(pow(x[0],2))+p1*(x[0])+yint;
return func;
}
//Gaussian
double gaussian_fit(double *x, double *param){

double amp = param[0];
double offset = param[1];
double sigma = param[2];

double func = amp*(exp(-0.5*pow((x[0]-offset)/sigma,2)));
return func;
}
//Exponential 
double expo_fit(double *x, double *param){

double amp = param[0];
double offset = param[1];
double str = param[2];

double func = amp*exp(offset+str*x[0]);
return func;

}


//Function to intialize fit parameters. Currently only supports SBS-4, SBS8, SBS9
vector<double> fit_Params(TString myKin,int sbs_field,TString targ){
vector<double> param (13);
//cout << (myKin == "SBS4") << " " << (sbs_field == 50) << " " << (targ == "LD2") << endl; 
//cout << "kin: " << myKin << " SBS_field:  " << sbs_field << " targ:  " << targ << endl;
if(myKin == "SBS4" && sbs_field == 30 && targ == "LD2"){
param[0] = 1086.83; //used for background p0
param[1] = -254.260;//used for background p1
param[2] = -291.156;//used for background p1
param[3] = 52.3022 ;//used for background p3
param[4] = 28.4984; //used for background p4
param[5] = 12943.1; //used for proton amplitude
param[6] = -0.645471; //used for proton mean
param[7] = 0.181710; //used for proton sigma
param[8] = 5111.02; //used for neutron amplitude
param[9] = 0.000555295; //used for neutron mean
param[10] = 0.174124; //used for neutron sigma
param[11] = -3.0; // min value for fit
param[12] = 2.0; // max value for fit
}
else if(myKin == "SBS4" && sbs_field == 50 && targ == "LD2"){
param[0] = 122.570; //used for background p0
param[1] = -44.2312;//used for background p1
param[2] = -38.1333;//used for background p2
param[3] = 9.83748 ;//used for background p3
param[4] = 4.55384; //used for background p4
param[5] = 1041.39; //used for proton amplitude
param[6] = -1.08320; //used for proton mean
param[7] = 0.184155; //used for proton sigma
param[8] = 429.021; //used for neutron amplitude
param[9] = 0.015794; //used for neutron mean
param[10] = 0.157336; //used for neutron sigma
param[11] = -3.0; // min value for fit
param[12] = 2.0; // max value for fit
}
else if(myKin == "SBS8" && sbs_field == 70 && targ == "LD2"){
param[0] = 2380.80; //used for background
param[1] = -688.796;//used for background
param[2] = -468.258;//used for background
param[3] = 103.072 ;//used for background
param[4] = 44.8213; //used for background
param[5] = 38169.9; //used for proton
param[6] = -0.769101; //used for proton
param[7] = 0.183487; //used for proton
param[8] = 13360.4; //used for neutron
param[9] = 0.0816066; //used for neutron
param[10] = 0.155811; //used for neutron
param[11] = -3.0; // min value for fit
param[12] = 2.0; // max value for fit
}
else if(myKin == "SBS8" && sbs_field == 100 && targ == "LD2"){
param[0] = 290.667; //used for background
param[1] = -93.4792;//used for background
param[2] = -55.806;//used for background
param[3] = 13.9312 ;//used for background
param[4] = 5.26703; //used for background
param[5] = 4274.15; //used for proton
param[6] = -1.12882; //used for proton
param[7] = 0.193342; //used for proton
param[8] = 1636.95; //used for neutron
param[9] = 0.0910029; //used for neutron
param[10] = 0.146692; //used for neutron
param[11] = -3.0; // min value for fit
param[12] = 2.0; // max value for fit
}
else if(myKin == "SBS8" && sbs_field == 50 && targ == "LD2"){
param[0] = 304.939; //used for background
param[1] = -73.2219;//used for background
param[2] = -54.0288;//used for background
param[3] = 9.19629 ;//used for background
param[4] = 4.34144; //used for background
param[5] = 5338.41; //used for proton
param[6] = -0.525136; //used for proton
param[7] = 0.173936; //used for proton
param[8] = 1895.43; //used for neutron
param[9] = 0.086927; //used for neutron
param[10] = 0.149202; //used for neutron
param[11] = -3.0; // min value for fit
param[12] = 2.0; // max value for fit
}

else if(myKin == "SBS9" && sbs_field == 70 && targ == "LD2"){
param[0] = 2724.99; //used for background
param[1] = -749.166;//used for background
param[2] = -680.79;//used for background
param[3] = 111.373 ;//used for background
param[4] = 57.697; //used for background
param[5] = 19177.4; //used for proton
param[6] = -0.82072; //used for proton
param[7] = 0.161423; //used for proton
param[8] = 6996.65; //used for neutron
param[9] = 0.0844196; //used for neutron
param[10] = 0.144542; //used for neutron
param[11] = -3.0; // min value for fit
param[12] = 2.0; // max value for fit
}

else if(myKin == "SBS4" && sbs_field == 0 && targ == "LH2"){
param[0] = 0.0; //used for background p0
param[1] = 1.0;//used for background p1
param[2] = 1.0;//used for background p2
param[3] = 1.0 ;//used for background p3
param[4] = 1.0; //used for background p4
param[5] = 15000; //used for proton amplitude
param[6] = 0.01886; //used for proton mean
param[7] = 0.0517; //used for proton sigma
}


else{
//cout << "kin: " << myKin << " SBS_field:  " << sbs_field << " targ:  " << targ << endl;
//Error message
cout << "Error: The kinematic setting you are analyzing does not have preset fit parameters. Plots will probably not make sense!";
}

return param;
}



//convert an int to a TString. Just a helper function
TString intToTString(int datInt){
stringstream ss;
ss << datInt;
std::string datInt_temp;
ss >> datInt_temp;
TString datInt_string(datInt_temp.c_str());
//cout << datInt_string << endl;
return datInt_string;
}
//convert from degrees to radians
double DegToRad(double myNum){
return (myNum * (TMath::Pi()/180.0));
}
double RadToDeg(double datRad){
return ((datRad * 180.0)/(TMath::Pi()));
}

TString getOutputDir(){
return output_directory;
}
TString makeOutputFileName(TString exp, TString Kin, int SBS_field,TString target){
TString outfile = Form("%s/yields_Zeke_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),Kin.Data(),target.Data(),SBS_field);
//cout << outfile << endl;
return outfile;
}

TString makeReportFileName(TString exp, TString Kin, int SBS_field,TString target){
TString outfile = Form("%s/efficiencyRep_Zeke_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(),Kin.Data(),target.Data(),SBS_field);
//cout << outfile << endl;
return outfile;
}

TString makeYieldReportFileName(TString exp, TString Kin, int SBS_field,TString target){
TString outfile = Form("%s/yieldRep_Zeke_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(),Kin.Data(),target.Data(),SBS_field);
//cout << outfile << endl;
return outfile;
}



//Might be rewriting some code here. But make a class which can then be used to parse the kinematic file info
class kinematic_obj{
//Define some private variables we will fill later.
private:
TString kinematic;
double Ebeam,bbtheta,bbdist,sbstheta,sbsdist,hcaltheta,hcaldist,Q2;
public:
//Let's make a constructor which will just directly parse the kinematic file and store the info 
kinematic_obj(const char *kinematic_file_name,TString Kin){
	kinematic = Kin;
	//Now read-in the kinematic information
	ifstream kinfile(kinematic_file_name);
	//check if there is a problem opening the file
	if(kinfile.fail()){
	cout << "There was a problem with the kinematic file " << kinematic_file_name << ". Figure it out nerd!" << endl;
	return;
	}
	TString datLine;
	TString datKin;
	bool gotKin = false;
	while(datLine.ReadLine(kinfile)){
		if(datLine.BeginsWith(kinematic)){
		//We found the right kinematic. Store the info
		TObjArray *myobjs = datLine.Tokenize(" ");
		//Assuming ordering is kinematic, Beam Energy, BB Angle, SBS angle, SBS dist, HCal dist, Expected Q^2
 		datKin = ((TObjString*) (*myobjs)[0])->GetString();
                Ebeam = (((TObjString*) (*myobjs)[1])->GetString()).Atof();
                bbtheta = (((TObjString*) (*myobjs)[2])->GetString()).Atof();
                bbdist = (((TObjString*) (*myobjs)[3])->GetString()).Atof();
                sbstheta = (((TObjString*) (*myobjs)[4])->GetString()).Atof();
                sbsdist = (((TObjString*) (*myobjs)[5])->GetString()).Atof();
                hcaltheta = (((TObjString*) (*myobjs)[6])->GetString()).Atof();
                hcaldist = (((TObjString*) (*myobjs)[7])->GetString()).Atof();
		Q2 = (((TObjString*) (*myobjs)[8])->GetString()).Atof();
                //cout << "Kinematic " << datKin << "  Beam Energy " << Ebeam << " BB Angle  " << bbtheta << " BB Dist " << bbdist << " SBS Angle  " << sbstheta << " SBS Dist  " << sbsdist << " HCal Angle " << hcaltheta <<  " HCal Dist  " << hcaldist << " Q^2 " << Q2  << endl;
		gotKin = true;
                //cout << gotRun << endl;
 		}
                else{
                //Where are still searching or its a comment
                //cout << "Cond 3" << endl; 
                 continue;
                }
        }
        if ((kinfile.eof()) && !gotKin){
         //Conditional that we checked the entire kinematic file and did not find the kinematic info
	 cout << "Did not find kinematic: " << datKin << " in the kinematic file! Quitting, figure it out!" << endl;
        return;
        }
 }
 TString getKinematic(){
        return kinematic;
 }
 double getBeamEnergy(){
        return Ebeam;
        }
 double getBBAngle_Deg(){
        return bbtheta;
        }
 double getBBAngle_Rad(){
        return DegToRad(bbtheta);
        }
 double getBBDist(){
        return bbdist;
        }
 double getSBSAngle_Deg(){
        return sbstheta;
        }
 double getSBSAngle_Rad(){
        return DegToRad(sbstheta);
        }
 double getSBSDist(){
        return sbsdist;
        }
 double getHCalAngle_Deg(){
        return hcaltheta;
        }
 double getHCalAngle_Rad(){
        return DegToRad(hcaltheta);
        }
 double getHCalDist(){
        return hcaldist;
        }
 double getQ2(){
        return Q2;
        }
 
 void printKinInfo(){
        cout << "------------------------"                              << endl
             << Form("Kinematic: %s,",(getKinematic()).Data())          << endl
             << Form("Beam Energy: %f,",getBeamEnergy())                << endl
             << Form("BB angle in Degrees: %f,",getBBAngle_Deg())       << endl
             << Form("BB angle in Radians: %f,",getBBAngle_Rad())       << endl
             << Form("BB Distance: %f,",getBBDist())                    << endl
             << Form("SBS angle in Degrees: %f,",getSBSAngle_Deg())     << endl
             << Form("SBS angle in Radians: %f,",getSBSAngle_Rad())     << endl
             << Form("SBS Distance: %f,",getSBSDist())                  << endl
             << Form("HCal angle in Degress: %f,",getHCalAngle_Deg())   << endl
             << Form("HCal angle in Radians: %f,",getHCalAngle_Rad())   << endl
             << Form("HCal Distance: %f,",getHCalDist())                << endl
             << "------------------------"                              << endl;
 }
};



//Make a class to store information about the data
class data_object{
//Define some private variables we will fill later.
private:
TString pass,kinematic,target,input_file;
double Ebeam,bbtheta,bbdist,sbstheta,sbsdist,hcaltheta,hcaldist;
int run,sbs_field;
string input_directory = "/work/halla/sbs/sbs-gmn";

TString makeInputFileName(){
 //All of this was to have a modular input directory. So let's make it
 const char *input_directory_char = input_directory.c_str();
 const char *pass_char = pass.Data();
 const char *kin_char = kinematic.Data();
 const char *tar_char = target.Data();
 TString inputfile = Form("%s/%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_char,pass_char,kin_char,tar_char,run);
 //cout << "File Location " << inputfile << endl;               
 return inputfile;
}

//constructor to search through the files we care about and store the information of use. Will handle some cases of unexpected behavior. Handle on a run by run basis. But one could store these in a vector and manipulate them
public:
 data_object(int runnum,const char *data_file_name,const char *kinematic_file_name,TString Kin, TString SBS_field, TString targ){
 
 //This part of the constructor reads in information about the data file itself
 ifstream datafile(data_file_name);
 //check if there is a problem opening the file
 if(datafile.fail()){
 cout << "There was a problem with the data file " << data_file_name << ". Figure it out nerd!" << endl;
 return;
 }
 TString currentLine;
 bool gotRun = false;
 TString runnum_string = intToTString(runnum);
 while(currentLine.ReadLine(datafile)){
	if(currentLine.BeginsWith(runnum_string)){
	 //We found the right run number. Store the info
	 TObjArray *tokens = currentLine.Tokenize(" ");
	 //Assuming ordering is runnum, pass, SBSKin, Target, sbsfield
	 run =  (((TObjString*) (*tokens)[0])->GetString()).Atoi();
	 pass = ((TObjString*) (*tokens)[1])->GetString();
	 kinematic =  ((TObjString*) (*tokens)[2])->GetString();
	 target =  ((TObjString*) (*tokens)[3])->GetString();
	 sbs_field = (((TObjString*) (*tokens)[4])->GetString()).Atoi();
	 //cout << "Run " << run  << " Pass " << pass << " Kin " << kinematic << " Target " << target << " SBS Field  " << sbs_field << endl;
	if(!(kinematic == Kin)){
	cout << "The run " << run << " has a mismatch in the kinematic, investigate what is going on!" << endl;
	return;
	}
	if(!(sbs_field == SBS_field)){
        cout << "The run "<< run << " has a mismatch in the sbs field, investigate what is going on!" << endl;
        return;
        }
 	if(!(target == targ )){
        cout << "The run " << run << " has a mismatch in the target, investigate what is going on!" << endl;
        return;
        }
	gotRun = true;
	 //cout << gotRun << endl;
	 }
	 else{
	 //Where are still searching but it's a comment denoted by #
	 // cout << "Cond 2" << endl; 
	 continue;
	 }
 }
 if ((datafile.eof()) && !gotRun){
 //Conditional that we checked the entire data file and did not find the runnum
 cout << "Did not find run number: " << runnum << " in the data file! Quitting, figure it out!" << endl;
 return;
 }  
 //Now read-in the kinematic information and store the info
 kinematic_obj datKin(kinematic_file_name, Kin);
 
 Ebeam = datKin.getBeamEnergy();
 bbtheta = datKin.getBBAngle_Deg();
 bbdist = datKin.getBBDist();
 sbstheta = datKin.getSBSAngle_Deg();
 sbsdist = datKin.getSBSDist();
 hcaltheta = datKin.getHCalAngle_Deg();
 hcaldist = datKin.getHCalDist();

 input_file = makeInputFileName();
}

 int getRun(){
	return run;
	}
 TString getPass(){
	return pass;
	}
 TString getKinematic(){
	return kinematic;
	}
 TString getTarget(){
	return target;
	}
 int getSBSField(){
	return sbs_field;
	}

 double getBeamEnergy(){
	return Ebeam;
	}		
 double getBBAngle_Deg(){
	return bbtheta;
	}
 double getBBAngle_Rad(){
	return DegToRad(bbtheta);
	}
 double getBBDist(){
	return bbdist;
	}
 double getSBSAngle_Deg(){
        return sbstheta;
        }
 double getSBSAngle_Rad(){
        return DegToRad(sbstheta);
        }
 double getSBSDist(){
        return sbsdist;
        }
 double getHCalAngle_Deg(){
        return hcaltheta;
        }
 double getHCalAngle_Rad(){
        return DegToRad(hcaltheta);
        }
 double getHCalDist(){
        return hcaldist;
        }
 string getInputDir(){
	return input_directory;
	}
 TString getInputFile(){

	return input_file;
	}
 void printRunInfo(){
        cout << "------------------------"                              << endl
             << Form("Run number: %i,",getRun())                        << endl
             << Form("Kinematic: %s,",(getKinematic()).Data())          << endl
             << Form("Target: %s,", (getTarget()).Data())               << endl
             << Form("SBS Field: %i,",getSBSField())       		<< endl
             << Form("Beam Energy: %f,",getBeamEnergy())                << endl
             << Form("BB angle in Degrees: %f,",getBBAngle_Deg())       << endl
             << Form("BB angle in Radians: %f,",getBBAngle_Rad())       << endl
             << Form("BB Distance: %f,",getBBDist())                    << endl
             << Form("SBS angle in Degrees: %f,",getSBSAngle_Deg())     << endl
             << Form("SBS angle in Radians: %f,",getSBSAngle_Rad())     << endl
             << Form("SBS Distance: %f,",getSBSDist())                  << endl
             << Form("HCal angle in Degress: %f,",getHCalAngle_Deg())   << endl
             << Form("HCal angle in Radians: %f,",getHCalAngle_Rad())   << endl
             << Form("HCal Distance: %f,",getHCalDist())                << endl
             << "------------------------"                              << endl;
 }
};
#endif
