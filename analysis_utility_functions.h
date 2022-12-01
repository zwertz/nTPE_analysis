#ifndef ANALYSIS_UTILITY_FUNCTIONS.H
#define ANALYSIS_UTILITY_FUNCTIONS.H

#include <iostream>
#include <fstream>
#include <sstream>

//This belongs to the header file, not the class. So that way an output file can come from multiple data objects
private string output_directory = "yields_output";
private TString output_file= "";

//Make a class to store information about the data
class data_object{
//Define some private variables we will fill later.
private TString pass,kinematic,target,input_file;
private double Ebeam,bbtheta,bbdist,sbstheta,sbsdist,hcaltheta,hcaldist;
int run,sbs_field;
private string input_directory = "/work/halla/sbs/sbs-gmn";
//constructor to search through the files we care about and store the information of use. Will handle some cases of unexpected behavior. Handle on a run by run basis. But one could store these in a vector and manipulate them
public data_object(int runnum,const char *data_file_name,const char *kinematic_file_name,TString Kin, TString SBS_field){
 
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
	 //cout << "Run " << datRun << " Pass " << pass << " Kin " << kinematic << " Target " << target << " SBS Field  " << sbs_field << endl;
	if(!(kinematic == Kin)){
	cout << "The run has a mismatch in the kinematic, investigate what is going on!" << endl
	return;
	}
	if(!(sbs_field == SBS_field)){
        cout << "The run has a mismatch in the sbs field, investigate what is going on!" << endl
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
 //Now read-in the kinematic information
 ifstream kinfile(kinematic_file_name);
 //check if there is a problem opening the file
 if(datafile.fail()){
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
 	//Assuming ordering is kinematic, Beam Energy, BB Angle, SBS angle, SBS dist, HCal dist
 	datKin = ((TObjString*) (*myobjs)[0])->GetString();
 	Ebeam = (((TObjString*) (*myobjs)[1])->GetString()).Atof();
 	bbtheta = (((TObjString*) (*myobjs)[2])->GetString()).Atof();
 	bbdist = (((TObjString*) (*myobjs)[3])->GetString()).Atof();
 	sbstheta = (((TObjString*) (*myobjs)[4])->GetString()).Atof();
 	sbsdist = (((TObjString*) (*myobjs)[5])->GetString()).Atof();
 	hcaltheta = (((TObjString*) (*myobjs)[6])->GetString()).Atof();
 	hcaldist = (((TObjString*) (*myobjs)[7])->GetString()).Atof();
 	//cout << "Kinematic " << datKin << "  Beam Energy " << Ebeam << " BB Angle  " << bbtheta << " BB Dist " << bbdist << " SBS Angle  " << sbstheta << " SBS Dist  " << sbsdist << " HCal Angle " << hcaltheta <<  " HCal Dist  " << hcaldist  << endl;
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
 input_file = makeInputFileName();


}


public int getRun(){
	return run;
	}
public TString getPass(){
	return pass;
	}
public TString getKinematic(){
	return kinematic;
	}
public TString getTarget(){
	return target;
	}
public int getSBSField(){
	return sbs_field;
	}

public double getBeamEnergy(){
	return Ebeam;
	}		
public double getBBAngle_Deg(){
	return bbtheta;
	}
public double getBBAngle_Rad(){
	return DegToRad(bbtheta);
	}
public double getBBDist(){
	return bbdist;
	}
public double getSBSAngle_Deg(){
        return sbstheta;
        }
public double getSBSAngle_Rad(){
        return DegToRad(sbstheta);
        }
public double getSBSDist(){
        return sbsdist;
        }
public double getHCalAngle_Deg(){
        return hcaltheta;
        }
public double getHCalAngle_Rad(){
        return DegToRad(hcaltheta);
        }
public double getHCalDist(){
        return hcaldist;
        }
public string getInputDir(){
	return input_directory;
	}
public TString getInputFile(){

	return input_file;
	}
private TString makeInputFileName(){
 //All of this was to have a modular input directory. So let's make it
 const char *input_directory_char = input_directory.c_str();
 const char *datRun_char = datRun.Data();
 const char *pass_char = pass.Data();
 const char *kin_char = kinematic.Data();
 const char *tar_char = target.Data();
 TString inputfile = Form("%s/%s/%s/%s/rootfiles/e1209019_fullreplay_%s_*.root",input_directory_char,pass_char,kin_char,tar_char,datRun_char);
 //cout << "File Location " << inputfile << endl;             	
	return inputfile;
	}
public void printRunInfo(){
	cout << "------------------------"	 			<< endl
	     << Form("Run number: %i,",getRun()) 			<< endl
	     << Form("Kinematic: %s,",getKinematic())			<< endl
	     << Form("Target: %s,", getTarget())			<< endl
	     << Form("SBS Field: %i,",getSBSField())			<< endl
	     << Form("Beam Energy: %d,",getBeamEnergy())		<< endl
	     << Form("BB angle in Degrees: %d,",getBBAngle_Deg())	<< endl 
	     << Form("BB angle in Radians: %d,",getBBAngle_Rad())	<< endl
	     << Form("BB Distance: %d,",getBBDist())		 	<< endl
	     << Form("SBS angle in Degrees: %d,",getSBSAngle_Deg())	<< endl
	     << Form("SBS angle in Radians: %d,",getSBSAngle_Rad())	<< endl
	     << Form("SBS Distance: %d,",getSBSDist())			<< endl
	     << Form("HCal angle in Degress: %d,",getHCalAngle_Deg())	<< endl
	     << Form("HCal angle in Radians: %d,",getHCalAngle_Rad())	<< endl
	     << Form("HCal Distance: %d,",getHCalDist())		<< endl
 	     << "------------------------"                              << endl;
	}

}
//convert an int to a TString. Just a helper function
public TString intToTString(int datInt){
stringstream ss;
 ss << datInt;
 std::string datInt_temp;
 ss >> datInt_temp;
 TString datInt_string(datInt_temp.c_str());
 //cout << datInt_string << endl;
 return datInt_string;
}

//convert from degrees to radians
public double DegToRad(myNum){
return (myNum * (TMath::Pi()/180.0));
}



public string getOutputDir(){
        return output_directory;
        }

public TString makeOutputFileName(TString exp, TString Kin, TString SBS_field,TString target){
	TString outfile = Form("%s/yields_Zeke_%s_%s_%s_%s.root",getOutputDir(),exp,Kin,target,SBS_field);

        return outfile;
        }

#endif
