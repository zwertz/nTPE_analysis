//Author Ezekiel Wertz
//03/15/2024
//Purpose: Parsing simulation information from a simc generator by nucleon, build histrograms for further analysis

#include "../../src/utility.C"
#include "../../include/physics_constants.h"
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/parse_config.C"
#include "../../src/data_object.C"
#include "../../src/cuts.C"
#include "../../src/physics.C"
#include "../../src/plots.C"

//Main
void mc_parse(const char *setup_file_name){

//Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//parse object to get in the information that The One Config file has and is manipulated
parse_config mainConfig(setup_file_name);
//mainConfig.printMCYields();

//store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
TString rootfile_dir = mainConfig.getRootFileDir();
TString histfile_dir = mainConfig.getHistFileDir();
TString replay_type = mainConfig.getReplayType();
TCut globalcut = mainConfig.getGlobalCut();
TString exp = mainConfig.getExp();
TString kin = mainConfig.getKin();
TString kinematic_file = mainConfig.getKinFileName();
int sbs_field = mainConfig.getSBSField();
TString partial_name_p = mainConfig.getPartialNameP();
TString partial_name_n = mainConfig.getPartialNameN();

//store all important kinematic info in local variables
kinematic_obj myKin(kinematic_file,kin);
double EBeam = myKin.getBeamEnergy();
double hcaldist = myKin.getHCalDist();
double sbsdist = myKin.getSBSDist();
double bbtheta = myKin.getBBAngle_Rad();
double hcaltheta = myKin.getHCalAngle_Rad();

//more setup stuff might need to go here

//setup output file
TString outfile = utility::makeOutputFileName_MCParse(exp,kin,sbs_field); 
TFile *fout = new TFile(outfile,"RECREATE");

//need to find the MC root and hist files 
//proton
vector<string> rootFileNames_p;
//find hist file first for proton
vector<string> histFileNames_p = utility::findHistFiles(replay_type,histfile_dir,partial_name_p);
//This function is a bit tricky as it will modify both the root file and hist file vectors. Not my favorite way to do that
//In the end we should have a set of matching hist files and root files
utility::matchMCFiles(replay_type,histFileNames_p,rootFileNames_p,rootfile_dir);

//neutron
vector<string> rootFileNames_n;
//find hist file first for neutron
vector<string> histFileNames_n = utility::findHistFiles(replay_type,histfile_dir,partial_name_n);
//Match hist and root files
utility::matchMCFiles(replay_type,histFileNames_n,rootFileNames_n,rootfile_dir);


// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
