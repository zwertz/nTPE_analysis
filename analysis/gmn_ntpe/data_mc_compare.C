
//Author Ezekiel Wertz
//04/09/2024
//Purpose of the script is to extract n/p yield ratio by comparing real with MC. This will be controlled by a particular kinematic and field setting. So it will not be implemented over multiple kinematics or field settings. The script will take in histograms and information from separate already parsed data and MC files. These files should already implement elastic cuts, an intime-clustering for HCal,and a fiducial cut. The goal is to then respectively fit the dx distributions to get a best fit. Then compose a sum of these fits allowing for scaling parameters which map to the total MC yields. Need to determine how best to handle the background, that is grounded in physics. Ultimately will apply the three floating parameter fit to the data and check residuals and chi-square. Should also figure out how to vary the fit function, ideally to yield better results. The ratio will then be extracted.

#include "TFile.h"
#include "../../src/utility.C"
#include "../../src/parse_config.C"

//Main
void data_mc_compare(const char *setup_file_name){

//Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//parse object to get in the information that The One Config file has and is manipulated
parse_config mainConfig(setup_file_name);
//mainConfig.printDataYields();

//store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
TString exp = mainConfig.getExp();
TString kin = mainConfig.getKin();
TString pass = mainConfig.getPass();
int sbs_field = mainConfig.getSBSField();
TString target = mainConfig.getTarg();

//setup input file info
TString Data_input_file_name = mainConfig.getDataFile();
TString MC_input_file_name = mainConfig.getMCFileName();

TFile *data_file = new TFile(Data_input_file_name.Data());
TFile *mc_file = new TFile(MC_input_file_name.Data());

//Setup output file
TString outfile = utility::makeOutputFileName_DataMCComp(exp,pass,kin,sbs_field,target);

//Get Histograms we will need from the respective data or MC files
TH1D *hdx_data = dynamic_cast<TH1D*>(data_file);










// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
