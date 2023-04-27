#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "analysis_utility_functions.h"
#include "BlindFactor.h"
#include <vector>


//Define max number of tracks per event
int MAXNTRACKS;

double tdiffmax, tdiff; // Max deviation from coin via tdctrig cut
TString Exp,kin,data_file_name,kinematic_file_name,targ;
int SBS_field,useAlshield;
double dx_low,dx_high,dy_low,dy_high,dxO_p,dyO_p,dxsig_p,dysig_p; //used for fitting proton information

//Use for Earm fit tuning
double W2_mean,W2_sigma; // min and max value in W2 plot range
double W2fitmax; //Max value in W2 plot range
double W2fitmaxwide; //Max value in wide W2 plot range
double binfac; //Number of bins/unit W2
double confac; //Number of sigma in delta x/y considered a "detection"
double W2confac; //Number of sigma in W2 considered elastic

//Use for Harm fit tuning
double hbinfac;//Number of bins per unit dx
double hcalemin; 
double dconfac; //Number of sigma outside of which inelastics only

//Define some fits to be used for analysis
//Background fit of a fourth-order polynomial
bool reject_bkgd;
double HCalBG_fit(double *x, double *param){

double e = param[0];
double d = param[1];
double c = param[2];
double b = param[3];
double a = param[4];

double func =a*(pow(x[0],4))+b*(pow(x[0],3))+ c*(pow(x[0],2))+d*(x[0])+e;
return func;

}

//Fit for protons using a Gaussian
double P_fit(double *x, double *param){
double amp = param[0];
double offset = param[1];
double sigma = param[2];

double func = amp*(exp(-0.5*pow((x[0]-offset)/sigma,2)));
return func;
}


//Total fit for HCal. Combine all fits
double HCalTot_fit(double *x, double *param){
double tot_func = HCalBG_fit(x,&param[0])+P_fit(x,&param[5]);
return tot_func;

}
//Background fit for EArm - use an expo 
double EBG_fit(double *x, double *param){
double amp = param[0];
double offset = param[1];
double str = param[2];

double func = amp*exp(offset+str*x[0]); 
return func;

}
//Fit for the E peak, assume a Gaussian
double E_fit(double *x, double *param){
double amp = param[0];
double offset = param[1];
double sigma = param[2];

double func = amp*(exp(-0.5*pow((x[0]-offset)/sigma,2)));
return func;
}

//Total fit for HCal. Combine all fits
double ETot_fit(double *x, double *param){
double tot_func = EBG_fit(x,&param[0])+E_fit(x,&param[3]);
return tot_func;

}


double poly4_fit(double *x, double *param){

double e = param[0];
double d = param[1];
double c = param[2];
double b = param[3];
double a = param[4];

double func =a*(pow(x[0],4))+b*(pow(x[0],3))+ c*(pow(x[0],2))+d*(x[0])+e;
return func;

}


TH1D *hW2elastic;
double W2total(double *x, double *par){ //get poly fit to bg with scaled fit to "pure elastics"
  double W2 = x[0];
  double sig_scale = par[0];
  double signal = sig_scale * hW2elastic->Interpolate(W2);
  return signal + poly4_fit(x,&par[1]);
}

TH1D *hW2elasticdxanticutresiduals;
Double_t W2totaldxanticutresiduals(double *x, double *par){ //get poly fit to bg with scaled fit to "pure elastics"
  double W2 = x[0];
  double sig_scale = par[0];
  double signal = sig_scale * hW2elasticdxanticutresiduals->Interpolate(W2);
  return signal + poly4_fit(x,&par[1]);
}


TCut globalcut = "";


//A vector to hold run numbers as TStrings eventually
vector<int> runnums;
//Make a function that parse the main config file
void parseMainConfig(const char *setup_file_name){
	ifstream setupfile(setup_file_name);
	//check if there is a problem opening the file
	 if(setupfile.fail()){
	cout << "There was a problem with the setup file " << setup_file_name << ". Figure it out nerd!" << endl;
	return;
	} 
	TString myline;
	//Look for the line where all the runs we care about are stored
	while(myline.ReadLine(setupfile) && !myline.BeginsWith("endrun")){
	 if(!myline.BeginsWith("#")){
	 TObjArray *demObjs = myline.Tokenize(" ");
		int numObjs = demObjs->GetEntries();
	 		for(int i=0;i < numObjs;i++){
	 		//store the run numbers in the vector as ints
	 		int runnum_temp = (((TObjString*) (*demObjs)[i])->GetString()).Atoi();
	 		runnums.push_back(runnum_temp);
			//cout << runnum_temp << " ";
			}
	
	 }
	}  
	while(myline.ReadLine(setupfile) && !myline.BeginsWith("endcut")){
        if(!myline.BeginsWith("#")){
         globalcut += myline;
	 cout << "Applying the following global cut to all data: " <<globalcut <<endl;
         }
        }
	while(myline.ReadLine(setupfile) && !myline.BeginsWith("#")){
        TObjArray *daObjs = myline.Tokenize(" ");
        int nObjs = daObjs->GetEntries();
		if(nObjs >1){
		TString key = ((TObjString*) (*daObjs)[0])->GetString();
		TString val = ((TObjString*) (*daObjs)[1])->GetString();
			if(key == "exp"){
			Exp = val;
			//cout << "Experiment " << Exp << endl;
			}
			else if(key == "kin"){
    			kin = val;            	
                        //cout << "Kinematic " << kin << endl;
			}
			else if(key == "data_map_name"){
                	data_file_name = val;
                        //cout << "Data File " << data_fiel_name << endl;
			}
			else if(key == "kinematic_name"){
                	kinematic_file_name = val;
                        //cout << "Kinematic File " << kinematic_file_name << endl;
			}
			else if(key == "SBS_field"){
                	SBS_field = val.Atoi();
                        //cout << "SBS Field " << SBS_field << endl;
			}
			else if(key == "W2_mean"){
                	W2_mean = val.Atof();
                        //cout << "W2 low " << W2_mean << endl;
			}
			else if(key == "W2_sigma"){
			W2_sigma = val.Atof();
                        //cout << "W2 high " << W2_sigma << endl;
			}
			else if(key == "targ"){
			targ = val;
			//cout << "Target " << targ << endl;
			}
			else if(key == "tdiff"){
			tdiff = val.Atof();
			//cout << "Time Diff " << tdiff << endl;
			}
			else if(key == "tdiffmax"){
			tdiffmax = val.Atof();
			//cout << "Time Cut " << tdiffmax << endl;
			}
			else if(key == "MAXNTRACKS"){
			MAXNTRACKS = val.Atoi();
			//cout << "Max Number of Tracks per event" << MAXNTRACKS << endl;
			}
			else if(key == "dxO_p"){
                        dxO_p = val.Atof();
                        //cout << "x-position of proton spot" << dxO_p << endl;
                        }
                        else if(key == "dyO_p"){
                        dyO_p = val.Atof();
                        //cout << "y-position of proton spot" << dyO_p << endl;
                        }
                        else if(key == "dxsig_p"){
                        dxsig_p = val.Atof();
                        //cout << "x sigma of proton spot" << dxsig_p << endl;
                        }                       
                        else if(key == "dysig_p"){
                        dysig_p = val.Atof();
                        //cout << "y sigma of proton spot" << dysig_p << endl;
                        }                        
			else if(key == "useAlshield"){
                        useAlshield = val.Atoi();
                        //cout << "Use Al shield" << useAlshield << endl;
                        }
			else if(key == "dx_low"){
                        dx_low = val.Atof();
                        //cout << "dx plot lower bound" << dx_low << endl;
                        }
			else if(key == "dx_high"){
                        dx_high = val.Atof();
                        //cout << "dx plot higher bound" << dx_high << endl;
                        }
			else if(key == "dy_low"){
                        dy_low = val.Atof();
                        //cout << "dy plot lower bound" << dy_low << endl;
                        }
                        else if(key == "dy_high"){
                        dy_high = val.Atof();
                        //cout << "dy plot higher bound" << dy_high << endl;
                        }
			else if(key == "binfac"){
			binfac = val.Atof();
			//cout << "binfac" << binfac << endl;
			}
			else if(key == "confac"){
			confac = val.Atof();
			//cout << "confac" << confac << endl;
			}
			else if(key == "W2confac"){
			W2confac = val.Atof();
			//cout << "W2confac" << W2confac << endl;
			}
			else if(key == "dconfac"){
			dconfac = val.Atof();
                        //cout << "dconfac" << dconfac << endl;
			}
			else if(key == "W2fitmax"){
                        W2fitmax = val.Atof();
			//cout << "W2fitmax" << W2fitmax << endl;
			}
			 else if(key == "W2fitmaxwide"){
                        W2fitmaxwide = val.Atof();
                        //cout << "W2fitmaxwide" << W2fitmaxwide << endl;
                        }
			else if(key == "hbinfac"){
                        hbinfac = val.Atof();
                        //cout << "hbinfac" << hbinfac << endl;
                        }
                        else if(key == "hcalemin"){
                        hcalemin = val.Atof();
                        //cout << "hcalemin" << hcalemin << endl;
                        }                        
			else{
			//We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
			cout << "Found a key: "<< key << " that this script can't handle. Fix that!" << endl;
			return;
			}
		//remove the objects to ensure a new comes through and no runaway
		delete daObjs; 
		}
		else{
		//We either got an empty line or 1 element.
		cout << "Line does not have the right number of elements. Look at the config file!" << endl;
		return;
		}
        }
	if(runnums.empty()){
	// if there are no run nums in the file, not good
	cout << "No run numbers in the config file, I can't do anything if I don't know where the data lives!" << endl;
	return;
	}
}

void HCalEfficiency( const char *setup_file_name){
 
 // Define a clock to check macro processing time
 TStopwatch *watch = new TStopwatch();
 watch->Start( kTRUE );     

 //Constructor for the data structure we will use to make plots
 TChain *C = new TChain("T");

 //parse function to get in the information that The One Config file has and is manipulated
 //this function will inialize the global parameters: globalcut,Exp,kin,data_file_name,kinematic_file_name,SBS_field,W2_low,W2_high,targ,runnums
 parseMainConfig(setup_file_name);

 //Load in information about the data files specified in The One config file. 
 //We will store the information in data_objects and put those in a vector
 vector<data_object> myData;
 for (int j=0; j < runnums.size();j++){
 data_object myObj(runnums[j],data_file_name,kinematic_file_name,kin,SBS_field,targ);
 myData.push_back(myObj);
//Use the class function getInputFile() to generate the file paths and add them to the TChain
 TString input_file_name(myObj.getInputFile());
 //cout << "File Location " << input_file_name << endl;
 C->Add(input_file_name);
 myObj.printRunInfo();
 }

  
  //variables we need  are BigBite track px,py,pz and sbs hcal x, y, e
  double Ebeam = myData[0].getBeamEnergy(); //Beam Energy, should be the same for all the data
  double hcaldist = myData[0].getHCalDist();// HCal Dist, should be the same for all the data
  double sbstheta = myData[0].getSBSAngle_Rad(); //SBS Angle, again should be the same for all data
  double hcaltheta = myData[0].getHCalAngle_Rad();//HCal Angle, same for all data being considered
  double bbtheta = myData[0].getBBAngle_Rad();//BB angle, same for all data
  int sbs_field = myData[0].getSBSField(); // SBS field
  double ntrack;
  
 //BBCal variables
  double BBtr_px[MAXNTRACKS], BBtr_py[MAXNTRACKS], BBtr_pz[MAXNTRACKS], BBtr_p[MAXNTRACKS];
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS],BBtr_chi2[MAXNTRACKS],BBtr_tgth[MAXNTRACKS],BBtr_tgph[MAXNTRACKS];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, GEMtr_hits, BBeoverp;  


  //HCal variables
  double atime[maxHCalChan],tdctime[maxHCalChan], row[maxHCalRows], col[maxHCalCols], cblkid[maxHCalChan], cblke [maxHCalChan];
  UInt_t TBits;
  double xhcal,yhcal,ehcal,nclus,nblk, hcal_atime,hcal_tdctime,bbcal_atime;
  double TDCT_id[maxTDCTrigChan], TDCT_tdc[maxTDCTrigChan], hodo_tmean[maxTDCTrigChan];
  int TDCTndata;
  double kineW2;
 

  //Cut on global parameters from setup config
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  //Setup leaves
  C->SetMakeClass(1); //Allows for viewing of general detector params from Event_Branch
  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.tg_th",1);
  C->SetBranchStatus("bb.tr.tg_ph",1);
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("e.kine.W2",1);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );
//need to change the bb time stuff
  C->SetBranchStatus( "sbs.hcal.atimeblk", 1 );
  C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 );

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.vz",vx);
  C->SetBranchAddress("bb.tr.vz",vy);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.tr.px",BBtr_px);
  C->SetBranchAddress("bb.tr.py",BBtr_py);
  C->SetBranchAddress("bb.tr.pz",BBtr_pz);
  C->SetBranchAddress("bb.tr.p",BBtr_p);
  C->SetBranchAddress("sbs.hcal.x",&xhcal);
  C->SetBranchAddress("sbs.hcal.y",&yhcal);
  C->SetBranchAddress("sbs.hcal.e",&ehcal);

  C->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", col );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );

  C->SetBranchAddress( "sbs.hcal.atimeblk", &hcal_atime );
  C->SetBranchAddress( "sbs.hcal.tdctimeblk", &hcal_tdctime );
  

  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "bb.tr.tg_th", BBtr_tgth );
  C->SetBranchAddress( "bb.tr.tg_ph", BBtr_tgph );
  C->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits ); 
  //need to change these bb times
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "bb.sh.atimeblk", &bbcal_atime );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  C->SetBranchAddress( "bb.gem.track.nhits", &GEMtr_hits );
  C->SetBranchAddress( "bb.etot_over_p", &BBeoverp );
  
  //Declare output file
  //Setup the name for the output file
  TString outfile = makeOutputFileName(Exp,kin,SBS_field,targ);
  TFile *fout = new TFile(outfile,"RECREATE");

 int nentries = C->GetEntries();
 //global variables for script
 
  double Emin = 0.0; // Minimum HCal energy to be considered as a hit, added for hard check on 1D distributions
  double efficiency_rel; // Relative detection efficiency of HCAL (elastic events detected by HCal) / (elastic events as defined by BB tracking)
  int hits_elBB = 0; // Count of total elastic hits detected in BB
  int hits_gHCAL = 0; // Count of all elastic events that HCal detected
  double pBeam = Ebeam/(1.0 +(Ebeam/M_p)*(1.0-cos(bbtheta)));
 //mean energy loss of the beam before scattering
 double Eloss_outgoing = (celldiameter/2)/sin(bbtheta)*rho_tgt*dEdx_tgt; //Approximately 1 MeV, one could correct more using the raster position
 if(useAlshield !=0){
 Eloss_outgoing += Alshieldthick * rho_Al*dEdx_Al;
 }
 //W2_low and W2_high. Should be equivalent to W2min_elastic and W2max_elastic. Also has a 2*W2confac range around some mean. This is tunable in the config file
 int hcalbins = 500;
 //Fit limits for HCal/W2 elastics
 //MC version
 double hcalfit_l_MC = posHCalXi_MC-2*HCalblk_l_v_MC;
 double hcalfit_h_MC = posHCalXf_MC+2*HCalblk_l_v_MC;
 double harmrange_MC = (hcalfit_h_MC) - (hcalfit_l_MC);
 //Detector version
 double hcalfit_l = posHCalXi-2*HCalblk_l;
 double hcalfit_h = posHCalXf+2*HCalblk_l;
 double harmrange = (hcalfit_h) - (hcalfit_l);
 double fit_l = 0.0;
 double fit_h = W2fitmax;

  //2sig 
  double dxmin = dxO_p - 2*dxsig_p;
  double dxmax = dxO_p + 2*dxsig_p;
  double dymin = dyO_p - 2*dysig_p;
  double dymax = dxO_p + 2*dxsig_p;

  //3sig dx
  double dxmin3 = dxO_p - 3*dxsig_p;
  double dxmax3 = dxO_p + 3*dxsig_p;

//Histograms
//HCal detections
  TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal X Expect vs Y Expect, no cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_hactivecut = new TH2D("hxy_hactivecut","HCal X Expect vs Y Expect, acceptance cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  
//No cuts
  TH1D *hW2_nocut = new TH1D( "hW2_nocut","W^{2}, no cut;W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_nocut_wide = new TH1D( "hW2_nocut_wide","W^{2}, no cut;W^{2} (GeV^{2})", binfac*W2fitmaxwide, 0.0, W2fitmaxwide );
  TH1D *hdx_nocut = new TH1D("hdx_nocut","dx, no cut;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_nocut_MC = new TH1D("hdx_nocut_MC","dx, no cut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hdy_nocut = new TH1D("hdy_nocut","dy, no cut;y_{HCAL}-y_{expect} (m)", hcalbins, posHCalYi-2*HCalblk_l, posHCalYf+2*HCalblk_l);
  TH1D *hdy_nocut_MC = new TH1D("hdy_nocut_MC","dy, no cut, MC param;y_{HCAL}-y_{expect} (m)", hcalbins, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC);
  TH2D *hdxdy_nocut = new TH2D("hdxdy_nocut","dx vs dy, no cut;dy_{HCAL} (m); dx_{HCAL} (m)", 250, posHCalYi-2*HCalblk_l, posHCalYf+2*HCalblk_l, 250, hcalfit_l, hcalfit_h );
  TH2D *hdxdy_nocut_MC = new TH2D("hdxdy_nocut_MC","dx vs dy, no cut, MC param;dy_{HCAL} (m); dx_{HCAL} (m)", 250, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC, 250, hcalfit_l_MC, hcalfit_h_MC );
  TH1D *hW2_allcut = new TH1D( "hW2_allcut","W^{2}, global/dy/dx/thetapq cut;W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );

//Cuts
  TH1D *hW2_BBcut = new TH1D( "hW2_BBcut","W^{2}, elastic cut;W^{2} (GeV^{2});", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_BBcut_norange = new TH1D( "hW2_BBcut_norange","W^{2}, elastic cut;W^{2} (GeV^{2});", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_BBcut_dycut = new TH1D( "hW2_BBcut_dycut","W^{2}, BBelas cut and dy cut;W^{2} (GeV)", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_BBcut_HCalcut = new TH1D( "hW2_BBcut_HCalcut","W^{2}, BBelas cut and HCal cut;W^{2} (GeV)", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hdx_BBcut = new TH1D("hdx_BBcut","dx, BBelas cut;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_BBcut_MC = new TH1D("hdx_BBcut_MC","dx, BBelas cut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hdy_BBcut = new TH1D("hdy_BBcut","dy, BBelas cut;y_{HCAL}-y_{expect} (m)", hcalbins, posHCalYi-2*HCalblk_l, posHCalYf+2*HCalblk_l);
  TH1D *hdy_BBcut_MC = new TH1D("hdy_BBcut_MC","dy, BBelas cut, MC param;y_{HCAL}-y_{expect} (m)", hcalbins, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC);
  TH2D *hdxdy_BBcut = new TH2D("hdxdy_BBcut","dxdy, BBelas cut;dy_{HCAL} (m); dx_{HCAL} (m)", 250, posHCalYi-2*HCalblk_l, posHCalYf+2*HCalblk_l, 250, hcalfit_l, hcalfit_h );
  TH2D *hdxdy_BBcut_MC = new TH2D("hdxdy_BBcut_MC","dxdy, BBelas cut, MC param;dy_{HCAL} (m); dx_{HCAL} (m)", 250, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC, 250, hcalfit_l_MC, hcalfit_h_MC );
  TH1D *hdx_BBcut_dycut = new TH1D("hdx_BBcut_dycut","dx, BBelas and dy cut;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_BBcut_dycut_MC = new TH1D("hdx_BBcut_dycut_MC","dx, BBelas and dy cut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hW2_dycut = new TH1D( "hW2_dycut","W^{2}, hcal dy cut;W^{2} (GeV^{2});", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hdx_allcut = new TH1D("hdx_allcut","dx, W^{2} and dy cut;x_{HCAL}-x_{expect} (m)",hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_allcut_MC = new TH1D("hdx_allcut_MC","dx, W^{2} and dy cut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  

  //Anticuts
  TH1D *hW2_badYcut = new TH1D( "hW2_badYcut","W^{2}, hcal dy anticut;W^{2} (GeV)", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hdx_badW2cut = new TH1D( "hdx_badW2cut","hcal dx, W^{2} anticut;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_badW2cut_MC = new TH1D( "hdx_badW2cut_MC","hcal dx, W^{2} anticut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hdx_badYcut = new TH1D( "hdx_badYcut","hcal dx, hcal dy anticut;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_badYcut_MC = new TH1D( "hdx_badYcut_MC","hcal dx, hcal dy anticut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hW2_fgcut = new TH1D( "hW2_fgcut","W^{2}, failed global anticut;W^{2} (GeV^{2});", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hdx_fgcut = new TH1D( "hdx_fgcut","hcal dx, failed global anticut;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_fgcut_MC = new TH1D( "hdx_fgcut_MC","hcal dx, failed global anticut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hW2_anticut = new TH1D( "hW2_anticut","W^{2} elastic anticut (global e-arm and hcal dy/dx);W^{2} (GeV^{2});", binfac*W2fitmax, 0.0, W2fitmax );

  //Diagnostic or useful quantities
  TH1D *hE_eloss = new TH1D("E_eloss","Scattered Electron Energy Loss in Target",500,0.0, Ebeam *0.001);
  TH1D *hE_ecorr = new TH1D("E_ecorr","Corrected Scattered Electron Energy",500,0.0,Ebeam*1.25);
  TH2D *hE_ecorr_vs_vertex = new TH2D("hE_ecorr_vs_vertex",";z_{vertex} (m); E_{corr} (GeV)", 250, -0.125, 0.125, 500, 0,0.001);
  TH1D *hE_ep = new TH1D("E_ep","Scattered Electron Energy",500,0.0,Ebeam*1.25);
  hE_ep->GetXaxis()->SetTitle("GeV");
  TH1D *hE_pp = new TH1D("E_pp","Scattered Proton Energy",500, 0.0, Ebeam*1.25);
  hE_pp->GetXaxis()->SetTitle("GeV");
  TH1D *hKE_p = new TH1D( "KE_p", "Scattered Proton Kinetic Energy", 500, 0.0, Ebeam*1.25 );
  hKE_p->GetXaxis()->SetTitle("GeV");

  TH1D *hQ2 = new TH1D("Q2","Q2",250,0.5,6.0);
  hQ2->GetXaxis()->SetTitle("GeV");
  TH1D *hvz = new TH1D( "hvz", "Vertex Position (m); m", 250, -0.2, 0.2 );
  TH1D *h_dpel_nocut = new TH1D("h_dpel_nocut","d_pel;p/p_{elastic}(#theta)-1;",250,-1.0,0.5);
  TH1D *h_dpel_BBcut = new TH1D("h_dpel_BBcut","d_pel;p/p_{elastic}(#theta)-1;",250,-1.0,0.5);
 
 //Accounting and diagnostic variables
 long nevent = 0;
 int treenum = 0, currenttreenum = 0,currentrun =0;
 double progress = 0.0;

 cout << "Processing events.." << endl;
 
 //create loop to display the progress bar
 while (progress < 1.0){
 //cout << progress<< " " << nevent<< " " <<nentries << endl;
 int barwidth = 70;
 int step =1;
 if(nevent >= nentries){
 //cout << "Hello" << endl;
 break;
 }
  while(C->GetEntry(nevent++)){
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
    //  cout << "My leaves" << endl;
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }
   bool failedglobal = GlobalCut->EvalInstance(0) == 0;

  
  //Correct the beam energy with energy loss in target using vertex position
  double Eloss = (vz[0]+(l_tgt/2))*rho_tgt*dEdx_tgt + uwallthick_LH2*rho_Al*dEdx_Al; //aproximately 3 MeV
  hE_eloss->Fill(Eloss);
  double Ecorr = Ebeam-Eloss;
  hE_ecorr->Fill(Ecorr);

  //Physics calculations
  double pcorr = BBtr_p[0]-Eloss_outgoing; //neglecting outgoing electron mass
  double p_ep = BBtr_p[0]; // Obtain the magnitude of scattered electron momentum
  double etheta = acos( BBtr_pz[0]/p_ep); //Use the uncorrected track momentum to reconstruct e' thetha
  double ephi = atan2( BBtr_py[0], BBtr_px[0] );
  TVector3 vertex( 0, 0, vz[0] ); //z location of vertex in hall coordinates
  TLorentzVector Pbeam(0,0,Ecorr,Ecorr);//Mass of e negligable
  TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
  TLorentzVector Ptarg( 0, 0, 0, M_p ); //Just use proton?
  TLorentzVector q = Pbeam - kprime; //Standard q-vector
  TLorentzVector PgammaN = Ptarg+q; // (-px, -py, Ebeam-pz,Mp+Ebeam-p)
  double E_ep = sqrt(pow(M_e,2)+pow(BBtr_p[0],2)); //Obtain the scattered electron energy
  hE_ep->Fill(E_ep);
  double pelastic = Ecorr /(1.0+(Ecorr/M_p)*(1.0-cos(etheta)));
  double nu = Ecorr-E_ep; // Obtain energy transfer
  double pp = sqrt(pow(nu,2)+ 2*M_p*nu);
  double phinucleon = ephi + TMath::Pi(); //assume coplanarity
  double thetanucleon = acos( (Ecorr - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics
  TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

   //Define HCal coordinate system
   TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
   TVector3 hcal_xaxis(0,-1,0);
   TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
   TVector3 hcal_origin = hcaldist *hcal_zaxis +HCALHeight*hcal_xaxis;
   TVector3 hcalpos = hcal_origin + xhcal * hcal_xaxis + yhcal * hcal_yaxis;
  //Define interesection points for hadron vector
  double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/pNhat.Dot( hcal_zaxis );
  TVector3 hcal_intersect = vertex + sintersect * pNhat;

  //Define the expected position of hadron on HCal from BB track
  double xhcal_expect = (hcal_intersect-hcal_origin).Dot( hcal_xaxis );
  double yhcal_expect = (hcal_intersect-hcal_origin).Dot( hcal_yaxis );

  //calculate quantities of interest
  double W2 = kineW2;//Get the invariant mass transfer W and four-momentum of the scattered nucleon
  //double W2c = PgammaN.M2();//M2 is magnitude squared of the 4-vector
  //double Wc = PgammaN.M();
  double Q2 = 2*Ecorr*E_ep*( 1-cos(etheta) );
  hQ2->Fill(Q2);
  //Use the electron kinematics to predict the proton momentum assuming elastic scattering on the free proton at rest (will need to correct for fermi motion);
  double E_pp = nu+M_p; //Get energy of proton
  double Enucleon = sqrt(pow(pp,2)+pow(M_p,2));
  hE_pp->Fill(E_pp);
  double KE_p = nu; //For elastics
  hKE_p->Fill(KE_p);
  double dpel = BBtr_p[0]/pelastic - 1.0;
  double precon = BBtr_p[0] + Eloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
  double nu_recon = Ecorr - precon;
  double Q2recon = 2.0*Ecorr*precon*(1.0-cos(etheta));
  double W2recon = pow(M_p,2) + 2.0*M_p*nu_recon - Q2recon;

    //Calculate q vector
    TLorentzVector Kprime_recon(precon*BBtr_px[0]/BBtr_p[0],precon*BBtr_py[0]/BBtr_p[0],precon*BBtr_pz[0]/BBtr_p[0],precon);
    TLorentzVector q_recon = Pbeam - Kprime_recon;
    TVector3 qvect_recon = q_recon.Vect();
    
    //Get theta pq for neutron and proton
    TVector3 NeutronDirection = (hcalpos - vertex).Unit();
    double BdL = SBSfield * maxSBSfield * Dgap * (sbs_field/100); //scale crudely by magnetic field
    double proton_deflection = tan( 0.3 * BdL / qvect_recon.Mag() ) * (hcaldist - (SBSdist + Dgap/2.0) );
    TVector3 ProtonDirection = (hcalpos + proton_deflection * hcal_xaxis - vertex).Unit();
    double thetapq = acos( ProtonDirection.Dot( qvect_recon.Unit() ) ); 
 
  //define dx,dy
  double dx = xhcal - xhcal_expect;
  double dy = yhcal - yhcal_expect;
 


  bool inboundsW2 = (W2>0.0) && (W2<W2fitmax);
  //Begin filling output histograms
  
  hxy_nocut->Fill(yhcal_expect,xhcal_expect);
  //No events which project hadrons off the face of HCal should be considered. Call this acceptance cut. This inherently relies on the MC DB values. Could also consider not that.
  bool offhcal =
      yhcal_expect > (posHCalYf_MC-HCalblk_l_h_MC) ||
      yhcal_expect < (posHCalYi_MC+HCalblk_l_h_MC) ||
      xhcal_expect > (posHCalXf_MC-HCalblk_l_v_MC) ||
      xhcal_expect < (posHCalXi_MC+HCalblk_l_v_MC);
 //Base level cut for all detection efficiency analysis
 if( offhcal ) continue;
 hxy_hactivecut->Fill(yhcal_expect,xhcal_expect);

 //all cut bools
 bool gW2 = abs(W2-W2_mean) < W2confac*W2_sigma; //W2 cut around elastic peak . One could implement this in W2_mean and W2_sigma form. But I am not sure this is better
 bool gdy = abs(dy-dyO_p)< confac*dysig_p; //dy cut around elastic peak

 
  //Fill histograms without most cuts, we always have an hcal-active-area/acceptance cut
  
  if(inboundsW2){
  hW2_nocut->Fill(W2); //Subtract dy anticut version from this
  }
  hW2_nocut_wide->Fill( W2 );
  hdxdy_nocut_MC->Fill( dy, dx );
  hdx_nocut_MC->Fill( dx );  //Subtract W2/dy anticut version from this
  hdy_nocut_MC->Fill( dy );
  h_dpel_nocut->Fill(dpel);
  hvz->Fill(vertex.Z());
 

 //anti-cuts
 if(!gW2){
  hdx_badW2cut_MC->Fill(dx); //Anticut dx histo
  }
  //having a continue causes a probelm somehow
  //anti-cuts
  if(!gdy ){
  hW2_badYcut->Fill(W2); //Anticut W2 histo
  hdx_badYcut_MC->Fill(dx); //Alt anticut dx histo
  }


  //e-arm
  if(gW2 && !failedglobal){
  hW2_BBcut->Fill(W2);
  hW2_BBcut_norange->Fill( W2 );
  hdx_BBcut_MC->Fill( dx );
  hdy_BBcut_MC->Fill( dy );
  hdxdy_BBcut_MC->Fill( dy, dx );
  h_dpel_BBcut->Fill(dpel);
  }
 
  //h-arm
  if(gdy && !failedglobal && (ehcal>hcalemin) ){
  hW2_dycut->Fill( W2 );

  }
  
  //both arms
  if(gW2 && gdy && !failedglobal  && (ehcal>hcalemin)){
  hdx_BBcut_dycut_MC->Fill( dx ); //H-arm elastics (efficiency numerator)
  hW2_BBcut_dycut->Fill( W2 );

  }
  
  
  //Make strongest possible hcal-only anticut to obtain background
  if( thetapq>0.04&&(dx<dxmin3||dx>dxmax3)&&(dy<dymin||dy>dymax)&&failedglobal==1 ){
  hW2_anticut->Fill( W2 );
  }

  //And tightest cut possible
  if( thetapq<0.04&& (dx>dxmin3&&dx<dxmax3)&&(dy>dymin&&dy<dymax)&&failedglobal==0 ){
  hW2_allcut->Fill( W2 );
  }                                


  cout << "[";
  int pos = barwidth*progress;
  
  for(int i=0; i<barwidth; ++i){
  //cout << step << endl;
  if(i<pos){
  cout << "_";
  }else if(i==pos){
   if(step%4==0){
   cout << "(>0_0)>";
   }
   if(step%4==1){
   cout << "<(0_0)>";
   }
   if(step%4==2){
   cout << "<(0_0<)";
   }
   if(step%4==3){
   cout << "<( ; )>";
   }
  }
  else{
  cout << " ";
  }
  }
 progress = (double) ((nevent+1.0)/nentries);
 cout << "]" << (int) (progress*100) << "%\r";
 cout.flush();
 if(nevent%10000==0){
 step++;
 }
 //cout << progress<< " " << nevent<< " " <<nentries << endl;
 
 }
}// This ends the double loop over the events

 //dy anti-cut BG subtraction methods
 TH1D *hdx_nocut_MC_clone = (TH1D*)(hdx_nocut_MC->Clone("hdx_nocut_MC_clone")); //for later fit option

 //Make canvas for W2 dy anticut analysis
 TCanvas *c1 = new TCanvas("c1",Form("hcal_sbs%s%smag%i_W2_dyanticut",kin.Data(), targ.Data(), sbs_field),1200,500);
 gStyle->SetPalette(53);
 c1->Divide(2,1);
 c1->cd(1);
        
 //Get difference between gY and !gY W2 after scaling to normalize
 hW2_nocut->SetLineColor(kBlack);
 hW2_nocut->SetTitle("W2");
 hW2_nocut->Draw();
 double gY_lastbin = hW2_nocut->GetBinContent(binfac*W2fitmax);
 double bY_lastbin = hW2_badYcut->GetBinContent(binfac*W2fitmax);

 double BGfac = gY_lastbin/bY_lastbin;
 //cout<< endl <<gY_lastbin << " " << bY_lastbin << " " << "e-arm dy anticut background factor: " << BGfac << endl;

 //Scale W2 "background" dy anticut histo up/down to match last bin of W2 nocut histo
 TH1D *hW2_badYcut_scaled = (TH1D*)(hW2_badYcut->Clone("hW2_badYcut_scaled"));
 hW2_badYcut_scaled->Scale(BGfac);
 hW2_badYcut_scaled->SetLineColor(kRed);
 hW2_badYcut_scaled->Draw("same hist");
 TH1D *hW2_BGsub = (TH1D*)(hW2_nocut->Clone("hW2_BGsub"));
 hW2_BGsub->Add(hW2_badYcut_scaled,-1);
 for( int i=1; i<=binfac*W2fitmax; i++ ){
 double binc = hW2_BGsub->GetBinContent(i);
  if( binc<0.0 ){ 
  hW2_BGsub->SetBinContent(i,0.0);
  }
 }
 //Add a legend to the canvas
 auto W2bgsublegend = new TLegend(0.1,0.8,0.5,0.9);
 W2bgsublegend->SetTextSize( 0.03 );
 W2bgsublegend->AddEntry( hW2_nocut, "No Cut", "l");
 W2bgsublegend->AddEntry( hW2_badYcut_scaled, Form("dy<%i*sigma anticut (scaled)",(int)confac), "l");
 W2bgsublegend->Draw();
           
 c1->cd(2);
 hW2_BGsub->SetLineColor(kBlue);
 hW2_BGsub->SetTitle("W2, BG Subtracted");
 hW2_BGsub->Draw("hist");

 double elastics_alt = hW2_BGsub->Integral(); //integrate all events in histogram after subtraction
 double eel_bincontents = 0.;

 for( int i=1; i<=binfac*W2fitmax; i++ ){ //bins start at 1
  eel_bincontents += hW2_BGsub->GetBinContent(i);
 }
 c1->Write();
 

 //Make canvas for hcal W2 anticut analysis
 TCanvas *c2 = new TCanvas("c2",Form("hcal_sbs%s%smag%i_dx_W2anticut",kin.Data(), targ.Data(), sbs_field),1200,500);
 gStyle->SetPalette(53);
 c2->Divide(2,1);
 c2->cd(1);

 //Get diff gW2 and !gW2 dx after scaling to normalize
 hdx_nocut_MC->SetLineColor(kBlack);
 hdx_nocut_MC->SetTitle("HCal dx");
 hdx_nocut_MC->Draw();

 double gW2_hlastbin = hdx_nocut_MC->GetBinContent(hcalbins);
 //Double_t gW2_hlastbin = hdx_HCAL_BBcut_dycut->GetBinContent(hcalbins);
 double bW2_hlastbin = hdx_badW2cut_MC->GetBinContent(hcalbins);
  
 double hcalBGfac = gW2_hlastbin/bW2_hlastbin;
 cout << endl << "h-arm W2 anticut background factor: " << hcalBGfac << endl;

 //Scale dx "background" W2 anticut histo up/down to match last bin of dx nocut histo
 TH1D *hdx_badW2cut_scaled = (TH1D*)(hdx_badW2cut_MC->Clone("hdx_badW2cut_MC_scaled"));
 hdx_badW2cut_scaled->Scale(hcalBGfac);
 hdx_badW2cut_scaled->SetLineColor(kRed);
 hdx_badW2cut_scaled->Draw("same hist");
 TH1D *hdx_BGsub = (TH1D*)(hdx_nocut_MC->Clone("hdx_BGsub"));
 //TH1D *hdx_BGsub = (TH1D*)(hdx_HCAL_BBcut_dycut->Clone("hdx_BGsub"));
 hdx_BGsub->Add(hdx_badW2cut_scaled,-1);
 for( int j=1; j<=hcalbins; j++ ){
 double binc = hdx_BGsub->GetBinContent(j);
  if( binc<0.0 ){
  hdx_BGsub->SetBinContent(j,0.0);
  }
 }

 //Add a legend to the canvas
 auto dxbgsublegend = new TLegend(0.1,0.8,0.5,0.9);
 dxbgsublegend->SetTextSize( 0.03 );
 dxbgsublegend->AddEntry( hdx_nocut_MC, "No Cut", "l");
 //dxbgsublegend->AddEntry( hdx_HCAL_BBcut_dycut, "BB elastic cut, dy cut", "l");
 dxbgsublegend->AddEntry( hdx_badW2cut_scaled, Form("W2<%i*sig anticut (scaled)",(int)W2confac), "l");
 dxbgsublegend->Draw();

  c2->cd(2);
  hdx_BGsub->SetLineColor(kBlue);
  hdx_BGsub->SetTitle("HCal dx, BG subtracted (W2 anticut)");
  hdx_BGsub->Draw("hist");

  double hcalelastics_antiW2 = hdx_BGsub->Integral();

  double hel_bincontents = 0.0;
  for( int i=1; i<=hcalbins; i++ ){ //bins start at 1
    hel_bincontents += hdx_BGsub->GetBinContent(i);
  }

  c2->Write();
  
  //Make canvas for hcal dy anticut analysis
  TCanvas *c0 = new TCanvas("c0",Form("hcal_sbs%s%smag%i_dx_dyanticut",kin.Data(), targ.Data(),sbs_field),1200,500);
  gStyle->SetPalette(53);
  c0->Divide(2,1);
  c0->cd(1);
  
  //Get diff gY and !gY dx after scaling to normalize
  hdx_nocut_MC->SetLineColor(kBlack);
  hdx_nocut_MC->SetTitle("HCal dx");
  hdx_nocut_MC->Draw();
         
  double gY_hlastbin = hdx_nocut_MC->GetBinContent(hcalbins);
  double bY_hlastbin = hdx_badYcut_MC->GetBinContent(hcalbins);
  
  double hcalYBGfac = gY_hlastbin/bY_hlastbin;
  cout << endl << "h-arm dy anticut background factor: " << hcalYBGfac << endl;      

  //Scale dx "background" W2 anticut histo up/down to match last bin of dx nocut histo
  TH1D *hdx_badYcut_scaled = (TH1D*)(hdx_badYcut_MC->Clone("hdx_badYcut_MC_scaled"));
  hdx_badYcut_scaled->Scale(hcalYBGfac);
  hdx_badYcut_scaled->SetLineColor(kRed);
  hdx_badYcut_scaled->Draw("same hist");
  TH1D *hdx_YBGsub = (TH1D*)(hdx_nocut_MC->Clone("hdx_YBGsub"));
  //TH1D *hdx_YBGsub = (TH1D*)(hdx_HCAL_BBcut_dycut->Clone("hdx_YBGsub"));
  hdx_YBGsub->Add(hdx_badYcut_scaled,-1);
  for( int i=1; i<=hcalbins; i++ ){
  double binc = hdx_YBGsub->GetBinContent(i);
   if( binc<0.0 ){
   hdx_YBGsub->SetBinContent(i,0.0);
   }
  }
 
  //Add a legend to the canvas
  auto dxabgsublegend = new TLegend(0.1,0.8,0.5,0.9);
  dxabgsublegend->SetTextSize( 0.03 );
  dxabgsublegend->AddEntry( hdx_nocut_MC, "No Cut", "l");
  //dxabgsublegend->AddEntry( hdx_HCAL_BBcut_dycut, "BB elastic cut, dy cut", "l");
  dxabgsublegend->AddEntry( hdx_badYcut_scaled, Form("dy<%i*sig anticut (scaled)",(int)confac), "l");
  dxabgsublegend->Draw();
  
  c0->cd(2);
  hdx_YBGsub->SetLineColor(kBlue);
  hdx_YBGsub->SetTitle("Hcal dx, BG subtracted (dy anticut)");
  hdx_YBGsub->Draw("hist");
 
  double hcalelastics_antidy = hdx_YBGsub->Integral();

  double hely_bincontents = 0.0;
  for( int i=1; i<=hcalbins; i++ ){ //bins start at 1
  hely_bincontents += hdx_YBGsub->GetBinContent(i);
  }

  c0->Write();

  //Make canvas for dxdy to show dy cut
  TCanvas *c3 = new TCanvas("c3",Form("hcal_sbs%s%smag%i",kin.Data(), targ.Data(), sbs_field),650,500);
  gStyle->SetPalette(53);
  //c3->cd();
  //hdxdy_nocut_MC->GetYaxis()->SetTitleOffset(1.3); //can use BBcut
  hdxdy_nocut_MC->Draw("colz");
  TLine l_pl; //left cut
  l_pl.SetLineColor(kBlue);
  l_pl.SetLineWidth(3);
  l_pl.DrawLine(dyO_p-confac*dysig_p,hcalfit_l,dyO_p-confac*dysig_p,hcalfit_h);
  TLine l_pr; //right cut
  l_pr.SetLineColor(kBlue);
  l_pr.SetLineWidth(3);
  l_pr.DrawLine(dyO_p+confac*dysig_p,hcalfit_l,dyO_p+confac*dysig_p,hcalfit_h);
  c3->Write();
  
  //Make canvas to hold fitted W2 and HCal dx histograms
  TCanvas *c4 = new TCanvas("c4",Form("W2_sbs%s%smag%i",kin.Data(), targ.Data(), sbs_field),1200,500);
  c4->Divide(2,1);
  gStyle->SetPalette(53);

  //HCal dx fit
  c4->cd(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
 
 
  vector<double> myParam (13),myFParam (13);
  myParam= fit_Params(kin,sbs_field,targ);
  
   //Declare total fit function
   TF1 *hcaltotalfit = new TF1("hcaltotalfit",HCalTot_fit,hcalfit_l_MC,hcalfit_h_MC,8); //total fit
   
   //Set the initial fit params
   hcaltotalfit->SetParameters(&myParam[0]);
   hcaltotalfit->SetParLimits(5,10000,100000);
   hcaltotalfit->SetParLimits(6,-1.,0.0);
   hcaltotalfit->SetParLimits(7,0.05,0.5);

   //Fit the dy cut HCal dx histogram with total fit
   hcaltotalfit->SetLineColor(kGreen);
   hcaltotalfit->SetLineWidth(4);
   hdx_BBcut_dycut_MC->Fit("hcaltotalfit","RBM"); //Consider dx distribution in HCal with full elastic selection including dy

   //Get parameters from total fit and generate identical bg and elastic functions to write on canvas
   double *hcaltotpar = hcaltotalfit->GetParameters();
   
   TF1 *hcalbg = new TF1("hcalbg",HCalBG_fit,hcalfit_l_MC,hcalfit_h_MC,5); //inelastic pol4 background
   hcalbg->SetParameters(&hcaltotpar[0]);
   hcalbg->SetLineColor(kRed);
   hcalbg->Draw("same");

   TF1 *p = new TF1("p",P_fit,hcalfit_l_MC,hcalfit_h_MC,3); //proton elastic gaussian
   p->SetParameters(&hcaltotpar[5]);
   p->SetLineColor(kBlack);
   p->Draw("same");
   
   //Get elastics/inelastics from normalized integrals of fit functions
   double hcalinelastics = (hcalbins/harmrange)*hcalbg->Integral(hcalfit_l,hcalfit_h);
   double hcalelastics = (hcalbins/harmrange)*p->Integral(hcalfit_l,hcalfit_h);
   
   //Bigbite W2 fit
   c4->cd(2);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
  
   //Declare initial fit params for W2
   vector<double> W2params;
   double glimh[3];
   double gliml[3];
   double earmEvents = hW2_BBcut->GetEntries();
  
   //Set some reasonable limits for gaussian fit to elastic peak
   gliml[0] = earmEvents/binfac/2; //Amplitude constrained by number of events and binfac
   glimh[0] = earmEvents;
   gliml[1] = 0.75; //Mean constrained by W2 peak for elastics
   glimh[1] = 1.25;
   gliml[2] = 0.05; //Empirical limits
   glimh[2] = 0.50;
  
   //W2 overall fits
   TF1 *W2_totalfit;
   TF1 *W2_bg;
   TF1 *elastic = new TF1("elastic",E_fit,fit_l,fit_h,3); //The elastic peak always gaussian
   
   //expo fit
   W2params.push_back(465.0);
   W2params.push_back(-3.5);
   W2params.push_back(3.5);
   W2params.push_back(3057.);
   W2params.push_back(0.8077);
   W2params.push_back(0.0961);
   
   W2_totalfit = new TF1("W2_totalfit",ETot_fit,fit_l,fit_h,6);
   W2_bg = new TF1("W2_bg",EBG_fit,fit_l,fit_h,3);
   W2_totalfit->SetParameter(0,W2params[0]);
   W2_totalfit->SetParLimits(0,0.,10000.);
   W2_totalfit->SetParameter(1,W2params[1]);
   W2_totalfit->SetParLimits(1,-10000.,0.);
   W2_totalfit->SetParameter(2,W2params[2]);
   W2_totalfit->SetParLimits(2,0.,10000.);
   W2_totalfit->SetParameter(3,W2params[3]);
   W2_totalfit->SetParameter(4,W2params[4]);
   W2_totalfit->SetParameter(5,W2params[5]);
   W2_totalfit->SetParLimits(3,gliml[0],glimh[0]);
   W2_totalfit->SetParLimits(4,gliml[1],glimh[1]);
   W2_totalfit->SetParLimits(5,gliml[2],glimh[2]);

   //Fit the W2 histogram with total fit
   W2_totalfit->SetLineColor(kGreen);
   W2_totalfit->SetLineWidth(4);
   hW2_nocut->Fit("W2_totalfit","RBM"); //Fit the W2 distribution with total fit
   
   //Get all parameters
   double *totpar = W2_totalfit->GetParameters();
   
   //Get fit parameters for bg function and draw identical function on canvas
   W2_bg->SetParameters(&totpar[0]);
   W2_bg->SetLineColor(kRed);
   W2_bg->Draw("same");
   
   //Get fit parameters for e gaussian and draw identical gaussian on canvas
   elastic->SetParameters(&totpar[3]);
   elastic->SetLineColor(kBlack);
   elastic->Draw("same");
   
   //Integrate the elastic peak to get total elastics detected in BB
   double inelastics = binfac*W2_bg->Integral(0.0,W2fitmax); //scale with binfac
   double elastics = binfac*elastic->Integral(0.0,W2fitmax);
   
   //Add a legend
   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->SetTextSize( 0.03 );
   legend->SetHeader( "W2, Electron Arm Cuts", "C" );
   legend->AddEntry( W2_bg, Form("Inelastic Fit: %i ev", (int)inelastics), "l");
   legend->AddEntry( elastic, Form("Elastic Fit: %i ev", (int)elastics), "l");
   legend->AddEntry( W2_totalfit, "Total Fit", "l");
   legend->Draw();
 
   c4->Write();

   TH1D *hW2_totfit = (TH1D*)(hW2_nocut->Clone("hW2_totfit"));

   TH1D *hdx_totfit;
   TH1D *hW2_totfitdxanticutresiduals = (TH1D*)(hW2_anticut->Clone("hW2_totfitdxanticutresiduals"));
    //Make canvas for hcal W2 anticut analysis
    TCanvas *c5 = new TCanvas("c5",Form("W2 dx anticut residuals with extracted BG sbs%s mag%i",kin.Data(), sbs_field),1200,500);
    gStyle->SetPalette(53);
    c5->cd();
    
    //Add background extraction
    hW2elasticdxanticutresiduals = (TH1D*)(hW2_allcut->Clone("hW2elasticdxanticutresiduals"));
    TF1 *totfitdxanticutresiduals = new TF1("totfitdxanticutresiduals",W2totaldxanticutresiduals,0,W2fitmax,6); //tfit npar = 1+pNfit_npar+1
    TF1 *bgdxanticutresiduals = new TF1("bgdxanticutresiduals",poly4_fit,0.,W2fitmax,5);
    totfitdxanticutresiduals->SetLineColor(kGreen);
    hW2_totfitdxanticutresiduals->Fit("totfitdxanticutresiduals","RBM");
    
    double *totpardxanticutresiduals = totfitdxanticutresiduals->GetParameters();
    
    //Get fit parameters for bg function and draw identical function on canvas
    bgdxanticutresiduals->SetParameters(&totpardxanticutresiduals[1]);
    bgdxanticutresiduals->SetLineColor(kRed);
    bgdxanticutresiduals->SetFillColor(kRed);
    bgdxanticutresiduals->SetFillStyle(3005);
    bgdxanticutresiduals->Draw("same");
    hW2elasticdxanticutresiduals->SetLineColor(kBlue);
    hW2elasticdxanticutresiduals->SetLineWidth(2);
    hW2elasticdxanticutresiduals->SetFillColor(kBlue);
    hW2elasticdxanticutresiduals->SetFillStyle(3003);
    hW2elasticdxanticutresiduals->Draw("same");
    

     //calculate some qtys (elastic divergence from bg begin: 0.5, end: 1.3)
     //Double_t bgint = bg->Integral(0.,W2fitmax)*binfac; //175715
     double bgint_dxanticutresiduals = bgdxanticutresiduals->Integral(0.4,1.3)*binfac; //175715
     //int W2elasb = 0.0*binfac;
     int W2elasb_dxanticutresiduals = 0.4*binfac;
     //Int_t W2elase = W2fitmax*binfac;
     int W2elase_dxanticutresiduals = 1.3*binfac;
     double W2nocutint_dxanticutresiduals = hW2_totfitdxanticutresiduals->Integral(W2elasb_dxanticutresiduals,W2elase_dxanticutresiduals);
     //Double_t tfitint = tfit->Integral(0.0,W2fitmax)*binfac;
     double totfitint_dxanticutresiduals = totfitdxanticutresiduals->Integral(0.4,1.3)*binfac;
     double W2elas_dxanticutresiduals = W2nocutint_dxanticutresiduals - bgint_dxanticutresiduals;
     double W2elasfits_dxanticutresiduals = totfitint_dxanticutresiduals - bgint_dxanticutresiduals;
    
    
    //Add a legend to the canvas
    auto interplegenddxanticutresiduals = new TLegend(0.1,0.6,0.5,0.9);
    interplegenddxanticutresiduals->SetTextSize( 0.03 );
    interplegenddxanticutresiduals->AddEntry( hW2elasticdxanticutresiduals, "Tight Elastic Cut", "l");
    interplegenddxanticutresiduals->AddEntry( bgdxanticutresiduals, "Interpolated Background (scaled)", "l");
    interplegenddxanticutresiduals->AddEntry( totfitdxanticutresiduals, "Total fit (background + signal)", "l");
    interplegenddxanticutresiduals->Draw();
      
    c5->Write();  



    //Make canvas for hcal W2 anticut analysis
    TCanvas *c6 = new TCanvas("c6",Form("W2 with extracted BG sbs%s mag%i",kin.Data(), sbs_field),1200,500);
    gStyle->SetPalette(53);
    c6->cd();
    
    //Add background extraction
    hW2elastic = (TH1D*)(hW2_allcut->Clone("hW2elastic"));
    TF1 *W2_totfit = new TF1("W2_totfit",W2total,0,W2fitmax,6); //tfit npar = 1+pNfit_npar+1
    TF1 *W2_bkgd = new TF1("W2_bkgd",poly4_fit,0.,W2fitmax,5);
    W2_totfit->SetLineColor(kGreen);
    hW2_totfit->Fit("W2_totfit","RBM");
    
    double *tpar = W2_totfit->GetParameters();
    
    //Get fit parameters for bg function and draw identical function on canvas
    W2_bkgd->SetParameters(&tpar[1]);
    W2_bkgd->SetLineColor(kRed);
    W2_bkgd->SetFillColor(kRed);
    W2_bkgd->SetFillStyle(3005);
    W2_bkgd->Draw("same");
    hW2elastic->SetLineColor(kBlue);
    hW2elastic->SetLineWidth(2);
    hW2elastic->SetFillColor(kBlue);
    hW2elastic->SetFillStyle(3003);
    hW2elastic->Draw("same");
   
 
    //calculate some qtys (elastic divergence from bg begin: 0.4, end: 1.3)
    //Double_t bgint = bg->Integral(0.,W2fitmax)*binfac; //175715
    double W2bgint = W2_bkgd->Integral(0.4,1.3)*binfac; //175715
    //Int_t W2elasb = 0.*binfac;
    int W2elasb = 0.4*binfac;
    //Int_t W2elase = W2fitmax*binfac;
    int W2elase = 1.3*binfac;
    double W2nocutint = hW2_totfit->Integral(W2elasb,W2elase);
    //Double_t tfitint = tfit->Integral(0.0,W2fitmax)*binfac;
    double totfitint = W2_totfit->Integral(0.4,1.3)*binfac;
    double  W2elas = W2nocutint - W2bgint;
    double W2elas_fits = totfitint - W2bgint;
    
    double effres = ((W2elas-W2elas_dxanticutresiduals) / W2elas)*100.;


  //Add a legend to the canvas
  auto interplegend = new TLegend(0.1,0.6,0.5,0.9);
  interplegend->SetTextSize( 0.03 );
  interplegend->AddEntry( hW2elastic, "Tight Elastic Cut", "l");
  interplegend->AddEntry( W2_bkgd, "Interpolated Background (scaled)", "l");
  interplegend->AddEntry( W2_totfit, "Total fit (background + signal)", "l");
  interplegend->AddEntry( (TObject*)0, "", "");
  interplegend->AddEntry( (TObject*)0, Form("HCal Detection Efficiency: %0.3f%%",effres), "");
  interplegend->Draw();
  
  c6->Write();  

     //Calculate error assuming fit error is negligable, might need to reevaluate that later
     //Binomial
     double eff = hcalelastics/elastics;
     double eff_alt = hcalelastics_antiW2/elastics_alt;
     double eff_alt2 = hcalelastics_antidy/elastics_alt;
     //Double_t eff_althybrid = hcalelastics/elastics_alt;
     double efferr = sqrt(eff*(1-eff)/elastics);
     double efferr_alt = sqrt(eff_alt*(1-eff_alt)/elastics_alt);
     double efferr_alt2 = sqrt(eff_alt2*(1-eff_alt2)/elastics_alt);
     //Double_t efferr_althybrid = sqrt(eff_althybrid*(1-eff_althybrid)/elastics_alt);
     
   TString reportfile = makeReportFileName(Exp,kin,SBS_field,targ);   


    
   //Declare outfile
   ofstream report;
   
   cout << "Chi^2 for total fit: " << W2_totalfit->GetChisquare() << endl;
   report << "Chi^2 for total fit: " << W2_totalfit->GetChisquare() << endl;
   
   report.open( reportfile );
   report << "HCal detection efficiency report for SBS" << kin << " LH2 data at " << sbs_field << " percent field" << endl << endl;
   
   cout << "Total elastics detected in HCal: " << hcalelastics << endl;
   report << "Total elastics detected in HCal: " << hcalelastics << endl << endl;
   
   cout << "Total elastics detected in Bigbite with e-arm cuts only: " << elastics << endl;
   report << "Total elastics detected in Bigbite with e-arm cuts only: " << elastics << endl << endl;
   
   cout << "Total elastics detected in HCal (antidy): " << hcalelastics_antidy << endl;
   report << "Total elastics detected in HCal (antidy): " << hcalelastics_antidy << endl << endl;
  
   cout << "Total elastics detected in HCal (alt): " << hcalelastics_antiW2 << endl;
   report << "Total elastics detected in HCal (alt): " << hcalelastics_antiW2 << endl << endl;
   
   cout << "Total elastics detected in Bigbite with e-arm cuts only (alt): " << elastics_alt << endl;
   report << "Total elastics detected in Bigbite with e-arm cuts only (alt): " << elastics_alt << endl << endl;
  
   cout << "Total entries in W2 plot: " << earmEvents << endl << endl;
   report << "Total entries in W2 plot: " << earmEvents << endl << endl;
   
   cout << "Total elastics from W2 Tight elastic cut (global, dx, dy, thetapq) " << W2elas << endl;
   report << "Total elastics from W2 Tight elastic cut (global, dx, dy, thetapq) " << W2elas << endl;

    cout << "Total elastics from W2 residual Tight elastic cut (global, dx, dy, thetapq) " << W2elas_dxanticutresiduals << endl;
   report << "Total elastics from W2 residual Tight elastic cut (global, dx, dy, thetapq) " << W2elas_dxanticutresiduals << endl;
 
  cout << "==============================================================" << endl;
  cout << "Total HCal elastic detection efficiency (fits): " << eff << " +/- " << efferr << endl;
  cout << "==============================================================" << endl << endl;

  cout << "===================================================================" << endl;
  cout << "Total HCal elastic detection efficiency (anticut W2): " << eff_alt << " +/- " << efferr_alt << endl;
  cout << "===================================================================" << endl << endl;

  cout << "===================================================================" << endl;
  cout << "Total HCal elastic detection efficiency (anticut dy): " << eff_alt2 << " +/- " << efferr_alt2 << endl;
  cout << "===================================================================" << endl << endl;

  // cout << "===================================================================" << endl;
  // cout << "Total HCal elastic detection efficiency (althybrid): " << eff_althybrid << " +/- " << efferr_althybrid << endl;
  // cout << "===================================================================" << endl << endl;
  report << "==============================================================" << endl;
  report << "Total HCal elastic detection efficiency (fits): " << eff << " +/- " << efferr <<  endl;
  report << "==============================================================" << endl;
  
  report << "===================================================================" << endl;
  report << "Total HCal elastic detection efficiency (anticut W2): " << eff_alt << " +/- " << efferr_alt << endl;
  report << "===================================================================" << endl << endl;
  
  report << "===================================================================" << endl;
  report << "Total HCal elastic detection efficiency (anticut dy): " << eff_alt2 << " +/- " << efferr_alt2 << endl;
  report << "===================================================================" << endl << endl;

  cout << "===================================================================" << endl;
  cout << "Total HCal elastic detection efficiency  Tight elastic cut (global, dx, dy, thetapq): " << effres << " +/- " << "" << endl;
  cout << "===================================================================" << endl << endl;

  report << "===================================================================" << endl;
  report << "Total HCal elastic detection efficiency  Tight elastic cut (global, dx, dy, thetapq): " << effres << " +/- " << "" << endl;
  report << "===================================================================" << endl << endl;



  cout << endl << endl << "raw bincontent efficiency anticut dy: " << hely_bincontents / eel_bincontents << endl;
  cout << endl << endl << "raw bincontent efficiency anticut W2: " << hel_bincontents / eel_bincontents << endl;

  report.close();


  //Print the plots to a pdf file for safe keeping
  TString plotname = outfile;
  plotname.ReplaceAll(".root",".pdf");
  TString start = Form("%s%s",plotname.Data(),"(");
  //middle is the same as the name
  TString end = Form("%s%s",plotname.Data(),")");
  c1->Print(start.Data(),"pdf");
  c2->Print(plotname.Data(),"pdf");
  c0->Print(plotname.Data(),"pdf");
  c3->Print(plotname.Data(),"pdf");  
  c4->Print(plotname.Data(),"pdf");
  c5->Print(plotname.Data(),"pdf");
  c6->Print(end.Data(),"pdf");  

  fout->Write();
//TH1D *hdx_residual = (TH1D*)hdx_cut->Clone("hdx_residual");
//hdx_residual->SetTitle("HCal dx residual (data - fit),All cuts");
//hdx_residual->GetXaxis()->SetTitle("m");
  





  watch->Stop();

  // Send time efficiency report to console
   cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;


} 
