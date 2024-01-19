//Author: E. Wertz
//Script using W2 anticut method to determine HCal proton efficiency
//Includes acceptance matching, tight elastic cut based on thetapq/dx/dy, in-time clustering algorithm, quadrature sum error for statistics

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
#include "TMatrixD.h"

//Define max number of tracks per event
int MAXNTRACKS;
//for in-time
double coin_mean,coin_sigma,coin_sig_fac;
TString Exp,kin,data_file_name,kinematic_file_name,targ;
int SBS_field,useAlshield,first_event,last_event;
double dx_low,dx_high,dy_low,dy_high,dxO_p,dyO_p,dxsig_p,dysig_p; //used for fitting proton information
double thetapq_low,thetapq_high;
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

TH1D *hW2elastic;
double W2total(double *x, double *par){ //get poly fit to bg with scaled fit to "pure elastics"
  double W2 = x[0];
  double sig_scale = par[0];
  double signal = sig_scale * hW2elastic->Interpolate(W2);
  return signal + poly4_fit(x,&par[1]);
}

TH1D *hW2elasticdxanticutresiduals;
double W2totaldxanticutresiduals(double *x, double *par){ //get poly fit to bg with scaled fit to "pure elastics"
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
			else if(key == "coin_mean"){
			coin_mean = val.Atof();
			//cout << "coin mean " << coin_mean << endl;
			}
			else if(key == "coin_sigma"){
			coin_sigma = val.Atof();
			//cout << "coin sigma " << coin_sigma << endl;
			}
			else if(key == "coin_sig_fac"){
			coin_sig_fac = val.Atof();
			//cout << "coin sig fac " << coin_sig_fac << endl;
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
			else if(key == "thetapq_low"){
                        thetapq_low = val.Atof();
                        //cout << "thetapq_low" << thetapq_low << endl;
                        }
			else if(key == "thetapq_high"){
                        thetapq_high = val.Atof();
                        cout << "thetapq_high " << thetapq_high << endl;
                        }
			else if(key == "first_event"){
                        first_event = val.Atoi();
                        cout << "first_event " << first_event << endl;
                        }    
			else if(key == "last_event"){
                        last_event = val.Atoi();
                        cout << "last_event " << last_event << endl;
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

void HCalEfficiency_bootstrap( const char *setup_file_name){
 
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
  double sbs_dist = myData[0].getSBSDist(); //SBS dist
  double ntrack;
  double hcal_offset = getHCalOffset(kin);  
 //BBCal variables
  double BBtr_px[MAXNTRACKS], BBtr_py[MAXNTRACKS], BBtr_pz[MAXNTRACKS], BBtr_p[MAXNTRACKS];
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS],BBtr_chi2[MAXNTRACKS],BBtr_tgth[MAXNTRACKS],BBtr_tgph[MAXNTRACKS];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, nhits, BBeoverp, atime_sh;  


  //HCal variables
  double hcal_clus_atime[maxClusters], hcal_clus_e[maxClusters], hcal_clus_x[maxClusters], hcal_clus_y[maxClusters],tdctime[maxHCalChan], row[maxHCalRows], col[maxHCalCols], cblkid[maxHCalChan], cblke [maxHCalChan];
  UInt_t TBits;
  double xhcal,yhcal,ehcal,nclus,nblk,hcal_tdctime,hcal_idx; 
  double TDCT_id[maxTDCTrigChan], TDCT_tdc[maxTDCTrigChan], hodo_tmean[maxTDCTrigChan];
  int TDCTndata, num_hcal_clusid;
  int hcal_clusid[maxClusters];
  double kineW2,kineQ2,kine_bbthetaprime,kine_nu,kine_omega,kine_phq,kine_thetaq;
 

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
  C->SetBranchStatus( "e.kine.Q2", 1 );
  C->SetBranchStatus( "e.kine.angle", 1 );
  C->SetBranchStatus( "e.kine.nu", 1 );
  C->SetBranchStatus( "e.kine.omega", 1 );
  C->SetBranchStatus( "e.kine.ph_q", 1 );
  C->SetBranchStatus( "e.kine.th_q", 1 );

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus("Ndata.sbs.hcal.clus.id", 1);
  C->SetBranchStatus("sbs.hcal.clus.atime",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("sbs.hcal.clus.x",1);
  C->SetBranchStatus("sbs.hcal.clus.y",1);
  C->SetBranchStatus("sbs.hcal.index",1);
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );
  C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 );

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.vx",vx);
  C->SetBranchAddress("bb.tr.vy",vy);
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
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );

  C->SetBranchAddress( "sbs.hcal.tdctimeblk", &hcal_tdctime );
  C->SetBranchAddress( "sbs.hcal.clus.atime", &hcal_clus_atime);
  C->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &num_hcal_clusid);
  C->SetBranchAddress( "sbs.hcal.clus.e", &hcal_clus_e);
  C->SetBranchAddress( "sbs.hcal.clus.x", &hcal_clus_x);
  C->SetBranchAddress( "sbs.hcal.clus.y", &hcal_clus_y);
  C->SetBranchAddress( "sbs.hcal.index", &hcal_idx);
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
  C->SetBranchAddress( "bb.sh.atimeblk", &atime_sh );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  C->SetBranchAddress( "e.kine.Q2", &kineQ2 );
  C->SetBranchAddress( "e.kine.angle", &kine_bbthetaprime );
  C->SetBranchAddress( "e.kine.omega", &kine_omega );
  C->SetBranchAddress( "e.kine.ph_q", &kine_phq );
  C->SetBranchAddress( "e.kine.th_q", &kine_thetaq );

  C->SetBranchAddress( "bb.gem.track.nhits", &nhits );
  C->SetBranchAddress( "bb.etot_over_p", &BBeoverp );
  
  //Declare output file
  //Setup the name for the output file
  TString outfile = makeOutputFileName_BS(Exp,kin,SBS_field,targ,first_event,last_event);
  TFile *fout = new TFile(outfile,"RECREATE");

 int nentries = C->GetEntries();
 cout << "Total events in TChain: " <<nentries << endl;
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
  double dymax = dyO_p + 2*dxsig_p;

  //3sig dx
  double dxmin3 = dxO_p - 3*dxsig_p;
  double dxmax3 = dxO_p + 3*dxsig_p;
  //cout << "max:" << dxmax3 << " min: " << dxmin3 << endl;

  //Parameters for acceptance matching using MC values. Includes info handle any magnetic field. Not just zero field data. These parameters should work for pass 2 data
  double HCal_right_MC = posHCalYf_MC-HCalblk_l_h_MC+dysig_p;
  double HCal_left_MC = posHCalYi_MC+HCalblk_l_h_MC-dysig_p;
  double HCal_bot_MC = posHCalXi_MC+HCalblk_l_v_MC+dxsig_p-dxO_p;
  double HCal_top_MC = posHCalXf_MC-HCalblk_l_v_MC-dxsig_p-dxO_p;

  //Only valid for pass 1 data
  double HCal_right = posHCalYf-HCalblk_l+dysig_p; 
  double HCal_left = posHCalYi+HCalblk_l-dysig_p;
  double HCal_bot = posHCalXi+HCalblk_l+dxsig_p-dxO_p;
  double HCal_top = posHCalXf-HCalblk_l-dxsig_p-dxO_p;

//Histograms
//HCal detections
  TH2D *hxyexp_nocut = new TH2D("hxyexp_nocut","HCal X Expect vs Y Expect, no cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxyexp_hactivecut = new TH2D("hxyexp_hactivecut","HCal X Expect vs Y Expect, acceptance cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  
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
  TH1D *hdx_cut_MC = new TH1D("hdx_cut_MC","dx, W^{2} and dy cut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hdx_anticut_MC = new TH1D("hdx_anticut_MC","dx, W^{2} and dy anticut, MC param;x_{HCAL}-x_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hdy_cut_MC = new TH1D("hdy_cut_MC","dy, W^{2} and dy cut, MC param;y_{HCAL}-y_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);
  TH1D *hdy_anticut_MC = new TH1D("hdy_anticut_MC","dy, W^{2} and dy anticut, MC param;y_{HCAL}-y_{expect} (m)", hcalbins, hcalfit_l_MC, hcalfit_h_MC);

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

  //Check global cut information
  TH1D *h_ntracks = new TH1D("ntracks","Number of Tracks;", 150, 0, 5);
  TH1D *h_PS_E = new TH1D("h_ps_e"," PS Cluster Energy (GeV);",250,0.0,2.2);
  TH1D *h_PS_E_cut = new TH1D("h_ps_e_cut"," PS Cluster Energy (GeV), Cuts;",250,0.0,2.2);
  TH1D *h_vert_z = new TH1D( "vert_z", "Vertex Position z-direction (m); m", 200, -0.2, 0.2 );
  TH1D *h_vert_z_cut = new TH1D( "vert_z_cut", "Vertex Position z-direction (m), Cuts; m", 200, -0.2, 0.2 );
  TH1D *h_HCal_E = new TH1D( "HCal_E", "HCal Cluster Energy (GeV); GeV", 250, 0, 0.4 );
  TH1D *h_HCal_E_cut = new TH1D( "HCal_E_cut", "HCal Cluster Energy (GeV), Cuts; GeV", 250, 0, 0.4 );
  TH1D *h_HCal_nclus = new TH1D("HCal_nclus","HCal number of clusters meeting threshold;", 250,0,10);
  TH1D *h_HCal_nclus_cut = new TH1D("HCal_nclus_cut","HCal number of clusters meeting threshold, Cuts;", 250,0,10);
  TH1D *h_TPS_SH = new TH1D("h_tps_sh","Total PS and SH cluster energy (GeV);",250,1.5,4.0);
  TH1D *h_TPS_SH_cut = new TH1D("h_tps_sh_cut","Total PS and SH cluster energy (GeV), All cuts;",250,1.5,4.0);
  TH1D *h_nhits = new TH1D("nhits","Number of hits on track;",150, 0, 6);
  TH1D *h_nhits_cut = new TH1D("nhits_cut","Number of hits on track;",150, 0, 6);
  TH1D *h_bbtrp_nocut = new TH1D("bbtrp_nocut","BigBite Track Momentum (GeV), no cut;",300, 0.0, 4.0);
  TH1D *h_bbtrp_cut = new TH1D("bbtrp_cut","BigBite Track Momentum (GeV), cuts;",300, 0.0, 4.0);
  TH1D *h_bbEoverp_nocut = new TH1D("bbEoverp_nocut","BigBite E over p, no cut;",100, 0.0, 2.0);
  TH1D *h_bbEoverp_cut = new TH1D("bbEoverp_cut","BigBite E over p, cuts;",100, 0.0, 2.0);


  //Diagnostic or useful quantities
  TH1D *hE_eloss = new TH1D("E_eloss","Scattered Electron Energy Loss in Target (GeV)",500,0.0, Ebeam *0.001);
  TH1D *hE_ecorr = new TH1D("E_ecorr","Corrected Scattered Electron Energy (GeV)",500,0.95*Ebeam,Ebeam*1.05);
  TH2D *hE_ecorr_vs_vertex = new TH2D("hE_ecorr_vs_vertex",";z_{vertex} (m); E_{corr} (GeV)", 250, -0.125, 0.125, 500, 0.95*Ebeam,1.05*Ebeam);
  TH1D *hE_ep = new TH1D("E_ep","Scattered Electron Energy",500,0.0,Ebeam*1.25);
  hE_ep->GetXaxis()->SetTitle("GeV");
  TH1D *hE_pp = new TH1D("E_pp","Scattered Proton Energy",500, 0.0, Ebeam*1.25);
  hE_pp->GetXaxis()->SetTitle("GeV");
  TH1D *hKE_p = new TH1D( "KE_p", "Scattered Proton Kinetic Energy", 500, 0.0, Ebeam*1.25 );
  hKE_p->GetXaxis()->SetTitle("GeV");
  TH1D *hE_ep_cut = new TH1D("E_ep_cut","Scattered Electron Energy, cut",500,0.0,Ebeam*1.25);
  hE_ep_cut->GetXaxis()->SetTitle("GeV");
  TH1D *hE_pp_cut = new TH1D("E_pp_cut","Scattered Proton Energy, cut",500, 0.0, Ebeam*1.25);
  hE_pp_cut->GetXaxis()->SetTitle("GeV");
  TH1D *hKE_p_cut = new TH1D( "KE_p_cut", "Scattered Proton Kinetic Energy, cut", 500, 0.0, Ebeam*1.25 );
  hKE_p_cut->GetXaxis()->SetTitle("GeV");

  TH1D *hQ2 = new TH1D("Q2","Q2",250,0.5,6.0);
  hQ2->GetXaxis()->SetTitle("GeV");
  TH1D *hQ2_cut = new TH1D("Q2_cut","Q2, cut",250,0.5,6.0);
  hQ2_cut->GetXaxis()->SetTitle("GeV");
  
  TH1D *h_dpel_nocut = new TH1D("h_dpel_nocut","d_pel;p/p_{elastic}(#theta)-1;",250,-1.0,0.5);
  TH1D *h_dpel_BBcut = new TH1D("h_dpel_BBcut","d_pel;p/p_{elastic}(#theta)-1;",250,-1.0,0.5);
  TH1D *h_Ep = new TH1D("h_Ep","Total E over p",400,0.0,2.0);
  h_Ep->GetXaxis()->SetTitle("GeV");
  TH1D *h_p_proton = new TH1D("h_p_proton","Scatter Proton momentum",500,1.0,8.0);
  h_p_proton->GetXaxis()->SetTitle("GeV");
  TH1D *h_p_proton_cut = new TH1D("h_p_proton_cut","Scatter Proton momentum, cut",500,1.0,8.0);
  h_p_proton_cut->GetXaxis()->SetTitle("GeV");
 
  TH1D *h_theta_pq = new TH1D("h_theta_pq","Theta pq for proton",120,0.0,0.6);
  TH1D *h_theta_pq_cut = new TH1D("h_theta_pq_cut","Theta pq for proton, cut",120,0.9*thetapq_low,1.1*thetapq_high);
  TH1D *h_theta_pq_anticut = new TH1D("h_theta_pq_anticut","Theta pq for proton, anticut",120,0.0,0.6);
  
   TH1D *hcoin = new TH1D( "hcoin", "HCal ADCt - BBCal ADCt, no cut; ns", 200, -100, 100 );
   TH1D *hcoin_cut = new TH1D( "hcoin_cut", "HCal ADCt - BBCal ADCt, cuts; ns", 200, -100, 100 );
//Accounting and diagnostic variables
 //long nevent = (long)first_event;
 int treenum = 0, currenttreenum = 0,currentrun =0;
 double progress = 0.0;

 cout << "Processing events.." << endl;
 

//Enable a loop over the entry range we care about instead
  for(long nevent= (long)first_event;nevent <= (long)last_event;nevent++){
  if(nevent%50000==0){
  cout << "nevent " << nevent << " last event " << last_event << endl;
  }
  C->GetEntry(nevent);
  //cout << "The loop " << nevent << endl;
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
    //  cout << "My leaves" << endl;
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }
   bool failedglobal = GlobalCut->EvalInstance(0) == 0;

  
  //Correct the beam energy with energy loss in target using vertex position
  //Convert the length component from meter to cm so that way units properly cancel
  double Eloss = ((vz[0]+(l_tgt/2))*100)*rho_tgt*dEdx_tgt + uwallthick_LH2*rho_Al*dEdx_Al; //aproximately 3 MeV
  hE_eloss->Fill(Eloss);
  double Ecorr = Ebeam-Eloss;
  hE_ecorr->Fill(Ecorr);
  hE_ecorr_vs_vertex->Fill(vz[0],Ecorr);
  //Physics calculations
  double pcorr = BBtr_p[0]-Eloss_outgoing; //neglecting outgoing electron mass
  double p_ep = BBtr_p[0]; // Obtain the magnitude of scattered electron momentum
  h_Ep->Fill(BBeoverp);
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
  h_p_proton->Fill(pp);
   //Define HCal coordinate system
   TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
   TVector3 hcal_xaxis(0,-1,0);
   TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
   TVector3 hcal_origin = hcaldist *hcal_zaxis +hcal_offset*hcal_xaxis;
  //   TVector3 hcalpos = hcal_origin + xhcal * hcal_xaxis + yhcal * hcal_yaxis;
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
   /* TVector3 NeutronDirection = (hcalpos - vertex).Unit();
    double BdL =  maxSBSfield * Dgap * (sbs_field/100); //scale crudely by magnetic field
    double proton_deflection = tan( 0.3 * BdL / qvect_recon.Mag() ) * (hcaldist - (sbs_dist + Dgap/2.0) );
    TVector3 ProtonDirection = (hcalpos + proton_deflection * hcal_xaxis - vertex).Unit();
    double thetapq = acos( ProtonDirection.Dot( qvect_recon.Unit() ) ); 
    h_theta_pq->Fill(thetapq);*/

  //define dx,dy
  double dx = xhcal - xhcal_expect;
  double dy = yhcal - yhcal_expect;
 

  //Fill histos for global cuts
  h_ntracks->Fill(ntrack);
  h_PS_E->Fill(BBps_e);
  h_vert_z->Fill(vz[0]);
  h_HCal_E->Fill(ehcal);
  h_HCal_nclus->Fill(nclus);
  h_TPS_SH->Fill(BBps_e+BBsh_e);
  h_nhits->Fill(nhits);
  h_bbtrp_nocut->Fill(BBtr_p[0]);


  bool inboundsW2 = (W2>0.0) && (W2<W2fitmax);
  //Begin filling output histograms
  
  hxyexp_nocut->Fill(yhcal_expect,xhcal_expect);
 
 ////////////////////////
 //HCAL CLUSTER ANALYSIS
 
 //clone cluster vector for selection analysis
 vector<double> clone_cluster_intime;
 

 //loop through all clusters and select without HCal position information
	//cout << num_hcal_clusid << endl; 
	for(int c = 0; c<num_hcal_clusid; c++){
	 
	 //calculate hcal physics quantities per cluster
	 double atime = hcal_clus_atime[c];
	 double atime_diff = atime - atime_sh; //Assuming best shower time on primary cluster
         double clus_energy = hcal_clus_e[c];

	 //use hcal atime till after pass 2, wide cut around 5 sigma
	 bool passCoin = abs(atime_diff - coin_mean)<coin_sig_fac*coin_sigma;
	 bool passE = clus_energy > hcalemin;

	 //in-time algorithm with new cluster, sort later
	 clone_cluster_intime.push_back(clus_energy);
		if(!passCoin){
		clone_cluster_intime[c]=0;
		}
	 }//end of loop over clusters
 
 //verify that the energy found is the one it should be
 /*if(ehcal != hcal_clus_e[(int)hcal_idx]){
 cout << "Error:Sorting failure. Debug and rerun." << endl;
 cout << "Index: " << (int)hcal_idx << " hcale:" << ehcal << endl;
 }*/
 //Get best intime indices from clone cluster
 int intime_idx = -1;
 double intime = 0.0;
 //sort clusters
 
for(int d=0; d<num_hcal_clusid;d++){

 	if(clone_cluster_intime[d]>intime){
	intime = clone_cluster_intime[d];
	intime_idx=d;
	}
 }

 //For now default intime algorithm for best cluster
 int clus_idx_best = intime_idx;
 
 //calculate important information from best cluster
 double xhcal_bestclus = hcal_clus_x[clus_idx_best];
 double yhcal_bestclus = hcal_clus_y[clus_idx_best];
 double dx_bestclus = xhcal_bestclus - xhcal_expect;
 double dy_bestclus = yhcal_bestclus - yhcal_expect;
 double hcal_atime_bestclus = hcal_clus_atime[clus_idx_best];
 double coin_bestclus = hcal_atime_bestclus - atime_sh;
 double hcal_e_bestclus = hcal_clus_e[clus_idx_best];

 //Need to do calculations here to account for best cluster version
 TVector3 hcalpos = hcal_origin + xhcal_bestclus * hcal_xaxis + yhcal_bestclus * hcal_yaxis;
 //Get theta pq for neutron and proton
 TVector3 NeutronDirection = (hcalpos - vertex).Unit();
 double BdL =  maxSBSfield * Dgap * (sbs_field/100); //scale crudely by magnetic field
 double proton_deflection = tan( 0.3 * BdL / qvect_recon.Mag() ) * (hcaldist - (sbs_dist + Dgap/2.0) );
 TVector3 ProtonDirection = (hcalpos + proton_deflection * hcal_xaxis - vertex).Unit();
 double thetapq = acos( ProtonDirection.Dot( qvect_recon.Unit() ) );
 h_theta_pq->Fill(thetapq);

 hcoin->Fill(coin_bestclus);

 //No events which project hadrons off the face of HCal should be considered. Call this acceptance cut. This inherently relies on the MC DB values. Could also consider not that.
 //This is the acceptance match for BB
 /*bool offhcal =
      yhcal_expect > (HCal_right_MC) ||
      yhcal_expect < (HCal_left_MC) ||
      xhcal_expect > (HCal_top_MC) ||
      xhcal_expect < (HCal_bot_MC);*/

  bool offhcal =
      yhcal_expect > (HCal_right) ||
      yhcal_expect < (HCal_left) ||
      xhcal_expect > (HCal_top) ||
      xhcal_expect < (HCal_bot);

 //Need to make Active Area cut for HCal face use. xhcal and yhcal for diagnostic later
 //Base level cut for all detection efficiency analysis
 if( offhcal ) continue;
 hxyexp_hactivecut->Fill(yhcal_expect,xhcal_expect);

 //all cut bools
 bool gW2 = abs(W2-W2_mean) < W2confac*W2_sigma; //W2 cut around elastic peak . One could implement this in W2_mean and W2_sigma form. But I am not sure this is better
 bool gdy = abs(dy-dyO_p)< confac*dysig_p; //dy cut around elastic peak

 bool gcoin = abs(coin_bestclus - coin_mean) < coin_sigma *coin_sig_fac;
  //Fill histograms without most cuts, we always have an hcal-active-area/acceptance cut
  
  if(inboundsW2){
  hW2_nocut->Fill(W2); //Subtract dy anticut version from this
  }
  
  hW2_nocut_wide->Fill( W2 );
  hdxdy_nocut_MC->Fill( dy_bestclus, dx_bestclus );
  hdx_nocut_MC->Fill( dx_bestclus );  //Subtract W2/dy anticut version from this
  hdy_nocut_MC->Fill( dy_bestclus );
  h_dpel_nocut->Fill(dpel);
  
 

 //anti-cuts
 if(!gW2){
  hdx_badW2cut_MC->Fill(dx_bestclus); //Anticut dx histo
  }
  //having a continue causes a probelm somehow
  //anti-cuts
  if(!gdy ){
  hW2_badYcut->Fill(W2); //Anticut W2 histo
  hdx_badYcut_MC->Fill(dx_bestclus); //Alt anticut dx histo
  }


  //e-arm
  if(gcoin && gW2 && !failedglobal){
  hW2_BBcut->Fill(W2);
  hW2_BBcut_norange->Fill( W2 );
  hdx_BBcut_MC->Fill( dx_bestclus );
  hdy_BBcut_MC->Fill( dy_bestclus );
  hdxdy_BBcut_MC->Fill( dy_bestclus, dx_bestclus );
  h_dpel_BBcut->Fill(dpel);
  }
 
  //h-arm
  if(gdy && !failedglobal && (ehcal>hcalemin) ){
  hW2_dycut->Fill( W2 );

  }
  
  //both arms
  if(gcoin && gW2 && gdy && !failedglobal  && (ehcal>hcalemin)){
  hdx_BBcut_dycut_MC->Fill( dx_bestclus ); //H-arm elastics (efficiency numerator)
  hW2_BBcut_dycut->Fill( W2 );

  }
  
  
  //Make strongest possible hcal-only anticut to obtain background
  if((thetapq>thetapq_high || thetapq<thetapq_low)&&(dx_bestclus<dxmin3||dx_bestclus>dxmax3)&&(dy_bestclus<dymin||dy_bestclus>dymax)&&failedglobal==1 ){
  hW2_anticut->Fill( W2 );
  h_theta_pq_anticut->Fill(thetapq);
  hdx_anticut_MC->Fill(dx_bestclus);
  hdy_anticut_MC->Fill(dy_bestclus);
  }

  //And tightest cut possible
  if( (thetapq<thetapq_high && thetapq>thetapq_low) && (dx_bestclus>dxmin3&&dx_bestclus<dxmax3) && (dy_bestclus>dymin&&dy_bestclus<dymax) && failedglobal==0 ){
  hW2_allcut->Fill( W2 );
  hE_ep_cut->Fill(E_ep);
  hE_pp_cut->Fill(E_pp);
  hKE_p_cut->Fill(KE_p);
  hQ2_cut->Fill(Q2);
  h_p_proton_cut->Fill(pp);
  h_theta_pq_cut->Fill(thetapq);
  hdx_cut_MC->Fill(dx_bestclus);
  hdy_cut_MC->Fill(dy_bestclus);
  
  h_vert_z_cut->Fill(vz[0]);
  h_PS_E_cut->Fill(BBps_e);
  h_HCal_E_cut->Fill(ehcal);
  h_HCal_nclus_cut->Fill(nclus);
  h_TPS_SH_cut->Fill(BBps_e+BBsh_e);
  h_nhits_cut->Fill(nhits);
  h_bbtrp_cut->Fill(BBtr_p[0]);
  hcoin_cut->Fill(coin_bestclus);
  }                              
 //cout << progress<< " " << nevent<< " " <<nentries << endl;
}// This ends the double loop over the events
 
 //dy anti-cut BG subtraction methods
 TH1D *hdx_nocut_MC_clone = (TH1D*)(hdx_nocut_MC->Clone("hdx_nocut_MC_clone")); //for later fit option

 //Make canvas for W2 dy anticut analysis
 TCanvas *c1 = new TCanvas("c1",Form("hcal_%s%smag%i_W2_dyanticut",kin.Data(), targ.Data(), sbs_field),1200,500);
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
 TCanvas *c2 = new TCanvas("c2",Form("hcal_%s%smag%i_dx_W2anticut",kin.Data(), targ.Data(), sbs_field),1200,500);
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
  TCanvas *c0 = new TCanvas("c0",Form("hcal_%s%smag%i_dx_dyanticut",kin.Data(), targ.Data(),sbs_field),1200,500);
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
  TCanvas *c3 = new TCanvas("c3",Form("hcal_%s%smag%i",kin.Data(), targ.Data(), sbs_field),650,500);
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
  TCanvas *c4 = new TCanvas("c4",Form("W2_%s%smag%i",kin.Data(), targ.Data(), sbs_field),1200,500);
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
  // hcaltotalfit->SetParLimits(5,10000,100000);
  // hcaltotalfit->SetParLimits(6,-1.,0.0);
  // hcaltotalfit->SetParLimits(7,0.05,0.5);

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
   TH1D *hW2_bkgdfit = (TH1D*)(hW2_nocut->Clone("hW2_bkgdfit"));
   TH1D *hdx_totfit;
   TH1D *hW2_totfitdxanticutresiduals = (TH1D*)(hW2_anticut->Clone("hW2_totfitdxanticutresiduals"));
    //Make canvas for hcal W2 anticut analysis
    TCanvas *c5 = new TCanvas("c5",Form("W2 dx anticut residuals with extracted BG %s mag%i",kin.Data(), sbs_field),1200,500);
    gStyle->SetPalette(53);
    c5->cd();
    
    //Add background extraction
    hW2elasticdxanticutresiduals = (TH1D*)(hW2_allcut->Clone("hW2elasticdxanticutresiduals"));
    TF1 *totfitdxanticutresiduals = new TF1("totfitdxanticutresiduals",W2totaldxanticutresiduals,0,W2fitmax,6); //tfit npar = 1+pNfit_npar+1
    TF1 *bgdxanticutresiduals = new TF1("bgdxanticutresiduals",poly4_fit,0.,W2fitmax,5);
    totfitdxanticutresiduals->SetLineColor(kGreen);
    //Need TFitResultPtr for error later
    TFitResultPtr W2totfitresidual_ptr = hW2_totfitdxanticutresiduals->Fit("totfitdxanticutresiduals","S");
    hW2_totfitdxanticutresiduals->Fit("totfitdxanticutresiduals","RBM");
    
    double *totpardxanticutresiduals = totfitdxanticutresiduals->GetParameters();
    
    //for(int i =0; i<6 ;i++){
    //cout << "totpardxanticutresiduals " << i << " " << totpardxanticutresiduals[i] << endl;
    //}

    //Get fit parameters for bg function and draw identical function on canvas
    bgdxanticutresiduals->SetParameters(&totpardxanticutresiduals[1]);
    double *bkgddxanticutresiduals = bgdxanticutresiduals->GetParameters();
    
    //for(int j =0; j<5 ;j++){
    //cout << "bkgddxanticutresiduals " << j << " " << bkgddxanticutresiduals[j] << endl;
    //}

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
     //double bgint_dxanticutresiduals = bgdxanticutresiduals->Integral(0.4,1.3)*binfac; //175715
     double bgint_dxanticutresiduals = bgdxanticutresiduals->Integral(0.4,W2fitmax-0.1)*binfac; //175715
     //int W2elasb = 0.0*binfac;
     int W2elasb_dxanticutresiduals = 0.4*binfac;
     //Int_t W2elase = W2fitmax*binfac;
     int W2elase_dxanticutresiduals = (W2fitmax-0.1)*binfac;
     double W2nocutint_dxanticutresiduals = hW2_totfitdxanticutresiduals->Integral(W2elasb_dxanticutresiduals,W2elase_dxanticutresiduals);
     //Double_t tfitint = tfit->Integral(0.0,W2fitmax)*binfac;
     //double totfitint_dxanticutresiduals = totfitdxanticutresiduals->Integral(0.4,1.3)*binfac;
     double totfitint_dxanticutresiduals = totfitdxanticutresiduals->Integral(0.4,W2fitmax-0.1,1.E-4)*binfac;
     double W2elas_dxanticutresiduals = W2nocutint_dxanticutresiduals - bgint_dxanticutresiduals;
     double W2elasfits_dxanticutresiduals = totfitint_dxanticutresiduals - bgint_dxanticutresiduals;
     //cout <<"W2nocutint_dxanticutresiduals: " << W2nocutint_dxanticutresiduals << " totfitint_dxanticutresiduals: " << totfitint_dxanticutresiduals << endl; 
    
    //Add a legend to the canvas
    auto interplegenddxanticutresiduals = new TLegend(0.1,0.6,0.5,0.9);
    interplegenddxanticutresiduals->SetTextSize( 0.03 );
    interplegenddxanticutresiduals->AddEntry( hW2elasticdxanticutresiduals, "Tight Elastic Cut", "l");
    interplegenddxanticutresiduals->AddEntry( bgdxanticutresiduals, "Interpolated Background (scaled)", "l");
    interplegenddxanticutresiduals->AddEntry( totfitdxanticutresiduals, "Total fit (background + signal)", "l");
    interplegenddxanticutresiduals->Draw();
      
    c5->Write();  



    //Make canvas for hcal W2 anticut analysis
    TCanvas *c6 = new TCanvas("c6",Form("W2 with extracted BG %s mag%i",kin.Data(), sbs_field),1200,500);
    gStyle->SetPalette(53);
    c6->cd();
    
    //Add background extraction
    hW2elastic = (TH1D*)(hW2_allcut->Clone("hW2elastic"));
    TF1 *W2_totfit = new TF1("W2_totfit",W2total,0.,W2fitmax,6); //tfit npar = 1+pNfit_npar+1
    TF1 *W2_bkgd = new TF1("W2_bkgd",poly4_fit,0.,W2fitmax,5);
    W2_totfit->SetLineColor(kGreen);
    //Need TFitResultPtr for error later
    TFitResultPtr W2totfit_ptr = hW2_totfit->Fit("W2_totfit","S");
    hW2_totfit->Fit("W2_totfit","RBM");
    double *tpar = W2_totfit->GetParameters();
    //for(int i =0; i<6 ;i++){
    //cout << "tpar " << i << " " << tpar[i] << endl;
    //}
    //Get fit parameters for bg function and draw identical function on canvas
    W2_bkgd->SetParameters(&tpar[1]);
    double *bkgd_par = W2_bkgd->GetParameters();
    //for(int j =0; j<5 ;j++){
    //cout << "bkgd_par " << j << " " << bkgd_par[j] << endl;
    //}
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
   // double W2bgint = W2_bkgd->Integral(0.4,1.3)*binfac; 
     double W2bgint = W2_bkgd->Integral(0.4,W2fitmax-0.1)*binfac; 
    //Int_t W2elasb = 0.*binfac;
    int W2elasb = 0.4*binfac;
    //Int_t W2elase = W2fitmax*binfac;
    int W2elase = (W2fitmax-0.1)*binfac;
    double W2nocutint = hW2_totfit->Integral(W2elasb,W2elase);
    //Double_t tfitint = tfit->Integral(0.0,W2fitmax)*binfac;
    double totfitint = W2_totfit->Integral(0.4,W2fitmax-0.1,1.E-4)*binfac;
    //cout << "W2nocutint: " << W2nocutint << " totfitint: " << totfitint << endl;
    double  W2elas = W2nocutint - W2bgint; 
    double W2elas_fits = totfitint - W2bgint;
    

    //The efficiency that I am using
    //Write eff = (tot-miss)/tot
    double effres = ((W2elas_fits-W2elasfits_dxanticutresiduals) / W2elas_fits)*100.;
    //Write eff = 1 - (miss/tot)
    double effres_simp = (1-(W2elasfits_dxanticutresiduals/W2elas_fits))*100.;

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

  //ACCEPTANCE MATCHING CHECK
  //Make canvas for hcal acceptance match from BB and HCal parameters
  TCanvas *c7 = new TCanvas("c7",Form("HDE accept match %s mag%i",kin.Data(), sbs_field),1200,500);
  gStyle->SetPalette(53);
  c7->cd();
  c7->Divide(2,1);
  c7->cd(1);
  
  //Again only valid for pass 1, need to update for pass 2
  //This is the physical HCal boundary
  double right = posHCalYf;
  double left = posHCalYi;
  double bot = posHCalXi;
  double top = posHCalXf;

  //Make the lines
  TLine* line_Top = new TLine(left,top,right,top);
  line_Top->SetLineColor(kGray);
  TLine* line_Bot = new TLine(left,bot,right,bot);
  line_Bot->SetLineColor(kGray);
  TLine* line_Left = new TLine(left,bot,left,top);
  line_Left->SetLineColor(kGray);
  TLine* line_Right = new TLine(right,bot,right,top);
  line_Top->SetLineColor(kGray);

  //Use the acceptance matching parameters to make lines for the actual acceptance we use
  TLine* line_Top_Acc = new TLine(HCal_left,HCal_top,HCal_right,HCal_top);
  line_Top_Acc->SetLineColor(kGreen);
  TLine* line_Bot_Acc = new TLine(HCal_left,HCal_bot,HCal_right,HCal_bot);
  line_Bot_Acc->SetLineColor(kGreen);
  TLine* line_Left_Acc = new TLine(HCal_left,HCal_bot,HCal_left,HCal_top);
  line_Left_Acc->SetLineColor(kGreen);
  TLine* line_Right_Acc = new TLine(HCal_right,HCal_bot,HCal_right,HCal_top);
  line_Right_Acc->SetLineColor(kGreen);
  
  //Draw histogram with no accept cut
  hxyexp_nocut->Draw("colz");
  
  //Draw the lines overlayed with the Histogram
  line_Top->Draw("SAME");
  line_Bot->Draw("SAME");
  line_Left->Draw("SAME");
  line_Right->Draw("SAME");
  line_Top_Acc->Draw("SAME");
  line_Bot_Acc->Draw("SAME");
  line_Left_Acc->Draw("SAME");
  line_Right_Acc->Draw("SAME");

  auto legend_acc = new TLegend(0.1,0.75,0.3,0.9);
  legend_acc->AddEntry(line_Top, "HCal Boundary", "l");
  legend_acc->AddEntry(line_Top_Acc, "Acceptance Region", "l");
  legend_acc->Draw();
  c7->cd(2);
  
  //Draw histogram with accept cut
  hxyexp_hactivecut->Draw("colz");

  //Draw the lines overlayed with the Histogram
  line_Top->Draw("SAME");
  line_Bot->Draw("SAME");
  line_Left->Draw("SAME");
  line_Right->Draw("SAME");
  line_Top_Acc->Draw("SAME");
  line_Bot_Acc->Draw("SAME");
  line_Left_Acc->Draw("SAME");
  line_Right_Acc->Draw("SAME");

  auto legend_acccut = new TLegend(0.1,0.75,0.3,0.9);
  legend_acccut->AddEntry(line_Top, "HCal Boundary", "l");
  legend_acccut->AddEntry(line_Top_Acc, "Acceptance Region", "l");
  legend_acccut->Draw();

  c7->Write();
     //other method errors, not needed really
     //Binomial
     double eff = hcalelastics/elastics;
     double eff_alt = hcalelastics_antiW2/elastics_alt;
     double eff_alt2 = hcalelastics_antidy/elastics_alt;
     //Double_t eff_althybrid = hcalelastics/elastics_alt;
     double efferr = sqrt(eff*(1-eff)/elastics);
     double efferr_alt = sqrt(eff_alt*(1-eff_alt)/elastics_alt);
     double efferr_alt2 = sqrt(eff_alt2*(1-eff_alt2)/elastics_alt);
     //Double_t efferr_althybrid = sqrt(eff_althybrid*(1-eff_althybrid)/elastics_alt);
     
     //Need to TFitResultPtr to get covariance matrix to make error work. Have the bkgdfit just be subset of that matrix.
     TMatrixD m_totfitresidual = W2totfitresidual_ptr->GetCovarianceMatrix();
     //m_totfitresidual.Print();
     TMatrixD m_bkgdfitresidual = m_totfitresidual.GetSub(1,5,1,5);
     //m_bkgdfitresidual.Print();
     //Get the error for W2elasfits_dxanticutresiduals from the fit integral
     double totfitint_dxanticutresiduals_error = totfitdxanticutresiduals->IntegralError(0.4,W2fitmax-0.1,totpardxanticutresiduals,m_totfitresidual.GetMatrixArray(),1.E-0)*binfac;
    // double totfitint_dxanticutresiduals_error = sqrt(totfitint_dxanticutresiduals);
     double bgint_dxanticutresiduals_error = bgdxanticutresiduals->IntegralError(0.4,W2fitmax-0.1,bkgddxanticutresiduals,m_bkgdfitresidual.GetMatrixArray(),1.E-0)*binfac;
     double W2elas_dxanticutresiduals_error = sqrt(pow(totfitint_dxanticutresiduals_error,2)+pow(bgint_dxanticutresiduals_error,2)) ;
     cout << "totfitint_dxanticutresiduals_error: " << totfitint_dxanticutresiduals_error << " bgint_dxanticutresiduals_error: " << bgint_dxanticutresiduals_error << " W2elas_dxanticutresiduals_error: " << W2elas_dxanticutresiduals_error << endl;
                              
     //Need to TFitResultPtr to get covariance matrix to make error work. Have the bkgdfit just be subset of that matrix.
     TMatrixD m_totfit = W2totfit_ptr->GetCovarianceMatrix();
     //m_totfit.Print();
     TMatrixD m_bkgdfit = m_totfit.GetSub(1,5,1,5);
     //m_bkgdfit.Print();
     //Get the error for W2elas_fits directly from the fit integral
     double totfitint_error = W2_totfit->IntegralError(0.4,W2fitmax-0.1,W2totfit_ptr->GetParams(),m_totfit.GetMatrixArray(),1.E-0)*binfac;
     //double totfitint_error = sqrt(totfitint);
     double W2bgint_error = W2_bkgd->IntegralError(0.4,W2fitmax-0.1,bkgd_par,m_bkgdfit.GetMatrixArray(),1.E-0)*binfac;
     double W2elas_error = sqrt(pow(totfitint_error,2)+pow(W2bgint_error,2));
     //cout << "totfitint_error: " << totfitint_error << " W2bgint_error: "<< W2bgint_error << " W2elas_error: "<< W2elas_error << endl;
   
     //Generalized error using partial derivatives and uncertainties
     //Consider the choice of eff = 1 - (miss/tot)
     // deff/dtot = miss/tot^2
     double partial_eff_tot = W2elasfits_dxanticutresiduals/pow(W2elas_fits,2) ;
     // deff/dmiss = -1/tot 
     double partial_eff_miss = -1.0/W2elas_fits;
     double effres_error_general = (sqrt(pow(partial_eff_tot,2)*pow(W2elas_error,2)+pow(partial_eff_miss,2)*pow(W2elas_dxanticutresiduals_error,2)))*100.;    
     //the top error makes sense, since we are subtractng 2 quantities we need to propagate error this way by taking the square root of the sum of the squares
     //double effres_top_error = sqrt(pow(W2elas_error,2) + pow(W2elas_dxanticutresiduals_error,2));
     //double effres_bot_error = W2elas_error; // This is true since the bottom is just W2elas_fits and we already calculated it
     //double effres_toterror = effres*sqrt((pow(effres_top_error/(W2elas_fits-W2elasfits_dxanticutresiduals),2)+pow(effres_bot_error/W2elas_fits,2)));
     double effres_error_binomial = (sqrt((effres/100.0)*(1.0-(effres/100.0))/W2elas_fits))*100.0;
     //cout << "effres_top_error: " << effres_top_error << " effres_bot_error: " << effres_bot_error << " effres_toterror: " << effres_toterror <<" error_binomial: "<< effres_error_binomial << endl;
     //cout << "partial_eff_tot: " << partial_eff_tot << " partial_eff_miss: " << partial_eff_miss << endl;

   TString reportfile = makeReportFileName_BS(Exp,kin,SBS_field,targ,first_event,last_event);   


    
   //Declare outfile
   ofstream report;
   report.open( reportfile );   

   cout << "Chi^2 for total fit: " << W2_totalfit->GetChisquare() << endl;
   report << "Chi^2 for total fit: " << W2_totalfit->GetChisquare() << endl;
   
   cout << "Number of degrees of freedom for total fit: " << W2_totalfit->GetNDF() << endl;
   report << "Number of degrees of freedom for total fit: " << W2_totalfit->GetNDF() << endl;

   report << "HCal detection efficiency report for " << kin << " LH2 data at " << sbs_field << " percent field" << endl << endl;
   

   //other inaccurate methods, might bring one or 2 back if I figure it out
   /*cout << "Total elastics detected in HCal: " << hcalelastics << endl;
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
   report << "Total entries in W2 plot: " << earmEvents << endl << endl; */
   
   cout << "Total elastics from W2 Tight elastic cut (global, dx, dy, thetapq) " << W2elas_fits << endl;
   report << "Total elastics from W2 Tight elastic cut (global, dx, dy, thetapq) " << W2elas_fits << endl;

    cout << "Total elastics from W2 anti Tight elastic cut (global, dx, dy, thetapq) " << W2elasfits_dxanticutresiduals << endl;
   report << "Total elastics from W2 anti Tight elastic cut (global, dx, dy, thetapq) " << W2elasfits_dxanticutresiduals << endl;
  
   /*cout << "Total error from adding in quadrature (as percentage) " << effres_toterror << endl;
   report << "Total error from adding in quadrature (as percentage) " << effres_toterror << endl;*/

   cout << "Total error from binomial error (as percentage) " << effres_error_binomial << endl;
   report << "Total error from binomial error (as percentage) " << effres_error_binomial << endl;

   cout << "Total error from quadrature sum (as percentage) " << effres_error_general << endl;
   report << "Total error from quadrature sum (as percentage) " << effres_error_general << endl;

  /*
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
  */

  cout << "===================================================================" << endl;
  cout << "Total HCal elastic detection efficiency  Tight elastic cut (global, dx, dy, thetapq): " << effres << " +/- " << "" << endl;
  cout << "===================================================================" << endl << endl;

  report << "===================================================================" << endl;
  report << "Total HCal elastic detection efficiency  Tight elastic cut (global, dx, dy, thetapq): " << effres << " +/- " << "" << endl;
  report << "===================================================================" << endl << endl;


  cout << "===================================================================" << endl;
  cout << "Total HCal elastic detection efficiency  Tight elastic cut (global, dx, dy, thetapq) simple ver: " << effres_simp << " +/- " << "" << endl;
  cout << "===================================================================" << endl << endl;

  report << "===================================================================" << endl;
  report << "Total HCal elastic detection efficiency  Tight elastic cut (global, dx, dy, thetapq) simple ver: " << effres_simp << " +/- " << "" << endl;
  report << "===================================================================" << endl << endl;

 // cout << endl << endl << "raw bincontent efficiency anticut dy: " << hely_bincontents / eel_bincontents << endl;
 // cout << endl << endl << "raw bincontent efficiency anticut W2: " << hel_bincontents / eel_bincontents << endl;

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
  c6->Print(plotname.Data(),"pdf");
  c7->Print(end.Data(),"pdf");
  fout->Write();
//TH1D *hdx_residual = (TH1D*)hdx_cut->Clone("hdx_residual");
//hdx_residual->SetTitle("HCal dx residual (data - fit),All cuts");
//hdx_residual->GetXaxis()->SetTitle("m");
  





  watch->Stop();

  // Send time efficiency report to console
   cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;


} 
