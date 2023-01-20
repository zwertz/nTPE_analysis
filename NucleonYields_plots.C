//Still need to implement fiducial cuts and fitting functions for background, protons, and neutrons. Also need to implement basic blinding via random number and yields




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
//Maybe some of these can be elevated to the setup config file
//Some other parameters that are potentially useful but need to understand them. Taken from Seeds Script
const int kNcell = 288; // Total number of HCal modules
const int kNrows = 24; // Total number of HCal rows
const int kNcols = 12; // Total number of HCal columns
const int kNtdc = 1000; // Reasonable max number of tdc signals per event
const double Xi = -2.20; // Distance from beam center to top of HCal in m
const double Xf = 1.47; // Distance from beam center to bottom of HCal in m
const double Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const double Yf = 0.853; // Distance from beam center to beam side of HCal in m
const double M_p = 0.938272; // Mass proton GeV
const double M_n = 0.939565; // Mass neutron GeV
const double sampfrac = 0.077; // Most recent estimate of the sampling fraction via MC

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double celldiameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch  

//Define max number of tracks per event
int MAXNTRACKS;

double tdiffmax; // Max deviation from coin via tdctrig cut



TCut globalcut = "";
TString Exp,kin,data_file_name,kinematic_file_name,targ;
int SBS_field,useAlshield;
double W_mean,W_sigma,tdiff,dxO_n,dyO_n,dxsig_n,dysig_n,dxO_p,dyO_p,dxsig_p,dysig_p,dxmax;
//double W_min = W_mean - W_sigma;
//double W_max = W_mean + W_sigma;
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
			else if(key == "W_mean"){
                	W_mean = val.Atof();
                        //cout << "W Mean " << W_mean << endl;
			}
			else if(key == "W_sigma"){
			W_sigma = val.Atof();
                        //cout << "W Sigma " << W_sigma << endl;
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
			else if(key == "dxO_n"){
                        dxO_n = val.Atof();
                        //cout << "x-position of neutron spot" << dxO_n << endl;
                        }
			else if(key == "dyO_n"){
                        dyO_n = val.Atof();
                        //cout << "y-position of neutron spot" << dyO_n << endl;
                        }
			else if(key == "dxsig_n"){
                        dxsig_n = val.Atof();
                        //cout << "x sigma of neutron spot" << dxsig_n << endl;
                        }
			else if(key == "dysig_n"){
                        dysig_n = val.Atof();
                        //cout << "y sigma of neutron spot" << dysig_n << endl;
                        }
			else if(key == "dxO_n"){
                        dxO_n = val.Atof();
                        //cout << "x-position of neutron spot" << dxO_n << endl;
                        }
                        else if(key == "dyO_n"){
                        dyO_n = val.Atof();
                        //cout << "y-position of neutron spot" << dyO_n << endl;
                        }
                        else if(key == "dxsig_n"){
                        dxsig_n = val.Atof();
                        //cout << "x sigma of neutron spot" << dxsig_n << endl;
                        }                       
                        else if(key == "dysig_n"){
                        dysig_n = val.Atof();
                        //cout << "y sigma of neutron spot" << dysig_n << endl;
                        }
			else if(key == "useAlshield"){
                        useAlshield = val.Atoi();
                        //cout << "Use Al shield" << useAlshield << endl;
                        }
			else{
			//We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
			cout << "Found a key that this script can't handle. Fix that!" << endl;
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

void NucleonYields_plots( const char *setup_file_name){
 //Constructor for the data structure we will use to make plots
 TChain *C = new TChain("T");

 //parse function to get in the information that The One Config file has and is manipulated
 //this function will inialize the global parameters: globalcut,Exp,kin,data_file_name,kinematic_file_name,SBS_field,W_mean,W_sigma,targ,runnums
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
  double ntrack;
  double BBtr_px[MAXNTRACKS], BBtr_py[MAXNTRACKS], BBtr_pz[MAXNTRACKS], BBtr_p[MAXNTRACKS];
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS];
  double xhcal,yhcal,ehcal,nblk,nclus,SHnclus,PSnclus;
  //Not sure exactly what these are usful for. Maybe describing energy in diff HCal blocks? Might not need all of these depending on the Histograms
  double atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke [kNcell];
  UInt_t TBits;
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
  double TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc];
  int TDCTndata;
 

  //Setup leaves
  C->SetMakeClass(1); //Allows for viewing of general detector params from Event_Branch
  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);

  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "fEvtHdr.fTrigBits", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "g.trigbits", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p",1);


  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.tr.px",BBtr_px);
  C->SetBranchAddress("bb.tr.py",BBtr_py);
  C->SetBranchAddress("bb.tr.pz",BBtr_pz);
  C->SetBranchAddress("bb.tr.p",BBtr_p);
  C->SetBranchAddress("sbs.hcal.x",&xhcal);
  C->SetBranchAddress("sbs.hcal.y",&yhcal);
  C->SetBranchAddress("sbs.hcal.e",&ehcal);

  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", col );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.sh.nclus", &SHnclus );
  C->SetBranchAddress( "bb.ps.nclus", &PSnclus );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits ); 
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  
  
  //Declare output file
  //Setup the name for the output file
  TString outfile = makeOutputFileName(Exp,kin,SBS_field,targ);
  TFile *fout = new TFile(outfile,"RECREATE");

   //Cut on global parameters from setup config
   TEventList *elist = new TEventList("elist","");
   C->Draw(">>elist",globalcut);

  double dx_low = -3;
  double dx_high = 3;
  double dy_low = -3;
  double dy_high = 3;
  // Initialize histograms
 // TH1D *h_atime = new TH1D( "atime", "HCal ADC Time, All Channels; ns", 160, 0, 160 );
 // TH2D *h_CvCh = new TH2D( "CvCh", "HCal Coeff Single Block Clusters; channel, GeV", kNcell, 0, kNcell, 200, 0, 1.0 );
  TH1D *h_E_all = new TH1D( "E_all", "HCal Cluster Energy (GeV), All Channels; GeV", 250, 0, 0.4 );
  TH1D *h_E_cut = new TH1D( "E_cut", "HCal Cluster Energy (GeV) All Cuts, All Channels; GeV", 250, 0, 0.4 );
 // TH1D *h_E_exp = new TH1D( "E_exp", "Expected Energy Dep in HCal; GeV", 100, 0, 0.2 );
  TH1D *h_vert = new TH1D( "vert", "Vertex Position (m); m", 200, -0.4, 0.4 );
 // TH2D *h_EvCh = new TH2D( "EvCh", "HCal Cluster E Single Block Clusters; channel, GeV", kNcell, 0, kNcell, 50, 0, 0.5 );
  TH1D *h_W2 = new TH1D( "W2", "W2 (GeV) No Cuts; GeV", 250, -1.0, 4.0 );
  TH1D *h_W2recon = new TH1D( "W2recon", "W2 Reconstructed (GeV) No Cuts; GeV", 250, -1.0, 4.0 );
  TH2D *hdxdy_cut = new TH2D("dxdy_cut","HCal dxdy All Cuts, All Channels;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, dx_low, dx_high );
  TH2D *hdxdy_all = new TH2D("dxdy_all","HCal dxdy All Channels ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",250,-1.25,1.25,250,dx_low,dx_high);
  TH1D *hdx = new TH1D( "dx", "HCal dx (m); m", 200, dx_low, dx_high );
  TH1D *hdy = new TH1D( "dy", "HCal dy (m); m", 200, dy_low, dy_high );
  TH1D *hKE_p = new TH1D( "KE_p", "Scattered Proton Kinetic Energy", 500, 0.0, 5.0 );
  TH2D *hdxVE = new TH2D("dxVE","dx vs VE;x_{HCAL}-x_{expect} (m); GeV", 100, -4.0, 2.0, 100, 0, 0.5 );
  TH1D *hKElow = new TH1D( "KElow", "Lowest Elastic E Sampled in HCal (GeV)", 500, 0.0, 0.2 );
  TH1D *htimeDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH2D *hrowcol = new TH2D( "hrowcol", "HCal Block Position Elastics, HCal; Col; Row", kNcols, 0, kNcols, kNrows, -kNrows, 0 );
  TH1D *hX = new TH1D( "X", "HCal X (m); m", 100, dx_low, dx_high );
  TH1D *hX_expect = new TH1D( "X_expect", "HCal X Expect (m); m", 100, dx_low, dx_high );
  TH1D *hY = new TH1D( "Y", "HCal Y (m); m", 100, dy_low, dy_high );
  TH1D *hY_expect = new TH1D( "Y_expect", "HCal Y Expect (m); m", 100, dy_low, dy_high );
  TH1D *hdx_cut = new TH1D( "dx_cut", "HCal dx (m), All cuts; m", 200, dx_low, dx_high );
  TH1D *hdy_cut = new TH1D( "dy_cut", "HCal dy (m), All cuts; m", 200, dy_low, dy_high );
  TH1D *hX_cut = new TH1D( "X_cut", "HCal X (m), All cuts; m", 100, dx_low, dx_high );
  TH1D *hX_expect_cut = new TH1D( "X_expect_cut", "HCal X Expect (m), All cuts; m", 100, dx_low, dx_high );
  TH1D *hY_cut = new TH1D( "Y_cut", "HCal Y (m),All cuts; m", 100, dy_low, dy_high );
  TH1D *hY_expect_cut = new TH1D( "Y_expect_cut", "HCal Y Expect (m), All cuts; m", 100, dy_low, dy_high );
  TH1D *h_W2recon_cut = new TH1D( "W2recon_cut", "W2 Reconstructed (GeV) with cuts; GeV", 250, -1.0, 4.0 );
  TH1D *htimeDiff_cut = new TH1D( "hDiff_cut","HCal time - BBCal time (ns), All Cuts", 150, 450, 600 );
  TH1D *h_Wrecon = new TH1D( "Wrecon", "W Reconstructed (GeV) no cuts; GeV", 250, -0.5, 3.0 );
  TH1D *h_Wrecon_cut = new TH1D( "Wrecon_cuts", "W Reconstructed (GeV) With cuts; GeV", 250, -0.5, 3.0 );
  TH1D *h_dpel = new TH1D("h_dpel","d_pel;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *h_TPS_SH = new TH1D("h_tps_sh","Total PS and SH cluster energy (GeV);",250,1.5,2.8);
  TH1D *h_PS_E = new TH1D("h_ps_e"," PS Cluster Energy (GeV);",250,0.0,2.2); 
  TH1D *h_dpel_cut = new TH1D("h_dpel_cut","d_pel,All Cuts;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *h_TPS_SH_cut = new TH1D("h_tps_sh_cut","Total PS and SH cluster energy (GeV), All cuts;",250,1.5,2.8);
  TH1D *h_PS_E_cut = new TH1D("h_ps_e_cut"," PS Cluster Energy (GeV), All cuts;",250,0.0,2.2);
  TH1D *hdr = new TH1D("dr","HCal dr (m), No Cuts; m", 200, -3, 3 );
  TH1D *hdr_cut = new TH1D("dr_cut","HCal dr (m), All Cuts; m", 200, -3, 3 );
  TH2D *hdxdy_pcut = new TH2D("hdxdy_pcut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",12,-0.9,0.9,24,-2.165,1.435);
  TH2D *hdxdy_ncut = new TH2D("hdxdy_ncut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",12,-0.9,0.9,24,-2.165,1.435);
  TH2D *hdxdy_fcut = new TH2D("hdxdy_fcut",";y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",250,-1.25,1.25,250,dx_low,dx_high); 
  TH1D *hdx_fcut = new TH1D( "dx_fcut","; x_{HCAL}-x_{expect} (m)", 200, dx_low, dx_high );
  TH1D *hdy_fcut = new TH1D( "dy_fcut","; y_{HCAL}-y_{expect} (m)", 200, dy_low, dy_high );
  TH1D *h_W2recon_fcut = new TH1D( "W2recon_fcut", "W2recon_fcut; GeV", 250, -1.0, 4.0 );
  TH1D *h_Wrecon_fcut = new TH1D( "Wrecon_fcut", "Wrecon_fcut; GeV", 250, -1.0, 4.0 );
  TH2D *hxy = new TH2D("hxy",";y_{HCAL} (m); x_{HCAL} (m)",12,-0.9,0.9,24,-2.165,1.435);
  TH2D *hxy_fcut = new TH2D("hxy_cut",";y_{HCAL} (m); x_{HCAL} (m)",12,-0.9,0.9,24,-2.165,1.435);
  TH2D *hxy_pcut = new TH2D("hxy_pcut",";y_{HCAL} (m); x_{HCAL} (m)",12,-0.9,0.9,24,-2.165,1.435);
  TH2D *hxy_ncut = new TH2D("hxy",";y_{HCAL} (m); x_{HCAL} (m)",12,-0.9,0.9,24,-2.165,1.435);

 //variables for script
  long nevents = elist -> GetN();
  int elastic_yield = 0;
  double pBeam = Ebeam/(1.0 +(Ebeam/M_p)*(1.0-cos(bbtheta)));

 //mean energy loss of the beam before scattering
 double Eloss_outgoing = (celldiameter/2)/sin(bbtheta)*rho_tgt*dEdx_tgt; //Approximately 1 MeV, one could correct more using the raster position
 if(useAlshield !=0){
 Eloss_outgoing += Alshieldthick * rho_Al*dEdx_Al;
 }

  
  for( int r =0; r<kNrows; r++){
    for (int c=0; c<kNcols; c++){
      hrowcol->Fill( (c+1), -(r+1) );
    }
  }
  
  //Correct the beam energy with energy loss in target using vertex position
  double Eloss = (vz[0]+(l_tgt/2))*rho_tgt*dEdx_tgt + uwallthick_LH2*rho_Al*dEdx_Al; //aproximately 3 MeV
  //fill Eloss in histogram

  double Ecorr = Ebeam-Eloss;
  //fill histo

  double pcorr = BBtr_p[0]-Eloss_outgoing; //neglecting outgoing electron mass

 


  //Beam and detector vars including definition of HCal coordinate axis
  TLorentzVector Pbeam(0,0,Ecorr,Ecorr);
  TLorentzVector Ptarg( 0, 0, 0, 0.5*(M_p+M_n) ); //Average of proton and neutron rest mass
  TVector3 hcal_origin( -hcaldist*sin(hcaltheta), 0, hcaldist*cos(hcaltheta) );
  //TVector3 hcal_origin( -hcaldist*sin(sbstheta), 0, hcaldist*cos(sbstheta) );
  TVector3 hcal_zaxis = hcal_origin.Unit();
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();


  double W_min = W_mean - W_sigma;
  double W_max = W_mean + W_sigma;
  //create loop with the number of events from elist
  for(long nevent = 0;nevent < nevents; nevent++){
  
  if( nevent%10000 == 0 ){
  cout << "Loop: " << nevent << "/" << nevents << ". Elastic yield: " << elastic_yield << endl;
  }
  cout.flush();
 
  //Get sequential event from elist (parsed chain with cuts from setup file)
  C->GetEntry( elist->GetEntry( nevent ) );
  
  //Define useful vectors
  TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
  TLorentzVector q = Pbeam - kprime; //Standard q-vector
  TVector3 vertex( 0, 0, vz[0] ); //Location of scattering event in target
  TVector3 qdir = q.Vect().Unit(); //q-vector direction
  
  double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qdir.Dot( hcal_zaxis );
  TVector3 hcal_intersect = vertex + sintersect * qdir;
  double xhcal_expect = hcal_intersect.Dot( hcal_xaxis );
  double yhcal_expect = hcal_intersect.Dot( hcal_yaxis );
  //cout << yhcal_expect << endl;
  double W2recon = (Ptarg + q).M2();//M2 is magnitude squared of the 4-vector
  double Wrecon;
  double E_ep = BBtr_p[0]; // Obtain the scattered electron energy, neglect mass e
  double p_ep = BBtr_p[0]; // Obtain the magnitude of scattered electron momentum
  double etheta = acos( BBtr_pz[0]/BBtr_p[0] );
  double Q2 = 2*Ecorr*E_ep*( 1-cos(etheta) );
  double nu = Ecorr-E_ep; // Obtain energy transfer
  double W2 = pow( M_p,2 )+2*M_p*nu-Q2; // Obtain W2 from Q2 and nu
  double ephi = atan2( BBtr_py[0], BBtr_px[0] );
  double phinucleon = ephi + TMath::Pi(); //assume coplanarity
  double thetanucleon = acos( (Ecorr - BBtr_pz[0])/p_ep ); //use elastic constraint on nucleon kinematics
  TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
  double KE_p = nu; //For elastics
  if(W2recon >0){
  Wrecon = sqrt(W2);
  
  }
  double pelastic = Ecorr /(1.0+(Ecorr/M_p)*(1.0-cos(etheta)));
  double dpel = BBtr_p[0]/pelastic - 1.0;
  hKE_p->Fill( KE_p );

  //define dx,dy, and dr 
  double dx = xhcal - xhcal_expect;
  double dy = yhcal - yhcal_expect;
  double dr;
  double r = sqrt(pow(dx,2)+ pow(dy,2)); 
  if(r != 0){
  //handle the case where the radius does not equal zero
  double myang = TMath::ATan2(dx,dy);
  double myang_deg = RadToDeg(myang);
  //cout << myang_deg << endl; 
  	if((myang_deg >= -180.00) && (myang_deg < 0.00) ){
  	dr = (-1)*r;
  	}else{
 	 dr = r;  
  	}
 }else{
 //handle if somehow r = 0
 dr = r;
 //cout << dr << endl;
 }
  //cout << dr << " " << myang_deg << endl;

  //Cut on BBCal and HCal trigger coincidence. 
  double bbcal_time=0., hcal_time=0., coin_time=0., rf_time=0.;
    bool cointrig = false;
    for(int ihit=0; ihit<TDCTndata; ihit++){
      if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==1) {
        coin_time=TDCT_tdc[ihit];
        cointrig=true;
      }
      if(TDCT_id[ihit]==4) rf_time=TDCT_tdc[ihit];
    }
 
  double diff = hcal_time - bbcal_time;
  //Fill some histograms
  //cout << yhcal - yhcal_expect << endl;
  htimeDiff->Fill( diff );
  hdxdy_all->Fill( dy, dx);
  h_vert->Fill(vz[0]);
  hxy->Fill(yhcal,xhcal);
  h_W2->Fill(W2);
  h_W2recon->Fill(W2recon);
  h_E_all->Fill(ehcal);
  h_Wrecon->Fill(Wrecon);
  
  hX_expect->Fill( xhcal_expect);
  hY_expect->Fill( yhcal_expect);
  hdx->Fill(dx );
  hdy->Fill(dy );
  hX->Fill( xhcal);
  hY->Fill( yhcal );
  h_dpel->Fill(dpel);
  h_TPS_SH->Fill(BBps_e+BBsh_e);
  h_PS_E->Fill(BBps_e);

  hdr->Fill(dr);
  
  ///////////////
  //    //PRIMARYCUTS//
  //Cut on vertex inside target. This cut is already handled by global cut in some capacity
 // if( abs(BBtr_vz[0])>0.27 ) continue;
  //Cut on coin
  //if( TBits==1 ) continue; //TBits 1 is BBCal trig, TBits 5 is coin
  //Atime cut hcal
  //if( atime[0]>60. && atime[0]<90. ) //Observed HCal ADC time peak
  //Cut on  W2
  if( abs(Wrecon-W_mean)>W_sigma ) continue; //Observed W2 peak
  //if( Wrecon >= W_min && Wrecon <= W_max) continue; //Observed W2 peak
  //Cut on BBCal HCal trig diff
  if( abs(diff-tdiff)>tdiffmax ) continue; //BBCal/HCal trigger difference time
  //Bigbite track / HCal angular correlation cut
  //ENDCUTS//
  ///////////
  //elastic_yield++;
  //Fill some histograms
  hKElow->Fill( KE_p*sampfrac );
  hdxVE->Fill(dx,ehcal);
  //Fill delta plots and others
  hdxdy_cut->Fill( dy, dx );
  h_E_cut->Fill(ehcal);
  hdx_cut->Fill( dx );
  hdy_cut->Fill( dy );
  //h_atime->Fill( atime[0] );
  hX_cut->Fill( xhcal);
  hY_cut->Fill( yhcal );
  h_W2recon_cut->Fill(W2recon);
  hX_expect_cut->Fill( xhcal_expect);
  hY_expect_cut->Fill( yhcal_expect);
  htimeDiff_cut->Fill( diff );
  h_Wrecon_cut->Fill(Wrecon);
  h_dpel_cut->Fill(dpel);
  h_TPS_SH_cut->Fill(BBps_e+BBsh_e);
  h_PS_E_cut->Fill(BBps_e);

  hdr_cut->Fill(dr);


  ///////////////
  //SECONDARYCUTS
  //Reject events where the primary block in the primary cluster is on the edge of the acceptance
 // if( row[0]==0 ||
 // row[0]==23 ||
 // col[0]==0 ||
 // col[0]==11 )
 //  continue;
 // elastic_yield++;
  //ENDCUTS
  
  
 //Probably the place for acceptance/fiducial cuts
 

 bool HCal_on = false, find_p = false, find_n = false;
 //Are these physical coordinates of the HCal in meters?
 if(yhcal>-0.75 && yhcal<0.75 && xhcal>-2.015 && xhcal<1.285){
 HCal_on = true;
 }
 //equation for a circle around proton spot. Maybe make ellipse later
 if(pow((dx-dxO_p)/dxsig_p,2)+pow((dy-dyO_p)/dysig_p,2)<= pow(3.5,2)){
 find_p = true;
 }
 //equation for a circle around neutron spot. Maybe make ellipse later
 //
 if(pow((dx-dxO_n)/dxsig_n,2)+pow((dy-dyO_n)/dysig_n,2)<= pow(2.5,2)){
 find_n = true;
 }
 //Fill 2D histograms if we find protons or neutrons 
 if(HCal_on == true && find_p == true){
 hdxdy_pcut->Fill(dy,dx);
 }
 
 if(HCal_on == true && find_n == true){
 hdxdy_ncut->Fill(dy,dx);
 }

 
 //Fill some histograms if we find protons or neutrons and they are within a range in the vertical direction
 if(HCal_on == true && find_n == true){
 	if((xhcal-dxmax)> -2.015 ){
	//need to understand why that value
	hdxdy_fcut->Fill(dy,dx);
        hdx_fcut->Fill(dx);
	hdy_fcut->Fill(dy);
	h_W2recon_fcut->Fill(W2recon);
	h_Wrecon_fcut->Fill(Wrecon);
	hxy_fcut->Fill(yhcal,xhcal);
	hxy_ncut->Fill(yhcal,xhcal);
        elastic_yield++;
	}
 
 }else if (HCal_on == true && find_p == true){
 	if((xhcal + dxmax)<1.285){
	//again need to know why the above number
	hdxdy_fcut->Fill(dy,dx);
        hdx_fcut->Fill(dx);
        hdy_fcut->Fill(dy);
        h_W2recon_fcut->Fill(W2recon);
        h_Wrecon_fcut->Fill(Wrecon);
        hxy_fcut->Fill(yhcal,xhcal);
        hxy_pcut->Fill(yhcal,xhcal);	
	elastic_yield++;
	}
 }
}
 //Declare Canvas
 TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
 c1->Divide(2,1);
 //place holder till I do fiducials
 c1->cd(1);
 hdxdy_cut->Draw("colz"); 

 //Sets to no fitting?
 gStyle->SetOptFit(0);
 //Set to no stat box?
 gStyle->SetOptStat(0);

 //Working on fitting the dx plot
 c1->cd(2);
 //make a clone of the dx plot
 TH1D *hdxcut_clone = (TH1D*)hdx_cut->Clone("hdxcut_clone");
 
 //Initialize fit parameters
 
 vector<Double_t> myParam (11);
 myParam= fit_Params(kin);
 Double_t fit_low = myParam[9];
 Double_t fit_high = myParam[10];


 //Make the total fit
 c1->cd(2);
 TF1 *totalFit = new TF1("totalfit",Tot_fit,fit_low,fit_high,9);
 totalFit->SetLineWidth(4);
 totalFit->SetLineColor(kMagenta);
 //Set or Fix parameter. Right now there are 9
 
 totalFit->SetParameter(0,myParam[0]);
 totalFit->SetParLimits(0,150,900);
 totalFit->SetParameter(1,myParam[1]);
 totalFit->SetParameter(2,myParam[2]);
 totalFit->SetParameter(3,myParam[3]);
 totalFit->FixParameter(4,myParam[4]); //fix the mean rather than fitting radiative tails
 totalFit->SetParameter(5,myParam[5]);
 totalFit->SetParameter(6,myParam[6]);
 totalFit->FixParameter(7,myParam[7]); //fix the mean
 totalFit->SetParameter(8,myParam[8]);


 //Q = minimum printing, R =fit using fitting range specified in the function range, B= Used when fixing or setting limits on one or more parameters in predefined function, + = add to list of fitted functions
 hdxcut_clone->Fit("totalfit","QRB+");
 TF1 *bkgd = new TF1("bkgd",BG_fit,fit_low,fit_high,3);
 TF1 *proton = new TF1("proton",P_fit,fit_low,fit_high,3);
 TF1 *neutron = new TF1("neutron",N_fit,fit_low,fit_high,3);
 bkgd->SetLineColor(kBlack);
 proton->SetLineColor(kRed);
 neutron->SetLineColor(kBlue);
 
 bkgd->SetParameters(&myParam[0]);
 proton->SetParameters(&myParam[3]);
 neutron->SetParameters(&myParam[6]);
 
 
 bkgd->Draw("same");
 proton->Draw("same");
 neutron->Draw("same");

 //Generate yields
 double p_yield = proton->Integral(dx_low,dx_high)/hdxcut_clone->GetBinWidth(1);
 double n_yield = neutron->Integral(dx_low,dx_high)/hdxcut_clone->GetBinWidth(1);

 //blind yields
 blind_factor *p_blind = new blind_factor("GetOffMyLawnYou");
 blind_factor *n_blind = new blind_factor("ILaughedTheToaster");
 double p_tot = p_blind->blind(p_yield);
 double n_tot = n_blind->blind(n_yield);

// cout << p_blind->getBlindFac() << " " << n_blind->getBlindFac() << endl;

 //make legend
 auto legend = new TLegend(0.5,0.7,0.9,0.9);
 legend->SetTextSize(0.03);
 legend->AddEntry(bkgd,"Background Fit","l");
 legend->AddEntry(proton,Form("Proton Fit, Yield: %d",(int)p_tot),"l");
 legend->AddEntry(neutron,Form("Neutron Fit, Yield: %d",(int)n_tot),"l");
 legend->AddEntry(totalFit,"Total Fit","l");
 legend->Draw();

 c1->Write();

 

  //output the info to some files
  cout << "Total calibration yield for run with current cuts: " << elastic_yield << "." << endl;
  c1->Write();
  TString plotname = outfile;
  plotname.ReplaceAll(".root",".pdf");
  c1->Print(plotname.Data(),"pdf");
  plotname.ReplaceAll(".pdf",".png");
  c1->Print(plotname.Data(),"png");
  fout->Write();


  cout << "Histograms populated and written to file." << endl;



} 
