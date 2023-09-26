//Author Ezekiel Wertz
//Includes Acceptance Cut based on MC parameters
//Includes acceptance matching and fiducial cuts
//Includes W2 cut
//Includes dy cut
//Makes ellipses around protons and neutrons in lieu of a thetapq cut
//Includes yield blinding
//Basic yield fitting using a 4th order poly background and gaussians for the proton and neutron data


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

double tdiffmax,dx_low,dx_high,dy_low,dy_high,dxsig_n_fac,dxsig_p_fac,dysig_fac; 
TString Exp,kin,data_file_name,kinematic_file_name,targ;
int SBS_field,useAlshield;
double W2_low,W2_high,tdiff,dxO_n,dyO_n,dxsig_n,dysig_n,dxO_p,dyO_p,dxsig_p,dysig_p,dxmax;

//Define some fits to be used for analysis
//Background fit of a fourth-order polynomial
bool reject_bkgd;

double BG_4fit(double *x, double *param){

if((kin == "SBS4") && (SBS_field == 30) && reject_bkgd && x[0]> -2.1  && x[0] < 1.15 ){
TF1::RejectPoint();
return 0;
}

if((kin == "SBS4") && (SBS_field == 50) && reject_bkgd && x[0]> -1.5  && x[0] < 0.7 ){
TF1::RejectPoint();
return 0;
}

if((kin == "SBS8") && (SBS_field == 100) && reject_bkgd && x[0]> -2.3  && x[0] < 0.7 ){
TF1::RejectPoint();
return 0;
}

if((kin == "SBS8") && (SBS_field == 50) && reject_bkgd && x[0]> -1.5  && x[0] < 0.7 ){
 TF1::RejectPoint();
 return 0;
}

if((kin == "SBS8") && (SBS_field == 70) && reject_bkgd && x[0]> -1.6  && x[0] < 0.75 ){
TF1::RejectPoint();
return 0;
}

if((kin == "SBS9") && (SBS_field == 70) && reject_bkgd && x[0]> -2.2  && x[0] < 0.75 ){
TF1::RejectPoint();
return 0;
}

return poly4_fit(x,param);
}

//Fit for protons using a Gaussian
double P_fit(double *x, double *param){
return gaussian_fit(x,param);
}

//Fit for neutrons using a Gaussian
double N_fit(double *x, double *param){
return gaussian_fit(x,param);
}


//Total fit to overlay. Combine all fits
//Scale a poly4 fit for background to the neutron and proton signal
//Does not totally work need a better signal sample. Probably from MC. Use old method till then

/*TH1D *hdx_justp;
TH1D *hdx_justn;
double Tot_fit(double *x, double *param){
double dx = x[0];
double sig_scale_p = param[0];
double sig_scale_n = param[1];
double signal = sig_scale_p*hdx_justp->Interpolate(dx) + sig_scale_n*hdx_justn->Interpolate(dx);
double tot_func = signal + poly4_fit(x,&param[2]);

return tot_func;

}*/
// Total fit to overlay.
// Currently 2 Gaussians and a 4th order poly
double Tot_fit(double *x, double *param){

double tot_func = BG_4fit(x, &param[0]) + P_fit(x, &param[5]) + N_fit(x, & param[8]);
return tot_func;
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
			else if(key == "W2_low"){
                	W2_low = val.Atof();
                        //cout << "W2 low " << W2_low << endl;
			}
			else if(key == "W2_high"){
			W2_high = val.Atof();
                        //cout << "W2 high " << W2_high << endl;
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
			else if(key == "dxmax"){
                        dxmax = val.Atof();
                        //cout << "max x difference between peaks" << dxmax << endl;
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
                        else if(key == "dxsig_n_fac"){
                        dxsig_n_fac = val.Atof();
                        //cout << "dx sigma factor for neutron" << dxsig_n_fac << endl;
                        }    
			else if(key == "dxsig_p_fac"){
                        dxsig_p_fac = val.Atof();
                        //cout << "dx sigma factor for proton" << dxsig_p_fac << endl;
                        } 
			else if(key == "dysig_fac"){
                        dysig_fac = val.Atof();
                        //cout << "dy sigma factor " << dysig_fac << endl;
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
  double sbs_dist = myData[0].getSBSDist(); //SBS dist

  double ntrack, nhits;
 //BBCal variables
  double BBtr_px[MAXNTRACKS], BBtr_py[MAXNTRACKS], BBtr_pz[MAXNTRACKS], BBtr_p[MAXNTRACKS];
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e,SHnclus,PSnclus, BB_E_over_p;

  //HCal variables
  double xhcal,yhcal,ehcal,nblk,nclus;
  double atime[maxHCalChan], row[maxHCalRows], col[maxHCalCols], cblkid[maxHCalChan], cblke [maxHCalChan];
  UInt_t TBits;
  double  hcal_atime,hcal_tdctime,bbcal_atime;
  double TDCT_id[maxTDCTrigChan], TDCT_tdc[maxTDCTrigChan], hodo_tmean[maxTDCTrigChan];
  int TDCTndata;
  double kineW2;

  // Cut on global parameters from setup config
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );
      
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

  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
 // C->SetBranchStatus( "sbs.hcal.clus.tdctime", 1 );
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
//need to change the bb time stuff
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus( "sbs.hcal.atimeblk", 1 );
  C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 );
  C->SetBranchStatus("e.kine.W2",1);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.gem.track.nhits",&nhits);
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
  //C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
  C->SetBranchAddress( "sbs.hcal.atimeblk", &hcal_atime );
  C->SetBranchAddress( "sbs.hcal.tdctimeblk", &hcal_tdctime );
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
  C->SetBranchAddress( "bb.etot_over_p", &BB_E_over_p );
  C->SetBranchAddress( "fEvtHdr.fTrigBits", &TBits ); 
  //need to change these bb times
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "bb.sh.atimeblk", &bbcal_atime );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  
  //Declare output file
  //Setup the name for the output file
  TString outfile = makeOutputFileName(Exp,kin,SBS_field,targ);
  TFile *fout = new TFile(outfile,"RECREATE");

  
   int nentries = C->GetEntries();

    
  // Initialize histograms
  
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

 // TH1D *h_atime = new TH1D( "atime", "HCal ADC Time, All Channels; ns", 160, 0, 160 );
 // TH2D *h_CvCh = new TH2D( "CvCh", "HCal Coeff Single Block Clusters; channel, GeV", maxHCalChan, 0, maxHCalChan, 200, 0, 1.0 );
 // TH1D *h_E_exp = new TH1D( "E_exp", "Expected Energy Dep in HCal; GeV", 100, 0, 0.2 );
  TH2D *hdxVE = new TH2D("dxVE","dx vs VE;x_{HCAL}-x_{expect} (m); GeV", 100, -4.0, 2.0, 100, 0, 0.5 );
  TH1D *hKElow = new TH1D( "KElow", "Lowest Elastic E Sampled in HCal (GeV)", 500, 0.0, 0.2 );
  TH1D *h_hcalatime = new TH1D("HCal_atime","HCal ADC Time; ns", 160,-50,200);
  TH1D *h_bbcalatime = new TH1D("BBCal_atime","BBCal ADC Time; ns", 160,-60,60); 
  TH1D *h_hcaltdctime = new TH1D("HCal_tdctime","HCal TDC Time; ns", 160,-220,80);

  TH1D *h_hcalatime_cut = new TH1D("HCal_atime_cut","HCal ADC Time,Cuts; ns", 160,-50,200);
  TH1D *h_bbcalatime_cut = new TH1D("BBCal_atime_cut","BBCal ADC Time,Cuts; ns", 160,-60,60);
  TH1D *h_hcaltdctime_cut = new TH1D("HCal_tdctime_cut","HCal TDC Time,Cuts; ns", 160,-220,80);

  TH1D *h_hcaltdctime_bbcalatime = new TH1D( "hHCalTDCBBCalA","HCal TDC time - BBCal time (ns)", 1300, -100, 75 );
  TH1D *h_hcalatime_bbcalatime = new TH1D( "hHCalABBCalA","HCal ADC time - BBCal time (ns)", 1300, -50, 100 );
 
  TH1D *h_hcaltdctime_bbcalatime_cut = new TH1D( "hHCalTDCBBCalA_cut","HCal TDC time - BBCal time (ns),Cuts", 1300, -100, 75 );
  TH1D *h_hcalatime_bbcalatime_cut = new TH1D( "hHCalABBCalA_cut","HCal ADC time - BBCal time (ns),Cuts", 1300, -50, 100 );

  

 //W2 and timing info
 // TH2D *h_EvCh = new TH2D( "EvCh", "HCal Cluster E Single Block Clusters; channel, GeV", maxHCalChan, 0, maxHCalChan, 50, 0, 0.5 );
  TH1D *h_W2 = new TH1D( "W2", "W2 (GeV) No Cuts; GeV", 250, -1.0, 6.0 );
  TH1D *h_W2recon = new TH1D( "W2recon", "W2 (GeV) No Cuts; GeV", 250, -1.0, 6.0 );
  TH1D *h_W2_cut = new TH1D( "W2_cut", "W2 (GeV) with cuts; GeV", 250, -1.0, 3.0 );
  TH1D *htimeDiff = new TH1D( "hDiff","HCal time - BBCal time (ns)", 1300, -500, 800 );
  TH1D *htimeDiff_cut = new TH1D( "hDiff_cut","HCal time - BBCal time (ns), All Cuts", 150, 450, 600 );
  TH1D *h_W2_pncut = new TH1D( "W2_pncut", "W2_pncut; GeV", 250, -1.0, 3.0 );
  

  //Various HCal related yield plots
  TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal X  vs Y, no cuts;HCal Y  (m); HCal X  (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
  TH2D *hxy_acceptancecut = new TH2D("hxy_accpetancecut","HCal X  vs Y, acceptance cut;HCal Y  (m); HCal X  (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
  TH2D *hxy_expect_nocut = new TH2D("hxy_expect_nocut","HCal X Expect vs Y Expect, no cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_cut = new TH2D("hxy_cut","HCal X  vs Y, cuts;HCal Y (m); HCal X (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
  TH2D *hxy_expect_cut = new TH2D("hxy_expect_cut","HCal X Expect vs Y Expect, cuts;HCal Y Expect (m); HCal X Expect (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
  
  TH2D *hdxdy_cut = new TH2D("dxdy_cut","HCal dxdy Cuts, All Channels;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 350, -1.25, 1.25, 350, dx_low, dx_high );
  TH2D *hdxdy_nocut = new TH2D("dxdy_nocut","HCal dxdy ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",350,-1.25,1.25,350,dx_low,dx_high);
  TH1D *hdx = new TH1D( "dx", "HCal dx (m); m", 250, dx_low, dx_high );
  TH1D *hdy = new TH1D( "dy", "HCal dy (m); m", 250, dy_low, dy_high );
  TH1D *hX = new TH1D( "X", "HCal X (m); m", 250, dx_low, dx_high );
  TH1D *hX_expect = new TH1D( "X_expect", "HCal X Expect (m); m", 250, dx_low, dx_high+1 );
  TH1D *hY = new TH1D( "Y", "HCal Y (m); m", 250, dy_low, dy_high );
  TH1D *hY_expect = new TH1D( "Y_expect", "HCal Y Expect (m); m", 250, dy_low, dy_high );
  TH1D *hdx_cut = new TH1D( "dx_cut", "HCal dx (m),  cuts; m", 250, dx_low, dx_high );
  TH1D *hdy_cut = new TH1D( "dy_cut", "HCal dy (m),  cuts; m", 250, dy_low, dy_high );
  TH1D *hX_cut = new TH1D( "X_cut", "HCal X (m),  cuts; m", 250, dx_low, dx_high );
  TH1D *hX_expect_cut = new TH1D( "X_expect_cut", "HCal X Expect (m), cuts; m", 250, dx_low, dx_high+1 );
  TH1D *hY_cut = new TH1D( "Y_cut", "HCal Y (m), cuts; m", 250, dy_low, dy_high );
  TH1D *hY_expect_cut = new TH1D( "Y_expect_cut", "HCal Y Expect (m),  cuts; m", 250, dy_low, dy_high );
  TH2D *hdxdy_pcut = new TH2D("hdxdy_pcut","HCal dxdy, proton fiducial;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",350,-1.25,1.25,350,dx_low,dx_high);
  TH2D *hdxdy_ncut = new TH2D("hdxdy_ncut","HCal dxdy, nuetron fiducial;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",350,-1.25,1.25,350,dx_low,dx_high);
  TH2D *hdxdy_pncut = new TH2D("hdxdy_pncut","HCal dxdy, just pn;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",350,-1.25,1.25,350,dx_low,dx_high);
  TH1D *hdx_pncut = new TH1D( "dx_pncut","HCal dx just pn; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
  TH1D *hdx_pcut = new TH1D( "dx_pcut","HCal dx just pcut; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
  TH1D *hdx_ncut = new TH1D( "dx_ncut","HCal dx just ncut; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
  TH1D *hdy_pncut = new TH1D( "dy_pncut","HCal dy just pn; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );
  TH2D *hxy_pncut = new TH2D("hxy_pncut","HCal X vs Y, just pn;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -2.5, 2.5);
  TH2D *hxy_pcut = new TH2D("hxy_pcut","HCal X vs Y, proton ;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -2.5, 2.5);
  TH2D *hxy_ncut = new TH2D("hxy_ncut","HCal X vs Y, neutron;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -2.5, 2.5);
  TH2D *hxy_expect_pcut = new TH2D("hxy_expect_pcut","HCal X Expect vs Y Expect, pcuts;HCal Y Expect (m); HCal X Expect (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
  TH2D *hxy_expect_ncut = new TH2D("hxy_expect_ncut","HCal X Expect vs Y Expect, ncuts;HCal Y Expect (m); HCal X Expect (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
  TH1D *hdx_anticut = new TH1D( "dx_anticut", "HCal dx (m),  anticut; m", 250, dx_low, dx_high );
  TH2D *hdxdy_anticut = new TH2D("dxdy_anticut","HCal dxdy anticut;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 350, -1.25, 1.25, 350, dx_low, dx_high );
  

  //Fiducial cut histograms
  TH2D *hdxdy_fidcut = new TH2D("hdxdy_fidcut","HCal dxdy, fiducial cut;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",350,-1.25,1.25,350,dx_low,dx_high);
  TH1D *hdx_fidcut = new TH1D( "dx_fidcut","HCal dx fiducial cut; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
  TH1D *hdy_fidcut = new TH1D( "dy_fidcut","HCal dy fiducial cut; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );
  TH2D *hxy_fidcut = new TH2D("hxy_fidcut","HCal X vs Y, fiducial cut;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -2.5, 2.5);
  TH1D *h_W2_fidcut = new TH1D( "W2_fidcut", "W2_fidcut; GeV", 250, -1.0, 3.0 );
  TH2D *hxy_expect_fidcut = new TH2D("hxy_expect_fidcut","HCal X Expect vs Y Expect, fiducial cut;HCal Y Expect (m); HCal X Expect (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
  TH1D *hdx_fidanticut = new TH1D( "dx_fidanticut", "HCal dx (m), Fid  anticut; m", 250, dx_low, dx_high );  
  TH2D *hdxdy_fidanticut = new TH2D("dxdy_fidanticut","HCal dxdy Fid anticut;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 350, -1.25, 1.25, 350, dx_low, dx_high );

  //From physics calculations
  TH1D *h_dpel = new TH1D("h_dpel","d_pel;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *h_dpel_cut = new TH1D("h_dpel_cut","d_pel, Cuts;p/p_{elastic}(#theta)-1;",250,-0.25,0.25);
  TH1D *hE_ep = new TH1D("E_ep","Scattered Electron Energy",500,0.0,Ebeam*1.25);
  hE_ep->GetXaxis()->SetTitle("GeV"); 
  TH1D *hE_eloss = new TH1D("E_eloss","Scattered Electron Energy Loss in Target",500,-0.001, Ebeam *0.003);
  hE_eloss->GetXaxis()->SetTitle("GeV");
  TH1D *hE_ecorr = new TH1D("E_ecorr","Corrected Scattered Electron Energy",500,Ebeam*0.995,Ebeam*1.005);
  hE_ecorr->GetXaxis()->SetTitle("GeV");
  TH2D *hE_ecorr_vs_vertex = new TH2D("hE_ecorr_vs_vertex",";z_{vertex} (m); E_{corr} (GeV)", 500, -0.15, 0.15, 500, Ebeam*0.995,Ebeam*1.005);
  TH2D *hE_eloss_vs_vertex = new TH2D("hE_eloss_vs_vertex",";z_{vertex} (m); E_{loss} (GeV)", 500, -0.15, 0.15, 500, -0.001, Ebeam *0.003);
  TH1D *hE_pp = new TH1D("E_pp","Scattered Proton Energy",500, 0.0, Ebeam*1.25);
  hE_pp->GetXaxis()->SetTitle("GeV");
  TH1D *hKE_p = new TH1D( "KE_p", "Scattered Proton Kinetic Energy", 500, 0.0, Ebeam*1.25 );
  hKE_p->GetXaxis()->SetTitle("GeV");

  TH1D *hQ2 = new TH1D("Q2","Q2",250,0.5,6.0);
  hQ2->GetXaxis()->SetTitle("GeV");
  TH1D *h_thetapq_nocut = new TH1D("h_thetapq_nocut","Theta_pq",150, -0.05 ,0.7);
  TH2D *h_W2_vs_thetapq_nocut = new TH2D("h_W2_theapq_nocut","W2 vs Theta pq, No Cuts;theta pq; W2 (GeV)",150, -0.05, 0.7, 250, 0.0, 6.0 );
  TH1D *h_thetapq_cut = new TH1D("h_thetapq_cut","Theta_pq",150, -0.05 ,0.5);
  TH2D *h_W2_vs_thetapq_cut = new TH2D("h_W2_theapq_cut","W2 vs Theta pq, Cuts;theta pq; W2 (GeV)",150, -0.05, 0.5, 250, 0.0, 1.75 );
 

  double pBeam = Ebeam/(1.0 +(Ebeam/M_p)*(1.0-cos(bbtheta)));

 //mean energy loss of the beam before scattering
 double Eloss_outgoing = (celldiameter/2)/sin(bbtheta)*rho_tgt*dEdx_tgt; //Approximately 1 MeV, one could correct more using the raster position
 if(useAlshield !=0){
 Eloss_outgoing += Alshieldthick * rho_Al*dEdx_Al;
 }

 

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
  //convert length to cm so units cancel properly
  double Eloss = ((vz[0]+(l_tgt/2))*100)*rho_tgt*dEdx_tgt + uwallthick_LH2*rho_Al*dEdx_Al; //aproximately 3 MeV
  hE_eloss->Fill(Eloss);
  hE_eloss_vs_vertex->Fill(vz[0],Eloss);

  double Ecorr = Ebeam-Eloss;
  hE_ecorr->Fill(Ecorr);
  hE_ecorr_vs_vertex->Fill(vz[0],Ecorr);

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

  //Define the expected position of hadron on HCal from BB track. This assumes a neutron hypothesis since the expect values for the HCal positions do not account for the deflection of the protons due to the magnetic field.
  double xhcal_expect = (hcal_intersect-hcal_origin).Dot( hcal_xaxis );
  double yhcal_expect = (hcal_intersect-hcal_origin).Dot( hcal_yaxis );

  //calculate quantities of interest
  double W2 = kineW2;//Get the invariant mass transfer W and four-momentum of the scattered nucleon
  //double W2recon = PgammaN.M2();//M2 is magnitude squared of the 4-vector
  //double Wrecon = PgammaN.M();
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
   double BdL =  maxSBSfield * Dgap * (SBS_field/100); //scale crudely by magnetic field
   double proton_deflection = tan( 0.3 * BdL / qvect_recon.Mag() ) * (hcaldist - (sbs_dist + Dgap/2.0) );
   TVector3 ProtonDirection = (hcalpos + proton_deflection * hcal_xaxis - vertex).Unit();
   double thetapq = acos( ProtonDirection.Dot( qvect_recon.Unit() ) );
   
 
  //define dx,dy
  double dx = xhcal - xhcal_expect;
  double dy = yhcal - yhcal_expect;
  
 
  //Cut on BBCal and HCal trigger coincidence. 
  double bbcal_time=0., hcal_time=0., rf_time=0., edtm_time = 0.,bbcalL0_time=0.,bbcalHiVeto_time=0.;
  double hcaltdc_bbcal = hcal_tdctime - bbcal_atime; 
  double hcala_bbcal = hcal_atime -bbcal_atime;
  /*  for(int ihit=0; ihit<TDCTndata; ihit++){
      if(TDCT_id[ihit]==5) bbcal_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==4) rf_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==3) edtm_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==2) bbcalL0_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==1) bbcalHiVeto_time=TDCT_tdc[ihit];
      if(TDCT_id[ihit]==0) hcal_time=TDCT_tdc[ihit];
      
     }*/

  //Fill histos for global cuts
  h_ntracks->Fill(ntrack);
  h_PS_E->Fill(BBps_e);
  h_vert_z->Fill(vz[0]);
  h_HCal_E->Fill(ehcal);
  h_HCal_nclus->Fill(nclus);
  h_TPS_SH->Fill(BBps_e+BBsh_e);
  h_nhits->Fill(nhits);
  h_bbtrp_nocut->Fill(BBtr_p[0]); 
  h_bbEoverp_nocut->Fill(BB_E_over_p);

 //cout << hcal_time << " " << bbcal_time << endl;
  //double diff = hcal_time - bbcal_time;
  //
  //Fill some histograms
  //cout << yhcal - yhcal_expect << endl;
  //htimeDiff->Fill( diff );
  hdxdy_nocut->Fill( dy, dx);
 
  h_W2->Fill(W2);
  h_W2recon->Fill(W2recon);
  hxy_nocut->Fill(yhcal,xhcal);
  hxy_expect_nocut->Fill(yhcal_expect,xhcal_expect);
  hX_expect->Fill( xhcal_expect);
  hY_expect->Fill( yhcal_expect);
  hdx->Fill(dx );
  hdy->Fill(dy );
  hX->Fill( xhcal);
  hY->Fill( yhcal );
  h_dpel->Fill(dpel);
  h_thetapq_nocut->Fill(thetapq); 
  h_W2_vs_thetapq_nocut->Fill(thetapq,W2);
 
  h_hcalatime->Fill(hcal_atime);
  h_bbcalatime->Fill(bbcal_atime);
  h_hcaltdctime->Fill(hcal_tdctime); 
  h_hcaltdctime_bbcalatime->Fill(hcaltdc_bbcal);
  h_hcalatime_bbcalatime->Fill(hcala_bbcal);

  //Define some HCal parameters
  //No events which project hadrons off the face of HCal should be considered. Call this acceptance cut. This inherently relies on the MC DB values. Could also consider not that.
  double HCal_left = posHCalYf_MC-HCalblk_l_h_MC;
  double HCal_right = posHCalYi_MC+HCalblk_l_h_MC;
  double HCal_top = posHCalXi_MC+HCalblk_l_v_MC;
  double HCal_bot = posHCalXf_MC-HCalblk_l_v_MC;
  
  //Primary Cuts
  bool offhcal =
  (yhcal < HCal_right) ||
  (yhcal > HCal_left) ||
  (xhcal < HCal_top) ||
  (xhcal > HCal_bot);

  //Base level acceptance cut for yield analysis
  if( offhcal ) continue;
  
  hxy_acceptancecut ->Fill(yhcal,xhcal);
   
  bool goodW2 = (W2 >= W2_low) && (W2 <= W2_high);  
  bool bad_dy = (abs(dy-dyO_n) > (dysig_fac*dysig_n)) || (abs(dy-dyO_p)>(dysig_fac*dysig_p));
  if(goodW2 && !failedglobal){  //Observed W2 peak
  
  h_bbEoverp_cut->Fill(BB_E_over_p);
  h_bbtrp_cut->Fill(BBtr_p[0]);
  hdxdy_cut->Fill( dy, dx );
  hY_cut->Fill( yhcal );
  hdy_cut->Fill( dy );
  hxy_cut->Fill(yhcal,xhcal);
  hxy_expect_cut->Fill(yhcal_expect,xhcal_expect);
  
  //cut on dy to make life easier
  if(!bad_dy){

  h_W2_cut->Fill(W2);
  h_vert_z_cut->Fill(vz[0]);
  h_PS_E_cut->Fill(BBps_e);
  h_HCal_E_cut->Fill(ehcal);
  h_HCal_nclus_cut->Fill(nclus);
  h_TPS_SH_cut->Fill(BBps_e+BBsh_e);
  h_nhits_cut->Fill(nhits);

  //Fill some histograms

  hKElow->Fill( KE_p*HCalSampFrac );
  hdxVE->Fill(dx,ehcal);
  //Fill delta plots and others
  hdx_cut->Fill( dx );
  //h_atime->Fill( atime[0] );
  hX_cut->Fill( xhcal);
  hX_expect_cut->Fill( xhcal_expect);
  hY_expect_cut->Fill( yhcal_expect);
  //htimeDiff_cut->Fill( diff );
  h_dpel_cut->Fill(dpel);
  h_thetapq_cut->Fill(thetapq);
  h_W2_vs_thetapq_cut->Fill(thetapq,W2);

 

 

  h_hcalatime_cut->Fill(hcal_atime);
  h_bbcalatime_cut->Fill(bbcal_atime);
  h_hcaltdctime_cut->Fill(hcal_tdctime);
  h_hcaltdctime_bbcalatime_cut->Fill(hcaltdc_bbcal);
  h_hcalatime_bbcalatime_cut->Fill(hcala_bbcal); 
   
  
  
 //Determine proton and neutron selection from the dx dy plots
 
 bool find_p = false, find_n = false, find_both = false;

 //equation for an ellipse around proton spot. 
 if(pow((dx-dxO_p)/(dxsig_p_fac*dxsig_p),2)+pow((dy-dyO_p)/(dysig_fac*dysig_p),2)<= pow(1,2)){
 find_p = true;
 }
 //equation for an ellipse around neutron spot.
 if(pow((dx-dxO_n)/(dxsig_n_fac*dxsig_n),2)+pow((dy-dyO_n)/(dysig_fac*dysig_n),2)<= pow(1,2)){
 find_n = true;
 }
 if(find_p && find_n){
 find_both = true;
 }
 //Fill 2D histograms if we find protons or neutrons. Fills only if true
 if(find_p){
 hdxdy_pcut->Fill(dy,dx);
 }
 if(find_n){
 hdxdy_ncut->Fill(dy,dx);
 }

 
 //We need to make sure we will see both neutrons and protons giving the kinematic information
 //fills only if true
 bool neutron_hyp = ((xhcal_expect - dxmax) >= HCal_top) && xhcal_expect <= HCal_bot ;
 bool proton_hyp = ((xhcal_expect + dxmax) <= HCal_bot) && xhcal_expect >= HCal_top ;


 /* if(find_both){
  if(neutron_hyp && proton_hyp){
        hdxdy_pncut->Fill(dy,dx);
        hdx_pncut->Fill(dx);
        hdy_pncut->Fill(dy);
        h_W2_pncut->Fill(W2);
        hxy_pncut->Fill(yhcal,xhcal);
       
  }*/
  if(find_n){
 	if(neutron_hyp){
	//HCal_bot =(posHCalXi_MC+HCalblk_l_v_MC) = 1.29434 which is the edge of HCal we need to worry about find neutrons at a kinematic
	hdxdy_pncut->Fill(dy,dx);
        hdx_pncut->Fill(dx);
	hdx_ncut->Fill(dx);
	hdy_pncut->Fill(dy);
	h_W2_pncut->Fill(W2);
	hxy_pncut->Fill(yhcal,xhcal);
	hxy_ncut->Fill(yhcal,xhcal);
        hxy_expect_ncut->Fill(yhcal_expect,xhcal_expect);     
	}
 
 }else if (find_p){
 	if(proton_hyp){
	//HCal_top =(posHCalXf_MC-HCalblk_l_v_MC) = -2.19435 which is the edge of HCal we need to worry about find protons at a kinematic
	hdxdy_pncut->Fill(dy,dx);
        hdx_pncut->Fill(dx);
        hdx_pcut->Fill(dx);
	hdy_pncut->Fill(dy);
        h_W2_pncut->Fill(W2);
        hxy_pncut->Fill(yhcal,xhcal);
        hxy_pcut->Fill(yhcal,xhcal);
        hxy_expect_pcut->Fill(yhcal_expect,xhcal_expect);	
	}
       }
 double HCal_fid_top = HCal_top - dxsig_p;
 double HCal_fid_bot = HCal_bot + dxsig_n;
 double HCal_fid_left = HCal_left + dysig_p; //currently dysig is same for both p and n
 double HCal_fid_right = HCal_right - dysig_p;
 //Actual fiducial cut. This assumes a neutron and proton hypothesis since the expect values for the HCal positions do not account for the deflection of the protons due to the magnetic field.
 bool neutron_hyp_fid = ((xhcal_expect - dxmax) >= HCal_fid_top) && xhcal_expect <= HCal_bot ;
 bool proton_hyp_fid = ((xhcal_expect + dxmax) <= HCal_fid_bot) && xhcal_expect >= HCal_top ;
 
 bool inFiducial = (neutron_hyp_fid || proton_hyp_fid);
      if(inFiducial){
        hdxdy_fidcut->Fill(dy,dx);
        hdx_fidcut->Fill(dx);
        hdy_fidcut->Fill(dy);
        h_W2_fidcut->Fill(W2);
        hxy_fidcut->Fill(yhcal,xhcal);
        hxy_expect_fidcut->Fill(yhcal_expect,xhcal_expect);
      }
     /* else{
      //else for fiducial cut
      //cout << "The Fiducial cut actually failed for once." << endl;
      hdx_fidanticut->Fill(dx);
      hdxdy_fidanticut->Fill(dy,dx);
       }*/
     }else{
     //else for dy cut    
      hdx_anticut->Fill(dx);
      hdxdy_anticut->Fill(dy,dx);
     }
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
 
 

//TH1D *hdx_fit = new TH1D( "dx_fit", "HCal dx fit (m), All cuts; m", 250, dx_low, dx_high );

//setup residual histo
TH1D *hdx_residual = (TH1D*)hdx_fidcut->Clone("hdx_residual");
hdx_residual->SetTitle("HCal dx residual (data - fit),All cuts");
hdx_residual->GetXaxis()->SetTitle("m");

//setup background subtracted histo
TH1D *hdx_BGSub = (TH1D*)hdx_fidcut->Clone("hdx_BGSub");
hdx_BGSub->SetTitle("HCal dx BG subtracted");
hdx_BGSub->GetXaxis()->SetTitle("m");


 //Declare Canvas
 TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
 c1->Divide(2,1);
 //place holder till I do fiducials
 c1->cd(1);
 hdxdy_fidcut->Draw("colz"); 

 TEllipse el_pro;
 el_pro.SetFillStyle(4005);
 el_pro.SetLineColor(2);
 el_pro.SetLineWidth(3);
 el_pro.DrawEllipse(dyO_p,dxO_p,sqrt(dysig_fac)*dysig_p,sqrt(dxsig_p_fac)*dxsig_p,0,360,0);

 
 TEllipse el_neu;
 el_neu.SetFillStyle(4005);
 el_neu.SetLineColor(4);
 el_neu.SetLineWidth(3);
 el_neu.DrawEllipse(dyO_n,dxO_n,sqrt(dysig_fac)*dysig_n,sqrt(dxsig_n_fac)*dxsig_n,0,360,0);


 //Sets to no fitting?
 gStyle->SetOptFit(0);
 //Set to no stat box?
 gStyle->SetOptStat(0);

 //Working on fitting the dx plot
 c1->cd(2);
 //make a clone of the dx plot
 // hdx_justp = (TH1D*)hdx_pcut->Clone("hdx_justp");
 // hdx_justn = (TH1D*)hdx_ncut->Clone("hdx_justn");
  TH1D *hdx_fidcut_clone = (TH1D*)hdx_fidcut->Clone("hdx_fidcut_clone");
 //Initialize fit parameters
 
 vector<double> myParam (13);
 //cout << "kin: " << kin << " SBS_field:  " << SBS_field << " targ:  " << targ << endl;
 myParam= fit_Params(kin,SBS_field,targ);
 double fit_low = myParam[11];
 double fit_high = myParam[12];
 


 //Make the total fit
 c1->cd(2);
 TF1 *totalFit = new TF1("totalfit",Tot_fit,fit_low,fit_high,11);
 totalFit->SetLineWidth(4);
 totalFit->SetLineColor(kMagenta);
 totalFit->SetNpx(250);
 //Set or Fix parameter. Right now there are 11
 
 totalFit->SetParameters(&myParam[0]);

 //Q = minimum printing, R =fit using fitting range specified in the function range, B= Used when fixing or setting limits on one or more parameters in predefined function, + = add to list of fitted functions, V = verbose
 hdx_fidcut_clone->Fit("totalfit","RB+");
 //totalFit->GetParameters(myFParam.data());
 double *myFParam = totalFit->GetParameters();
 TH1D *hdx_fit = (TH1D*) (totalFit->GetHistogram())->Clone("hdx_fit"); 
 hdx_residual->Add(hdx_fit,-1);
 /*for(int i=0; i<myParam.size();i++){
 cout << myParam[i] << " " << myFParam[i] << endl;
 
 }*/
 

 TF1 *bkgd = new TF1("bkgd",BG_4fit,fit_low,fit_high,5);
 TF1 *proton = new TF1("proton",P_fit,fit_low,fit_high,3);
 TF1 *neutron = new TF1("neutron",N_fit,fit_low,fit_high,3);
 bkgd->SetLineColor(kBlack);
 proton->SetLineColor(kRed);
 neutron->SetLineColor(kBlue);

 bkgd->SetNpx(500);
 proton->SetNpx(500);
 neutron->SetNpx(500);
 
/* bkgd->SetParameters(&myFParam[0]);
 proton->SetParameters(&myFParam[5]);
 neutron->SetParameters(&myFParam[8]);*/
 
 bkgd->SetParameters(&myParam[0]);
 proton->SetParameters(&myParam[5]);
 neutron->SetParameters(&myParam[8]);



 bkgd->Draw("same");
 proton->Draw("same");
 neutron->Draw("same");
 
 //Generate yields
 double p_yield = proton->Integral(dx_low,dx_high)/hdx_fidcut_clone->GetBinWidth(1);
 double n_yield = neutron->Integral(dx_low,dx_high)/hdx_fidcut_clone->GetBinWidth(1);

 //blind yields
 blind_factor *p_blind = new blind_factor("GetOffMyLawnYou");
 blind_factor *n_blind = new blind_factor("ILaughedTheToaster");
 double p_tot = p_blind->blind(p_yield);
 double n_tot = n_blind->blind(n_yield);
 double ratio = n_tot/p_tot;
// cout << p_blind->getBlindFac() << " " << n_blind->getBlindFac() << endl;

 //make legend
 auto legend = new TLegend(0.5,0.7,0.9,0.9);
 legend->SetTextSize(0.03);
 legend->AddEntry(bkgd,"Background Fit","l");
 legend->AddEntry(proton,Form("Proton Fit, Yield: %d",(int)p_tot),"l");
 legend->AddEntry(neutron,Form("Neutron Fit, Yield: %d",(int)n_tot),"l");
 legend->AddEntry(totalFit,"Total Fit","l");
 legend->AddEntry((TObject*)0,Form("Ratio, R: %f",ratio),"");
 legend->Draw();


 //Subtract the background from the data
 hdx_BGSub->Add(bkgd,-1);
 for( int i=1; i<=hdx_BGSub->GetNbinsX(); i++ ){
 double binc = hdx_BGSub->GetBinContent(i);
  if( binc<0.0 ){
  hdx_BGSub->SetBinContent(i,0.0);
  }
 } 


 c1->Write();

 TString reportfile = makeYieldReportFileName(Exp,kin,SBS_field,targ);
 //Declare outfile
 ofstream report;
 report.open( reportfile );
 
 report << "Yield report for SBS" << kin << " LD2 data at " << SBS_field << " percent field" << endl << endl;

 cout << "Chi^2 for total fit: " << totalFit->GetChisquare() << endl;
 report << "Chi^2 for total fit: " << totalFit->GetChisquare() << endl;

 cout << "Blinded neutron yield: " << (int)n_tot << endl;
 report << "Blinded neutron yield: " << (int)n_tot << endl;

 cout << "Blinded proton yield: " << (int)p_tot << endl;
 report << "Blinded proton yield: " << (int)p_tot << endl;

 cout << "Ratio of neutron over proton yields: " << ratio << endl;
 report << "Ratio of neutron over proton yields: " << ratio << endl;

 report.close();

  //output the info to some files
  TString plotname = outfile;
  plotname.ReplaceAll(".root",".pdf");
  c1->Print(plotname.Data(),"pdf");
  plotname.ReplaceAll(".pdf",".png");
  c1->Print(plotname.Data(),"png");
  fout->Write();


  watch->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;



} 
