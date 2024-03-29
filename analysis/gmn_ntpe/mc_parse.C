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
bool sync_jobs = mainConfig.get_syncJobs();
double W2_mean = mainConfig.getW2Mean();
double W2_sigma = mainConfig.getW2Sigma();
double W2_sigfac = mainConfig.getW2SigFac();
double dxO_n = mainConfig.get_dxOn();
double dyO_n = mainConfig.get_dyOn();
double dxsig_n = mainConfig.get_dxsign();
double dysig_n = mainConfig.get_dysign();
double dxO_p = mainConfig.get_dxOp();
double dyO_p = mainConfig.get_dyOp();
double dxsig_p = mainConfig.get_dxsigp();
double dysig_p = mainConfig.get_dysigp();
double dx_pn = mainConfig.get_dxpn();
double dx_low = mainConfig.get_dxLow();
double dx_high = mainConfig.get_dxHigh();
double dy_low = mainConfig.get_dyLow();
double dy_high = mainConfig.get_dyHigh();
double dxsig_n_fac = mainConfig.get_dxSignFac();
double dxsig_p_fac = mainConfig.get_dxSigpFac();
double dysig_n_fac = mainConfig.get_dySignFac();
double dysig_p_fac = mainConfig.get_dySigpFac();
double W2fitmax = mainConfig.getW2FitMax();
double binfac = mainConfig.getBinFac();
double hbinfac = mainConfig.getHBinFac();
int e_method = mainConfig.get_emethod();

double W2_low = W2_mean - W2_sigfac*W2_sigma;
double W2_high = W2_mean + W2_sigfac*W2_sigma;

double hcalfit_low = exp_constants::hcalposXi_mc; //lower fit/bin limit for hcal dx plots. 
double hcalfit_high = exp_constants::hcalposXf_mc; //higher fit/bin limit for hcal dx plots.
double hcal_fitrange = exp_constants::hcal_vrange; //Full range of hcal dx plots


//store all important kinematic info in local variables
kinematic_obj myKin(kinematic_file,kin);
double EBeam = myKin.getBeamEnergy();
double hcaldist = myKin.getHCalDist();
double sbsdist = myKin.getSBSDist();
double bbtheta = myKin.getBBAngle_Rad();
double hcaltheta = myKin.getHCalAngle_Rad();

//setup hcal physical bounds that match database
vector<double> hcalpos = cuts::hcal_Position_MC();

//setup hcal active area with bounds that match database
vector<double> hcalaa = cuts::hcal_ActiveArea_MC(1,1);

//setup fiducial region based on dx and dy spot information
vector<double> hcalfid = cuts::hcalfid(dxsig_p,dxsig_n,dysig_p,hcalaa,dxsig_p_fac,dysig_p_fac);

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

	//check to make number the number of hist and root files matches for protons and neutrons
	if(rootFileNames_p.size() != histFileNames_p.size() || rootFileNames_n.size() != histFileNames_n.size()){
	cerr << "Error: File Matching failure, vector size mismatch!" << endl;
	}

	//Continuing to make sure there exist both protons and neutrons. Make sure the job has a both a proton and neutron file. If not get rid of it
	if(sync_jobs){
	//Function that searchs, finds, and then removes any files that are not in both proton and neutron vectors
	utility::syncJobNumbers(rootFileNames_p,rootFileNames_n);
	}	
   
//double check that we have the same number of files again
int pFiles = rootFileNames_p.size();
int nFiles = rootFileNames_n.size();
cout << endl << "Number of proton files " << pFiles << " , Number of neutron files " << nFiles << endl;
	if(pFiles != nFiles){
	cout << endl << "Warning: Number proton and neutron root files loaded does not match. Investigate!" << endl << endl;
		if(sync_jobs){
		return;
		}
	}	
//Histograms///////	

//global cuts
  TH1D *h_ntracks = new TH1D("ntracks","Number of Tracks;", 150, 0, 5);
  TH1D *h_ntracks_globcut = new TH1D("ntracks_globcut","Number of Tracks,global cut;", 150, 0, 5);
  TH1D *h_ntracks_cut = new TH1D("ntracks_cut","Number of Tracks, cuts;", 150, 0, 5);
  TH1D *h_PS_E = new TH1D("h_ps_e"," PS Cluster Energy (GeV);",250,0.0,2.2);
  TH1D *h_PS_E_globcut = new TH1D("h_ps_e_globcut"," PS Cluster Energy (GeV),global cut;",250,0.0,2.2);
  TH1D *h_PS_E_cut = new TH1D("h_ps_e_cut"," PS Cluster Energy (GeV), Cuts;",250,0.0,2.2);
  TH1D *h_vert_z = new TH1D( "vert_z", "Vertex Position z-direction (m); m", 200, -0.2, 0.2 );
  TH1D *h_vert_z_globcut = new TH1D( "vert_z_globcut", "Vertex Position z-direction (m), global cut; m", 200, -0.2, 0.2 );
  TH1D *h_vert_z_cut = new TH1D( "vert_z_cut", "Vertex Position z-direction (m), Cuts; m", 200, -0.2, 0.2 );
  TH1D *h_HCal_E = new TH1D( "HCal_E", "HCal Cluster Energy (GeV); GeV", 250, 0, 0.4 );
  TH1D *h_HCal_E_globcut = new TH1D( "HCal_E_globcut", "HCal Cluster Energy (GeV), global cut; GeV", 250, 0, 0.4 );
  TH1D *h_HCal_E_cut = new TH1D( "HCal_E_cut", "HCal Cluster Energy (GeV), Cuts; GeV", 250, 0, 0.4 );
  TH1D *h_HCal_nclus = new TH1D("HCal_nclus","HCal number of clusters meeting threshold;", 250,0,10);
  TH1D *h_HCal_nclus_globcut = new TH1D("HCal_nclus_globcut","HCal number of clusters meeting threshold, global cut;", 250,0,10);
  TH1D *h_HCal_nclus_cut = new TH1D("HCal_nclus_cut","HCal number of clusters meeting threshold, Cuts;", 250,0,10);
  TH1D *h_SH_nclus = new TH1D("SH_nclus","SH number of clusters meeting threshold;", 250,0,10);
  TH1D *h_SH_nclus_globcut = new TH1D("SH_nclus_globcut","SH number of clusters meeting threshold, global cut;", 250,0,10);
  TH1D *h_SH_nclus_cut = new TH1D("SH_nclus_cut","SH number of clusters meeting threshold, Cuts;", 250,0,10);
  TH1D *h_TPS_SH_globcut = new TH1D("h_tps_sh_globcut","Total PS and SH cluster energy (GeV),global cut;",250,0.5,4.0);
  TH1D *h_TPS_SH_cut = new TH1D("h_tps_sh_cut","Total PS and SH cluster energy (GeV), All cuts;",250,0.5,4.0);
  TH1D *h_nhits = new TH1D("nhits","Number of hits on track;",150, 0, 6);
  TH1D *h_nhits_globcut = new TH1D("nhits_globcut","Number of hits on track, global cut;",150, 0, 6);
  TH1D *h_nhits_cut = new TH1D("nhits_cut","Number of hits on track, all cuts;",150, 0, 6);
  TH1D *h_bbtrp_nocut = new TH1D("bbtrp_nocut","BigBite Track Momentum (GeV), no cut;",300, 0.0, 4.0);
  TH1D *h_bbtrp_globcut = new TH1D("bbtrp_globcut","BigBite Track Momentum (GeV), global cut;",300, 0.0, 4.0);
  TH1D *h_bbtrp_cut = new TH1D("bbtrp_cut","BigBite Track Momentum (GeV), cuts;",300, 0.0, 4.0);
  TH1D *h_bbEoverp_globcut = new TH1D("bbEoverp_globcut","BigBite E over p, global cut;",100, 0.0, 2.0);
  TH1D *h_bbEoverp_cut = new TH1D("bbEoverp_cut","BigBite E over p, cuts;",100, 0.0, 2.0);

//basic H-arm
  TH2D *hxy_globcut = new TH2D("hxy_globcut","HCal X  vs Y, global cut;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_glob_W2_cut = new TH2D("hxy_glob_W2_cut","HCal X  vs Y, global & W2 cuts;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_acceptancecut = new TH2D("hxy_accpetancecut","HCal X  vs Y,global & W2 & acceptance cuts, ;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_cut = new TH2D("hxy_cut","HCal X vs Y, all cuts, best cluster;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -3.0, 3.0);

//E-arm
  TH1D *h_W2_globcut = new TH1D( "W2_globcut", "W2 (GeV) global cut; GeV", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_glob_W2_cut = new TH1D( "W2_glob_W2_cut", "W2 (GeV) global & W2 cuts; GeV", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_cut = new TH1D( "W2_cut", "W2 (GeV) all cuts; GeV", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_Q2_globcut = new TH1D( "Q2_globcut", "Q2 (GeV) global cut; GeV", 300, 0.0, 6.0 );
  TH1D *h_Q2_cut = new TH1D( "Q2_cut", "Q2 (GeV) all cuts; GeV", 300, 0.0, 6.0 );
  TH2D *hxy_expect_globcut = new TH2D("hxy_expect_globcut","HCal X Expect vs Y Expect, global cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_glob_W2_cut = new TH2D("hxy_expect_glob_W2_cut","HCal X Expect vs Y Expect, global & W2 cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_n = new TH2D("hxy_expect_n","HCal X Expect vs Y Expect, elastic cuts neutron;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_p = new TH2D("hxy_expect_p","HCal X Expect vs Y Expect, elastic cuts proton;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_fidcutp = new TH2D("hxy_expect_fidcutp","HCal X Expect vs Y Expect, proton passed fiducial;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
 TH2D *hxy_expect_fidcutn = new TH2D("hxy_expect_fidcutn","HCal X Expect vs Y Expect, neutron passed fiducial;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_failedfid = new TH2D("hxy_expect_failedfid","HCal X Expect vs Y Expect, failed fiducial cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_cut = new TH2D("hxy_expect_cut","HCal X Expect vs Y Expect, all cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

//Both arms
  TH2D *hdxdy_globcut = new TH2D("dxdy_globcut","HCal dxdy, global cut ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_glob_W2_cut = new TH2D("dxdy_glob_W2_cut","HCal dxdy, global & W2 cuts ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","HCal dxdy, all cuts;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_nofid = new TH2D("dxdy_nofid","HCal dxdy, all cuts but fiducial ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_globcut_p = new TH2D("dxdy_globcut_p","HCal dxdy proton, global cut ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_cut_p = new TH2D("hdxdy_cut_p","HCal dxdy proton, all cuts;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);  
  TH2D *hdxdy_globcut_n = new TH2D("dxdy_globcut_n","HCal dxdy neutron, global cut ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_cut_n = new TH2D("hdxdy_cut_n","HCal dxdy neutron, all cuts;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxvE = new TH2D("dxvE","dx vs HCal E, all cuts;E_{HCAL} (GeV);x_{HCAL}-x_{expect} (m)", 400, 0.0, 4.0, 600, dx_low, dx_high );
  TH2D *hdxvW2 = new TH2D("dxvW2", "dx vs W2, all cuts; W2 (GeV); x_{HCAL}-x_{expect} (m)",binfac*W2fitmax, 0.0, W2fitmax, 600, dx_low,dx_high);
  TH1D *hdx_nocut = new TH1D( "dx_nocut", "HCal dx (m); m",hbinfac*hcal_fitrange ,hcalfit_low ,hcalfit_high);
  TH1D *hdx_cut = new TH1D( "dx_cut","HCal dx, all cuts; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_nofid = new TH1D( "dx_cut_nofid","HCal dx, all cuts but fiducial; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_failfid = new TH1D( "dx_cut_failfid","HCal dx, all cuts fail fiducial; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_globcut = new TH1D( "dx_globcut","HCal dx, global cut; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_glob_W2_cut = new TH1D( "dx_glob_W2_cut","HCal dx, global & W2 cuts; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );

  TH1D *hdx_globcut_p = new TH1D( "dx_globcut_p","HCal dx proton, global cut; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_p = new TH1D( "dx_cut_p","HCal dx proton, all cuts; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  
  TH1D *hdx_globcut_n = new TH1D( "dx_globcut_n","HCal dx neutron, global cut; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_n = new TH1D( "dx_cut_n","HCal dx neutron, all cuts; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );

  TH1D *hdy_nocut = new TH1D( "dy_nocut", "HCal dy (m), no cuts; m", 250, dy_low, dy_high );
  TH1D *hdy_globcut = new TH1D( "dy_globcut", "HCal dy (m), global cut; m", 250, dy_low, dy_high );
  TH1D *hdy_glob_W2_cut = new TH1D( "dy_glob_W2_cut", "HCal dy (m), global & W2 cuts; m", 250, dy_low, dy_high );
  TH1D *hdy_cut = new TH1D( "dy_cut","HCal dy, all cuts; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );




// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
