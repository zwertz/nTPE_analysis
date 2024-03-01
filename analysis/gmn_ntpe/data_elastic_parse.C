//Author Ezekiel Wertz
//02/01/2024
//Purpose: Parsing Script for LD2 data to produce output histograms for later analysis


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
#include "TString.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "../../src/utility.C"
#include "../../include/physics_constants.h"
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/parse_config.C"
#include "../../src/data_object.C"
#include "../../src/cuts.C"
#include "../../src/physics.C"

//Main
void data_elastic_parse(const char *setup_file_name){

  //Define a clock to check macro processing time
  TStopwatch *watch = new TStopwatch(); 
  watch->Start( kTRUE );

  //parse object to get in the information that The One Config file has and is manipulated
  parse_config mainConfig(setup_file_name); 
  //mainConfig.printDataYields();
  
  //store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
  vector<int> runNums = mainConfig.getRunNums();
  TCut globalcut = mainConfig.getGlobalCut();
  TString exp = mainConfig.getExp();
  TString kin = mainConfig.getKin();
  TString kinematic_file = mainConfig.getKinFileName();
  TString data_map = mainConfig.getDataFileName();
  TString pass = mainConfig.getPass();
  int sbs_field = mainConfig.getSBSField();
  int maxtracks = mainConfig.getMAXNTRACKS();
  TString target = mainConfig.getTarg();
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
  int useAlshield = mainConfig.getAlshield();
  double dxsig_n_fac = mainConfig.get_dxSignFac();
  double dxsig_p_fac = mainConfig.get_dxSigpFac();
  double dysig_n_fac = mainConfig.get_dySignFac();
  double dysig_p_fac = mainConfig.get_dySigpFac();
  double W2fitmax = mainConfig.getW2FitMax();
  double binfac = mainConfig.getBinFac();
  double hbinfac = mainConfig.getHBinFac();
  double hcal_offset = exp_constants::getHCalOffset(kin,pass);
  

  int W2_low = W2_mean - W2_sigfac*W2_sigma;
  int W2_high = W2_mean + W2_sigfac*W2_sigma;

  //store all important kinematic info in local variables
  kinematic_obj myKin(kinematic_file,kin);
  double EBeam = myKin.getBeamEnergy();
  double hcaldist = myKin.getHCalDist();
  double sbsdist = myKin.getSBSDist();
  double bbtheta = myKin.getBBAngle_Rad();
  double hcaltheta = myKin.getHCalAngle_Rad();
  
  //setup hcal active area with bounds that match database depending on pass
  vector<double> hcalaa = cuts::hcal_ActiveArea_data(1,1,pass);

  int num_runs = runNums.size();
  double hcalfit_low = exp_constants::hcalposXi_mc; //lower fit/bin limit for hcal dx plots. 
  double hcalfit_high = exp_constants::hcalposXf_mc; //higher fit/bin limit for hcal dx plots.
  double hcal_fitrange = exp_constants::hcal_vrange; //Full range of hcal dx plots

  //store all the data information we care about in a vector of data_objects for further use
  vector<data_object> myData;
  for(int i=0; i < num_runs; i++ ){
  //Run Num, date map name, kinematic map name, Kinematic, SBS Field, Target, Pass
  data_object myObj(runNums[i],data_map,kinematic_file,kin,sbs_field,target,pass);
  myData.push_back(myObj);
  myObj.printRunInfo();
  }

  //setup output file
  TString outfile = utility::makeOutputFileNameParse(exp,pass,kin,sbs_field,target);
  TFile *fout = new TFile(outfile,"RECREATE");


  //Histograms///////

  //global cuts
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

  //basic H-arm
  TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal X  vs Y, no cuts;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_acceptancecut = new TH2D("hxy_accpetancecut","HCal X  vs Y, acceptance cut, ;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_cut = new TH2D("hxy_cut","HCal X vs Y, all cuts, best cluster;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -3.0, 3.0);

  //E-arm
  TH1D *h_W2 = new TH1D( "W2", "W2 (GeV) No Cuts; GeV", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_cut = new TH1D( "W2_cut", "W2 (GeV) with cuts; GeV", binfac*W2fitmax, 0.0, W2fitmax );  
  TH2D *hxy_expect_nocut = new TH2D("hxy_expect_nocut","HCal X Expect vs Y Expect, no cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_n = new TH2D("hxy_expect_n","HCal X Expect vs Y Expect, elastic cuts neutron;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_p = new TH2D("hxy_expect_p","HCal X Expect vs Y Expect, elastic cuts proton;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_fidcutp = new TH2D("hxy_expect_fidcutp","HCal X Expect vs Y Expect, proton passed fiducial;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
 TH2D *hxy_expect_fidcutn = new TH2D("hxy_expect_fidcutn","HCal X Expect vs Y Expect, neutron passed fiducial;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_failedfid = new TH2D("hxy_expect_failedfid","HCal X Expect vs Y Expect, failed fiducial cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_cut = new TH2D("hxy_expect_cut","HCal X Expect vs Y Expect, all cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  //Both arms
  TH2D *hdxdy_nocut = new TH2D("dxdy_nocut","HCal dxdy ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","HCal dxdy, all cuts;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_nofid = new TH2D("dxdy_nofid","HCal dxdy, all cuts but fiducial ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxvE = new TH2D("dxvE","dx vs HCal E, all cuts;E_{HCAL} (GeV);x_{HCAL}-x_{expect} (m)", 400, 0.0, 4.0, 600, dx_low, dx_high );
  TH2D *hdxvW2 = new TH2D("dxvW2", "dx vs W2, all cuts; W2 (GeV); x_{HCAL}-x_{expect} (m)",binfac*W2fitmax, 0.0, W2fitmax, 600, dx_low,dx_high);
  TH1D *hdx_nocut = new TH1D( "dx_nocut", "HCal dx (m); m",hbinfac*hcal_fitrange ,hcalfit_low ,hcalfit_high);
  TH1D *hdx_cut = new TH1D( "dx_cut","HCal dx, all cuts; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_nofid = new TH1D( "dx_cut_nofid","HCal dx, all cuts but fiducial; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_failfid = new TH1D( "dx_cut_failfid","HCal dx, all cuts fail fiducial; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_global = new TH1D( "dx_cut_global","HCal dx, global cuts only; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_W2 = new TH1D( "dx_cut_W2","HCal dx, W2 cut only; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_earm = new TH1D( "dx_cut_earm","HCal dx, earm cuts only; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );

  TH1D *hdy = new TH1D( "dy", "HCal dy (m); m", 250, dy_low, dy_high );
  TH1D *hdy_cut = new TH1D( "dy_cut","HCal dy, all cuts; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );

  TH1D *hcoin = new TH1D( "hcoin", "HCal ADCt - BBCal ADCt, no cut; ns", 200, 0, 100 );
  TH1D *hcoin_cut = new TH1D( "hcoin_cut", "HCal ADCt - BBCal ADCt, cuts; ns", 200, 0, 100 );
  TH1D *hcoin_pclus = new TH1D( "hcoin_pclus", "HCal ADCt - BBCal ADCt, pclus, no cut; ns", 200, 0, 100 );
  TH1D *hcoin_pclus_cut = new TH1D( "hcoin_pclus_cut", "HCal ADCt - BBCal ADCt,pclus, cuts; ns", 200, 0, 100 );

  //general
  TH1D *hMott_cs = new TH1D( "hMott_cs", "Mott Cross Section, no cut; (GeV/c)^{-2}", 200, 0, 0.0001 );
  
  //allocate memory at each run
  TChain *C = nullptr;

  //create output tree
  TTree *Parse = new TTree("Parse","Analysis Tree");

  //uncut output tree variables
  double dx_out;
  double dy_out;
  double xexp_out;
  double yexp_out;
  double xhcal_out;
  double yhcal_out;
  double W2_out;
  double Q2_out;
  double mott_out;
  double ehcal_out;
  double BBtot_e_out;
  double BBsh_e_out;
  double BBps_e_out;
  double hcal_atime_out;
  double BBsh_atime_out;
  double BBps_atime_out;
  double BBgem_nhits_out;
  double BBtr_x_out;
  double BBtr_y_out;
  double BBtr_p_out;
  double coin_mean;
  double coin_sigma;
  double dyO_p_out;
  double dyO_n_out; 
  double dysig_p_out;
  double dysig_n_out;
  double dysig_n_fac_out;
  double dysig_p_fac_out;
  double dxO_p_out;
  double dxO_n_out;
  double dxsig_p_out;
  double dxsig_n_out;
  double dxsig_n_fac_out;
  double dxsig_p_fac_out;
  double W2mean_out;
  double W2sig_out;
  double W2sig_fac_out;

  int num_hcal_clusid_out;
  int hcal_clus_blk_out;
  int BBtr_n_out;
  int passFid_out;
  int run_out;
  int mag_out;

  //setup new output tree branches
  Parse->Branch("dx", &dx_out, "dx/D");
  Parse->Branch("dy", &dy_out, "dy/D");
  Parse->Branch("xexp", &xexp_out, "xexp/D");
  Parse->Branch("yexp", &yexp_out, "yexp/D");
  Parse->Branch("xhcal", &xhcal_out, "xhcal/D");
  Parse->Branch("yhcal", &yhcal_out, "yhcal/D");
  Parse->Branch("W2", &W2_out, "W2/D");
  Parse->Branch("Q2", &Q2_out, "Q2/D");
  Parse->Branch("mott", &mott_out, "mott/D");
  Parse->Branch("ehcal", &ehcal_out, "ehcal/D");
  Parse->Branch("BBtot_e", &BBtot_e_out, "BBtot_e/D");
  Parse->Branch("BBsh_e", &BBsh_e_out, "BBsh_e/D");
  Parse->Branch("BBps_e", &BBps_e_out, "BBps_e/D");
  Parse->Branch("hcal_atime", &hcal_atime_out, "hcal_atime/D");
  Parse->Branch("BBsh_atime", &BBsh_atime_out, "BBsh_atime/D");
  Parse->Branch("BBps_atime", &BBps_atime_out, "BBps_atime/D");
  Parse->Branch("BBgem_nhits", &BBgem_nhits_out, "BBgem_nhits/D");
  Parse->Branch("BBtr_x", &BBtr_x_out, "BBtr_x/D");
  Parse->Branch("BBtr_y", &BBtr_y_out, "BBtr_y/D");
  Parse->Branch("BBtr_p", &BBtr_p_out, "BBtr_p/D");
  Parse->Branch("coin_mean", &coin_mean, "coin_mean/D");
  Parse->Branch("coin_sigma", &coin_sigma, "coin_sigma/D");
  Parse->Branch("dyO_p", &dyO_p_out, "dyO_p/D");
  Parse->Branch("dyO_n", &dyO_n_out, "dyO_n/D");
  Parse->Branch("dysig_p", &dysig_p_out, "dysig_p/D");
  Parse->Branch("dysig_n", &dysig_n_out, "dysig_n/D");
  Parse->Branch("dysig_n_fac", &dysig_n_fac_out, "dysig_n_fac/D");
  Parse->Branch("dysig_p_fac", &dysig_p_fac_out, "dysig_p_fac/D");
  Parse->Branch("dxO_p", &dxO_p_out, "dxO_p/D");
  Parse->Branch("dxO_n", &dxO_n_out, "dxO_n/D");
  Parse->Branch("dxsig_p", &dxsig_p_out, "dxsig_p/D");
  Parse->Branch("dxsig_n", &dxsig_n_out, "dxsig_n/D");
  Parse->Branch("dxsig_n_fac", &dxsig_n_fac_out, "dxsig_n_fac/D");
  Parse->Branch("dxsig_p_fac", &dxsig_p_fac_out, "dxsig_p_fac/D");
  Parse->Branch("W2mean", &W2mean_out, "W2mean/D");
  Parse->Branch("W2sig", &W2sig_out, "W2sig/D");
  Parse->Branch("W2sig_fac", &W2sig_fac_out, "W2sig_fac/D");

  Parse->Branch("num_hcal_clusid", &num_hcal_clusid_out, "num_hcal_clusid/I");
  Parse->Branch("hcal_clus_blk", &hcal_clus_blk_out, "hcal_clus_blk/I");
  Parse->Branch("BBtr_n", &BBtr_n_out, "BBtr_n/I");
  Parse->Branch("passFid", &passFid_out, "passFid/I");
  Parse->Branch("run", &run_out, "run/I");
  Parse->Branch("mag", &mag_out, "mag/I");


  //loop over the run numbers
  for(int j = 0; j<num_runs; j++){

  //get run specific information
  data_object datData = myData[j];
  int run = datData.getRun();
  int field = datData.getSBSField();
  double Ebeam = datData.getBeamEnergy();
  double hcaldist = datData.getHCalDist();
  double hcaltheta = datData.getHCalAngle_Rad();
  double bbtheta = datData.getBBAngle_Rad();
  TString input_file_name = datData.getInputFile();
  
  //add the file to the TChain
  C = new TChain("T");
  
  C->Add(input_file_name);

  // setting up ROOT tree branch addresses
  C->SetBranchStatus("*",0);

  //HCal general branches
  double x_hcal,y_hcal,e_hcal,nclus_hcal,idx_hcal;

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.nclus",1); 
  C->SetBranchStatus("sbs.hcal.index",1);

  C->SetBranchAddress("sbs.hcal.x", &x_hcal);
  C->SetBranchAddress("sbs.hcal.y", &y_hcal);
  C->SetBranchAddress("sbs.hcal.e", &e_hcal);
  C->SetBranchAddress("sbs.hcal.nclus", &nclus_hcal);
  C->SetBranchAddress("sbs.hcal.index", &idx_hcal);

  //HCal cluster branches
  double hcal_clus_atime[exp_constants::maxclus], hcal_clus_e[exp_constants::maxclus], hcal_clus_x[exp_constants::maxclus], hcal_clus_y[exp_constants::maxclus];
  int hcal_clusid[exp_constants::maxclus];

  C->SetBranchStatus("sbs.hcal.clus.atime",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("sbs.hcal.clus.x",1);
  C->SetBranchStatus("sbs.hcal.clus.y",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("Ndata.sbs.hcal.clus.id", 1);

  C->SetBranchAddress("sbs.hcal.clus.atime", &hcal_clus_atime);
  C->SetBranchAddress("sbs.hcal.clus.e", &hcal_clus_e);
  C->SetBranchAddress("sbs.hcal.clus.x", &hcal_clus_x);
  C->SetBranchAddress("sbs.hcal.clus.y", &hcal_clus_y);
  C->SetBranchAddress("Ndata.sbs.hcal.clus.id", &hcal_clusid);

  //BBCal shower
  double atime_sh, e_sh, nclus_sh, nblk_sh;

  C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.sh.nblk", 1 );

  C->SetBranchAddress("bb.sh.atimeblk", &atime_sh);
  C->SetBranchAddress("bb.sh.e", &e_sh);  
  C->SetBranchAddress("bb.sh.nclus", &nclus_sh);
  C->SetBranchAddress("bb.sh.nblk", &nblk_sh);

  //BBCal preshower
  
  double atime_ps, e_ps, nclus_ps, nblk_ps;

  C->SetBranchStatus( "bb.ps.atimeblk", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nblk", 1 );

  C->SetBranchAddress("bb.ps.atimeblk", &atime_ps);       
  C->SetBranchAddress("bb.ps.e", &e_ps);
  C->SetBranchAddress("bb.ps.nclus", &nclus_ps);
  C->SetBranchAddress("bb.ps.nblk", &nblk_ps);

  //BBGEM hits

  double gem_hits[maxtracks];  

  C->SetBranchStatus("bb.gem.track.nhits", 1);

  C->SetBranchAddress("bb.gem.track.nhits",&gem_hits);

  // track branches

  double ntrack, tr_px[maxtracks], tr_py[maxtracks], tr_pz[maxtracks], tr_p[maxtracks], tr_x[maxtracks], tr_y[maxtracks], tr_vx[maxtracks], tr_vy[maxtracks], tr_vz[maxtracks];
  
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);  
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.x",1);  
  C->SetBranchStatus("bb.tr.y",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.px",&tr_px);
  C->SetBranchAddress("bb.tr.py",&tr_py);
  C->SetBranchAddress("bb.tr.pz",&tr_pz);
  C->SetBranchAddress("bb.tr.p",&tr_p);
  C->SetBranchAddress("bb.tr.x",&tr_x);
  C->SetBranchAddress("bb.tr.y",&tr_y);
  C->SetBranchAddress("bb.tr.vx",&tr_vx);
  C->SetBranchAddress("bb.tr.vy",&tr_vy);
  C->SetBranchAddress("bb.tr.vz",&tr_vz);

  //ekine branches

  double ekine_Q2, ekine_W2, ekine_eps, ekine_nu, ekine_qx, ekine_qy, ekine_qz;

  C->SetBranchStatus("e.kine.Q2",1);
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("e.kine.epsilon",1);
  C->SetBranchStatus("e.kine.nu",1);
  C->SetBranchStatus("e.kine.q_x",1);
  C->SetBranchStatus("e.kine.q_y",1);
  C->SetBranchStatus("e.kine.q_z",1);

  C->SetBranchAddress("e.kine.Q2", &ekine_Q2);
  C->SetBranchAddress("e.kine.W2", &ekine_W2); 
  C->SetBranchAddress("e.kine.epsilon", &ekine_eps);
  C->SetBranchAddress("e.kine.nu", &ekine_nu);
  C->SetBranchAddress("e.kine.q_x", &ekine_qx);
  C->SetBranchAddress("e.kine.q_y", &ekine_qy);
  C->SetBranchAddress("e.kine.q_z", &ekine_qz);

  //global cut branches
  //already handled above

  //setup global cut formula
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  //setup hcal coordinate system with hcal angle with respect to exist beamline
  TVector3 hcal_zaxis = physics::getHCal_zaxis(hcaltheta);
  TVector3 hcal_xaxis = physics::getHCal_xaxis();
  TVector3 hcal_yaxis = physics::getHCal_yaxis(hcal_xaxis,hcal_zaxis);
  //define HCal origin
  TVector3 hcal_origin = physics::getHCal_origin(hcaldist,hcal_offset,hcal_xaxis,hcal_zaxis);  
  double BdL = physics::getBdL(sbs_field);
  //mean energy loss of the beam before scattering
  double Eloss_outgoing = physics::getEloss_outgoing(bbtheta,target);

  //Accounting event variables
  long nevent = 0, nentries = C->GetEntries();

  //ttree formula variables
  int treenum = 0, currenttreenum = 0;

  	//event loop 
  	while(C->GetEntry(nevent++)){
	
	//progress tracker
	cout << "Processing run " <<  j << "/" << num_runs << " run number " << run << " event " << nevent << "/" << nentries << "\r";
	cout.flush();

	//single loop global cut
	currenttreenum = C->GetTreeNumber();
    	if( nevent == 1 || currenttreenum != treenum ){
    	treenum = currenttreenum;
    	GlobalCut->UpdateFormulaLeaves();
    	}
        //Is true if failed global cut
    	bool failedglobal = GlobalCut->EvalInstance(0) == 0;

	///////////
	//Electron-arm physics calculations
	
	//correct beam energy from vertex information
	double Eloss = physics::getEloss(tr_vz[0],target);
	double Ecorr = physics::getEcorr(Ebeam,Eloss);

	//make the vertex, assuming only beam direction
	TVector3 vertex = physics::getVertex(tr_vz[0]);

	//reconstructed momentum, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
	double pcorr = physics::getp_recon_corr(tr_p[0],Eloss_outgoing); 

        //four momentum vector for electron beam with correted Energy value
        TLorentzVector pbeam = physics::getpBeam(Ecorr);

	//four momentum for scattered electron based on reconstruction
	TLorentzVector p_eprime = physics::getp_eprime(tr_px[0],tr_py[0],tr_pz[0],tr_p[0],pcorr);

	}//end event loop
  }//end loop over the run numbers




  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;	

}//end Main
