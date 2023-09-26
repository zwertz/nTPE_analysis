//Script to extract HCal proton and neutron detection efficiency from replayed MC data by Juan Carlos Cornejo detailed by both Provakar and Sebastian
//Ezekiel Wertz Aug. 2023

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

TString proton_root_file, neutron_root_file, Exp,kin, kinematic_name;
double proton_thresh_fac, neutron_thresh_fac;
int num_bin;
double proton_Efit_p0, proton_Efit_p1, neutron_Efit_p0, neutron_Efit_p1, pmin, pmax, Emin, Emax,dx_low,dx_high,dy_low,dy_high;


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
	while(myline.ReadLine(setupfile) && !myline.BeginsWith("endfile")){
         if(!myline.BeginsWith("#")){
         TObjArray *demObjs = myline.Tokenize(" ");
         int numObjs = demObjs->GetEntries();
         	if(numObjs>1){
		TString key = ((TObjString*) (*demObjs)[0])->GetString();
                TString val = ((TObjString*) (*demObjs)[1])->GetString();
			if(key == "proton"){
			proton_root_file = val;
			//cout << "Proton File: " << proton_root_file << endl;		
			}
			else if(key == "neutron"){
			neutron_root_file=val;
			//cout << "Neutron File: " << neutron_root_file << endl;
			}
			else{
                	//We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
			cout << "Found a key that this script can't handle. Fix that!" << endl;
                	return;
                	}
	    }
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
                        //cout << "Experiment: " << Exp << endl;
                        }
                        else if(key == "kin"){
                        kin = val;
                        //cout << "Kinematic: " << kin << endl;
                        }
			else if(key == "proton_thresh_fac"){
			proton_thresh_fac = val.Atoi();
			//cout << "Proton thresh factor: " << proton_thresh_fac << endl;
			}
			else if(key == "neutron_thresh_fac"){
			neutron_thresh_fac = val.Atoi();
			//cout << "Neutron thresh factor: " << neutron_thresh_fac << endl;
			}
			else if(key == "proton_Efit_p0"){
			proton_Efit_p0 = val.Atof();
			//cout << "Proton p0: " << proton_Efit_p0 << endl;
			}
			else if(key == "neutron_Efit_p0"){
			neutron_Efit_p0 = val.Atof();
			//cout << "Neutron p0: " << neutron_Efit_p0 << endl;
			}
			else if(key == "proton_Efit_p1"){
			proton_Efit_p1 = val.Atof();
                        //cout << "Proton p1: " << proton_Efit_p1 << endl;
			}
			else if(key == "neutron_Efit_p1"){
			neutron_Efit_p1 = val.Atof();
                        //cout << "Neutron p1: " << neutron_Efit_p1 << endl;
                        }
			else if(key == "pmin"){
			pmin = val.Atof();
			//cout << "pmin: " << pmin << endl;
                        }
                        else if(key == "pmax"){
			pmax = val.Atof();
			//cout << "pmax: " << pmax << endl;
                        }
                        else if(key == "Emin"){
			Emin = val.Atof();
			//cout << "Emin: " << Emin << endl;
                        }
			else if(key == "Emax"){
			Emax = val.Atof();
			//cout << "Emax: " << Emax << endl;
                        }
                        else if(key == "num_bin"){
			num_bin = val.Atoi();
			//cout << "Num bins: " << num_bin << endl;
                        }
			else if(key == "kinematic_name"){
                        kinematic_name = val;
                        //cout << "Kinematic File " << kinematic__name << endl;
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
}

//Will rely on replayed g4sbs simulated data. Should contain pgun/ngun, zero field, will do for SBS4, SBS8, SBS9
void HCalEfficiency_MC( const char *setup_file_name){


// Define a clock to check macro processing time
 TStopwatch *watch = new TStopwatch();
 watch->Start( kTRUE );


//parse function to get in the information that The One Config file has and is manipulated
parseMainConfig(setup_file_name);


kinematic_obj my_kin(kinematic_name,kin);

//Amount of momentum transferred per bin
double p_step = (pmax-pmin)/num_bin;

//data structure to hold HCal energy histograms. What we will use to determine efficiency
vector<TH1D*> HCalE_proton (num_bin);
vector<TH1D*> HCalE_cut_proton (num_bin);
vector<TH1D*> HCalE_neutron (num_bin);
vector<TH1D*> HCalE_cut_neutron (num_bin);
vector<double> threshold_proton;
vector<double> threshold_neutron;

//Reporting
gStyle->SetOptFit();
gStyle->SetEndErrorSize(0);
TCanvas *c0 = new TCanvas("c0","HCal E vs Proton P",1600,1200);
TCanvas *c1 = new TCanvas("c1","HCal E vs Neutron P",1600,1200);
//Iteration 0 arrays
vector<double> bin_p_pro (num_bin);
vector<double> bin_p_neu (num_bin);

//For proton
vector<double> Emean_p (num_bin);
vector<double> Esig_p (num_bin);
vector<double> binerror_p (num_bin);

//For neutron
vector<double> Emean_n (num_bin);
vector<double> Esig_n (num_bin);
vector<double> binerror_n (num_bin);

//Iteration 1 arrays
vector<double> HCalEff_proton (num_bin);
vector<double> HCalEff_neutron (num_bin);
vector<double> HCalEff_np_ratio (num_bin);


//Will be used for HCal E vs nucleon momentum
auto multi_graph = new TMultiGraph();
//move the fits out here so they can be used anytime
TF1 *proton_gaus_fit, *neutron_gaus_fit;

//setup output file
TString outfile_name = Form("outfiles/mc_HCalEff_%s.root",kin.Data());
TFile *file_out = new TFile(outfile_name,"RECREATE");

//setup relevant histograms
  TH2D *hEdepvP_p = new TH2D("hEdepvP_p","HCal E dep vs proton momentum; p_{proton} (GeV); E_{hcal} (GeV)", num_bin, pmin, pmax, 200, Emin, Emax);
  TH2D *hEdepvP_n = new TH2D("hEdepvP_n","HCal E dep vs neutron momentum; p_{neutron} (GeV); E_{hcal} (GeV)", num_bin, pmin, pmax, 200, Emin, Emax);
  TH2D *hEdepvP_p_Ecut = new TH2D("hEdepvP_p_Ecut","HCal E dep vs proton momentum, Mean E / 4 cut; p_{proton} (GeV); E_{hcal} (GeV)", num_bin, pmin, pmax, 200, Emin, Emax);
  TH2D *hEdepvP_n_Ecut = new TH2D("hEdepvP_n_Ecut","HCal E dep vs neutron momentum, Mean E / 4 cut; p_{neutron} (GeV); E_{hcal} (GeV)", num_bin, pmin, pmax, 200, Emin, Emax);
 // TH1D *hdx_p = new TH1D("hdx_p","dx proton (sd track);x_{HCAL}-x_{expect} (m)", 400, dx_low, dx_high);
 // TH1D *hdy_p = new TH1D("hdy_p","dy proton (sd track);y_{HCAL}-y_{expect} (m)", 400, 3.8, 7.8);
 // TH1D *hdx_n = new TH1D("hdx_n","dx neutron (sd track);x_{HCAL}-x_{expect} (m)", 400, dx_low, dx_high);
 // TH1D *hdy_n = new TH1D("hdy_n","dy neutron (sd track);y_{HCAL}-y_{expect} (m)", 400, 3.8, 7.8);
  TH1D *hdx_p = new TH1D("hdx_p","dx proton (angles);x_{HCAL}-x_{expect} (m)", 400, dx_low, dx_high);
  TH1D *hdy_p = new TH1D("hdy_p","dy proton (angles);y_{HCAL}-y_{expect} (m)", 400, dy_low, dy_high);
  TH1D *hdx_n = new TH1D("hdx_n","dx neutron (angles);x_{HCAL}-x_{expect} (m)", 400, dx_low, dx_high);
  TH1D *hdy_n = new TH1D("hdy_n","dy neutron (angles);y_{HCAL}-y_{expect} (m)", 400, dy_low, dy_high);

  TH1D *hxexp_p = new TH1D("hxexp_p","x exp proton (angles);x_{expect} (m)", 400, dx_low, dx_high);
  TH1D *hyexp_p = new TH1D("hyexp_p","y exp proton (angles);y_{expect} (m)", 400, dy_low, dy_high);
  TH2D *hxyexp_p = new TH2D("hxyexp_p","x exp vs y exp proton; y exp (m); x exp (m)", 400, dy_low, dy_high,400, dx_low, dx_high);
  TH1D *hx_p = new TH1D("hx_p","x_{HCAL} proton; x_{HCAL} (m)", 400, dx_low, dx_high);
  TH1D *hy_p = new TH1D("hy_p","y_{HCAL} proton; y_{HCAL} (m)", 400, dy_low, dy_high);
  TH2D *hxy_p = new TH2D("hxy_p","x_{HCAL} vs y_{HCAL} proton; y_{HCAL} (m); x_{HCAL} (m)", 400, dy_low, dy_high,400, dx_low, dx_high);
  TH2D *hxexp_p_vp = new TH2D("hxexp_p_vp","x exp vs proton p; p_{p} (GeV); x_{expect} (m)", num_bin, pmin, pmax, 400, dx_low, dx_high);
  TH2D *hyexp_p_vp = new TH2D("hyexp_p_vp","y exp vs proton p; p_{p} (GeV); y_{expect} (m)", num_bin, pmin, pmax, 400, dy_low, dy_high);
  TH2D *hx_p_vp = new TH2D("hx_p_vp","x_{HCAL} vs proton p; p_{p} (GeV); x_{HCAL} (m)", num_bin, pmin, pmax, 400, dx_low, dx_high);
  TH2D *hy_p_vp = new TH2D("hy_p_vp","y_{HCAL} vs proton p; p_{p} (GeV); y_{HCAL} (m)", num_bin, pmin, pmax, 400, dy_low, dy_high);

  TH2D *hxexp_n_vp = new TH2D("hxexp_n_vp","x exp vs neutron p; p_{p} (GeV); x_{expect} (m)", num_bin, pmin, pmax, 400, dx_low, dx_high);
  TH2D *hyexp_n_vp = new TH2D("hyexp_n_vp","y exp vs neutron p; p_{p} (GeV); y_{expect} (m)", num_bin, pmin, pmax, 400, dy_low, dy_high);
  TH2D *hx_n_vp = new TH2D("hx_n_vp","x_{HCAL} vs neutron p; p_{p} (GeV); x_{HCAL} (m)", num_bin, pmin, pmax, 400, dx_low, dx_high);
  TH2D *hy_n_vp = new TH2D("hy_n_vp","y_{HCAL} vs neutron p; p_{p} (GeV); y_{HCAL} (m)", num_bin, pmin, pmax, 400, dy_low, dy_high);

  TH1D *hxexp_n = new TH1D("hxexp_n","x exp neutron (angles);x_{expect} (m)", 400, dx_low, dx_high);
  TH1D *hyexp_n = new TH1D("hyexp_n","y exp neutron (angles);y_{expect} (m)", 400, dy_low, dy_high);
  TH2D *hxyexp_n = new TH2D("hxyexp_n","x exp vs y exp neutron; y exp (m); x exp (m)", 400, dy_low, dy_high,400, dx_low, dx_high);
  TH1D *hx_n = new TH1D("hx_n","x_{HCAL} neutron; x_{HCAL} (m)", 400, dx_low, dx_high);
  TH1D *hy_n = new TH1D("hy_n","y_{HCAL} neutron; y_{HCAL} (m)", 400, dy_low, dy_high);
  TH2D *hxy_n = new TH2D("hxy_n","x_{HCAL} vs y_{HCAL} neutron ; y_{HCAL} (m); x_{HCAL} (m)", 400, dy_low, dy_high,400, dx_low, dx_high);

  //TH1D *h_W2 = new TH1D( "W2", "W2 (GeV); GeV", 250, 0.0, 100 );

  TH2D *hdxvp_p = new TH2D("hdxvp_p","dx vs proton p; p_{p} (GeV); x_{HCAL}-x_{expect} (m)", num_bin, pmin, pmax, 400, dx_low, dx_high);
  TH2D *hdyvp_p = new TH2D("hdyvp_p","dy vs proton p; p_{p} (GeV); y_{HCAL}-y_{expect} (m)", num_bin, pmin, pmax, 400, dy_low, dy_high);
  TH2D *hdxvp_n = new TH2D("hdxvp_n","dx vs neutron p; p_{n} (GeV); x_{HCAL}-x_{expect} (m)", num_bin, pmin, pmax, 400, dx_low, dx_high);
  TH2D *hdyvp_n = new TH2D("hdyvp_n","dy vs neutron p; p_{n} (GeV); y_{HCAL}-y_{expect} (m)", num_bin, pmin, pmax, 400, dy_low, dy_high);
 
  //Get relevant kinematic information for HCal
  double hcaltheta = my_kin.getHCalAngle_Rad();
  double hcaldist = my_kin.getHCalDist();
  
  //Define HCal coordinate system
 // TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
 // TVector3 hcal_xaxis(0,-1,0);
 // TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit(); 
 // TVector3 hcal_origin = hcaldist *hcal_zaxis;
  //cout << "HCal Theta:" << hcaltheta << " "<< my_kin.getHCalAngle_Deg() << " "<< (TMath::RadToDeg()*hcaltheta) << endl ;
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  TString nuc;


  //loop over the bins to initialize some histograms
  for(int bin=0; bin< num_bin; bin++){
  HCalE_proton[bin] = new TH1D(Form("HCalE_proton_%i",bin),Form("HCalE_proton_bin_%i",bin),num_bin,Emin,Emax);
  HCalE_cut_proton[bin] = new TH1D(Form("HCalE_cut_proton_%i",bin),Form("HCalE_cut_proton_bin_%i",bin),num_bin,Emin,Emax);

  HCalE_neutron[bin] = new TH1D(Form("HCalE_neutron_%i",bin),Form("HCalE_neutron_bin_%i",bin),num_bin,Emin,Emax);
  HCalE_cut_neutron[bin] = new TH1D(Form("HCalE_cut_neutron_%i",bin),Form("HCalE_cut_neutron_bin_%i",bin),num_bin,Emin,Emax);
  }






  //loop over nucleons, n== 0 is proton, n==1 is neutron
  for (int n=0; n<2; n++) {
  //now intialize TChain
  C = new TChain("T");
  //cout << "How many loops " << n << endl;
  if( n==0 ){
      nuc = "proton";
      //cout << proton_root_file << endl;
      C->Add(proton_root_file);
  }else if( n==1 ){
      nuc = "neutron";
      //cout << neutron_root_file << endl;
      C->Add(neutron_root_file);
  }else{
    //We should never get here. Something has gone horribly wrong with the loop
    cout << "Error: Run loop has a third nucleon" << endl;
    break;
  }
  // setting up ROOT tree branch addresses
  C->SetBranchStatus("*",0);

  // HCal general branches
  double hcalid, ehcal, xhcal, yhcal, row_hcal, col_hcal, hcal_tdc, hcal_atime, kineW2;
  
  C->SetBranchStatus( "sbs.hcal.idblk", 1 );
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus( "sbs.hcal.rowblk", 1 );
  C->SetBranchStatus( "sbs.hcal.colblk", 1 );
  C->SetBranchStatus( "sbs.hcal.atimeblk", 1 );
  C->SetBranchStatus( "sbs.hcal.tdctimeblk", 1 );
  
  C->SetBranchAddress("sbs.hcal.idblk",&hcalid);
  C->SetBranchAddress("sbs.hcal.x",&xhcal);
  C->SetBranchAddress("sbs.hcal.y",&yhcal);
  C->SetBranchAddress("sbs.hcal.e",&ehcal);
  C->SetBranchAddress("sbs.hcal.rowblk",&row_hcal);
  C->SetBranchAddress("sbs.hcal.colblk",&col_hcal);
  C->SetBranchAddress( "sbs.hcal.atimeblk", &hcal_atime );
  C->SetBranchAddress( "sbs.hcal.tdctimeblk", &hcal_tdc );

  //MC nucleon branches
  double mc_p, mc_px, mc_py, mc_pz, mc_vx, mc_vy, mc_vz, mc_nucl, mc_posx, mc_posy, mc_fnucl;
  C->SetBranchStatus("MC.mc_ep",1);
  C->SetBranchStatus("MC.mc_epx",1);
  C->SetBranchStatus("MC.mc_epy",1);
  C->SetBranchStatus("MC.mc_epz",1);
  C->SetBranchStatus("MC.mc_vx",1);
  C->SetBranchStatus("MC.mc_vy",1);
  C->SetBranchStatus("MC.mc_vz",1);
  C->SetBranchStatus("MC.mc_nucl",1);
  C->SetBranchStatus("MC.sdtrack_posx",1);
  C->SetBranchStatus("MC.sdtrack_posy",1);
  C->SetBranchStatus("e.kine.W2",1);

  C->SetBranchAddress("MC.mc_ep", &mc_p);
  C->SetBranchAddress("MC.mc_epx", &mc_px);
  C->SetBranchAddress("MC.mc_epy", &mc_py);
  C->SetBranchAddress("MC.mc_epz", &mc_pz);
  C->SetBranchAddress("MC.mc_vx", &mc_vx);
  C->SetBranchAddress("MC.mc_vy", &mc_vy);
  C->SetBranchAddress("MC.mc_vz", &mc_vz);
  C->SetBranchAddress("MC.mc_nucl", &mc_nucl);
  C->SetBranchAddress("MC.sdtrack_posx", &mc_posx);
  C->SetBranchAddress("MC.sdtrack_posy", &mc_posy);
  C->SetBranchAddress( "e.kine.W2", &kineW2 );

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

  //cout << n << " " << nuc << endl;
  //Loop indices
  long nevent = 0;
  long nevents= C->GetEntries();  
  //cout << "Total events " << nevents << endl;	
 	 while(C->GetEntry(nevent++)){
	//Define HCal coordinate system
	TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
	TVector3 hcal_xaxis(0,-1,0);
	TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
	TVector3 hcal_origin = hcaldist *hcal_zaxis;
	
	//Quick print out to report progress
		if(nevent%100000 == 0){
 		cout << "Processing first loop " << nuc << " MC data, event " << nevent << " / " << nevents << endl; 
		std::cout.flush();
		}
		
		//We got a proton that should have been a neutron
		//if(n==0 &&(((int) (mc_fnucl))==0) ){
		//continue;
		//}
		//We got a nuetron that should have been a proton
		//if(n==1 && (((int) (mc_fnucl))==1)){
		//continue;
		//}

	
  	//Calculate using MC nucleon momenta
  	TVector3 vertex(mc_vx,mc_vy,mc_vz); //Target Vertex vector in Hall Coordinates
  	double theta_nucleon_exp = acos(mc_pz/mc_p); //scattered nucleon theta
	double phi_nucleon_exp = atan2(mc_py,mc_px); //scattered phi nucleon  
	//projected q-vector
	TVector3 pN_hat(sin(theta_nucleon_exp)*cos(phi_nucleon_exp),sin(theta_nucleon_exp)*sin(phi_nucleon_exp) ,cos(theta_nucleon_exp));

	//Intersection of q-vector at plane of HCal
	double s_intersect = (hcal_origin-vertex).Dot(hcal_zaxis)/(pN_hat.Dot(hcal_zaxis));

	TVector3 hcal_intersect = vertex + s_intersect*pN_hat;

	double xhcal_exp = (hcal_intersect - hcal_origin).Dot(hcal_xaxis);
	double yhcal_exp = (hcal_intersect - hcal_origin).Dot(hcal_yaxis);

	double dx = xhcal - xhcal_exp;
	double dy = yhcal - yhcal_exp;
	
  		//Base level acceptance cut 
		if( offhcal ) continue;

		//populate relevant histograms
		if(n==0){//proton
		hEdepvP_p->Fill(mc_p,ehcal);
                hdxvp_p->Fill(mc_p,dx);
                hdyvp_p->Fill(mc_p,dy);
                hdx_p->Fill(dx);
                hdy_p->Fill(dy);

		hxexp_p_vp->Fill(mc_p,xhcal_exp);
		hyexp_p_vp->Fill(mc_p,yhcal_exp);
		hx_p_vp->Fill(mc_p,xhcal);
		hy_p_vp->Fill(mc_p,yhcal);

		hxexp_p->Fill(xhcal_exp);
		hyexp_p->Fill(yhcal_exp);
		hxyexp_p->Fill(yhcal_exp,xhcal_exp);

		hx_p->Fill(xhcal);
		hy_p->Fill(yhcal);
		hxy_p->Fill(yhcal,xhcal);

		double bin_width = hEdepvP_p->GetXaxis()->GetBinWidth(1);
		int bin_i = int((mc_p-pmin)/bin_width);		
		//int bin_i = int((mc_p-pmin)/p_step);
			if(bin_i>num_bin){
			cout << "Warning nucleon momentum bin out of bounds at " << bin_i << endl;
                        continue;
			}
		HCalE_proton[bin_i]->Fill(ehcal);
		}            

		if(n==1){//neutron
		
		hEdepvP_n->Fill(mc_p,ehcal);
                hdxvp_n->Fill(mc_p,dx);
                hdyvp_n->Fill(mc_p,dy);
                hdx_n->Fill(dx);
                hdy_n->Fill(dy);

		hxexp_n->Fill(xhcal_exp);
		hyexp_n->Fill(yhcal_exp);
		hxyexp_n->Fill(yhcal_exp,xhcal_exp);

		hxexp_n_vp->Fill(mc_p,xhcal_exp);
                hyexp_n_vp->Fill(mc_p,yhcal_exp);
                hx_n_vp->Fill(mc_p,xhcal);
                hy_n_vp->Fill(mc_p,yhcal);

		hx_n->Fill(xhcal);
		hy_n->Fill(yhcal);
		hxy_n->Fill(yhcal,xhcal);
		double bin_width = hEdepvP_n->GetXaxis()->GetBinWidth(1);
                int bin_i = int((mc_p-pmin)/bin_width);
		//int bin_i = int((mc_p-pmin)/p_step);
                        if(bin_i>num_bin){
                        cout << "Warning nucleon momentum bin out of bounds at " << bin_i << endl;
                        continue;
                        }

		HCalE_neutron[bin_i]->Fill(ehcal);
		}

  	} //end of event loop

//Calculate the E mean values from scratch. Store that info in files
  if(n==0){//proton
  proton_gaus_fit = new TF1("proton_gaus_fit","gaus");
  //loop over the bins and fit the proton histograms
	for(int bin=0; bin < num_bin; bin++){
        double p = bin*p_step+pmin;
	//Get the energy info from the histogram
 	int max_bin_p = HCalE_proton[bin]->GetMaximumBin();
	double max_bin_center_p = HCalE_proton[bin]->GetXaxis()->GetBinCenter(max_bin_p);
	double max_count_p = HCalE_proton[bin]->GetMaximum();
	double bin_width_p = HCalE_proton[bin]->GetBinWidth(max_bin_p);
	double stdDev_p = HCalE_proton[bin]->GetStdDev();
		//Handle the cases if we have improper data in the first couple of bins. This will mess with fitting the overall histogram
		//Reject lower energy peaks
		if(max_bin_p < 2){
			while(HCalE_proton[bin]->GetBinContent(max_bin_p+1) <= HCalE_proton[bin]->GetBinContent(max_bin_p)){
			cout << HCalE_proton[bin]->GetBinContent(max_bin_p+1) << " "<< HCalE_proton[bin]->GetBinContent(max_bin_p)<<" " << bin << " " << max_bin_p << endl;
			//cout << "Part 3" << endl;
			max_bin_p++;
			}
		HCalE_proton[bin]->GetXaxis()->SetRange(max_bin_p+1, HCalE_proton[bin]->GetNbinsX());
		max_bin_p = HCalE_proton[bin]->GetMaximumBin();
		max_bin_center_p = HCalE_proton[bin]->GetXaxis()->GetBinCenter(max_bin_p);
		max_count_p = HCalE_proton[bin]->GetMaximum();
		bin_width_p = HCalE_proton[bin]->GetBinWidth(max_bin_p);
		stdDev_p = HCalE_proton[bin]->GetStdDev();
		} //end of if condition
	//Fit the histogram
	double lower_bin_p = Emin + max_bin_p*bin_width_p - stdDev_p;
	double upper_bin_p = Emin + max_bin_p*bin_width_p + stdDev_p;
	proton_gaus_fit->SetLineWidth(4);
	proton_gaus_fit->SetParameters(max_count_p,max_bin_center_p);
	proton_gaus_fit->SetRange(lower_bin_p,upper_bin_p);
	HCalE_proton[bin]->Fit(proton_gaus_fit,"RMQ");
	HCalE_proton[bin]->Draw();
        Emean_p[bin] = proton_gaus_fit->GetParameter(1);
        Esig_p[bin]= proton_gaus_fit->GetParameter(2);
        HCalE_proton[bin]->SetTitle(Form("Loop:%d Np:%f Nuc:%d Mean:%f Sigma:%f",bin,p,n,Emean_p[bin],Esig_p[bin]));
	bin_p_pro[bin]=p;
	threshold_proton.push_back(proton_gaus_fit->GetParameter(1)/proton_thresh_fac);
	}// end of bin for loop
  }else if(n==1){//neutron
  neutron_gaus_fit = new TF1("neutron_gaus_fit","gaus");
  //loop over the bins and fit the neutron histograms
        for(int bin=0; bin < num_bin; bin++){
        double p = bin*p_step+pmin;
        //Get the energy info from the histograms
        int max_bin_n = HCalE_neutron[bin]->GetMaximumBin();
        double max_bin_center_n = HCalE_neutron[bin]->GetXaxis()->GetBinCenter(max_bin_n);
        double max_count_n = HCalE_neutron[bin]->GetMaximum();
        double bin_width_n = HCalE_neutron[bin]->GetBinWidth(max_bin_n);
        double stdDev_n = HCalE_neutron[bin]->GetStdDev();

		//Handle the cases if we have improper data in the first couple of bins. This will mess with fitting the overall histogram
		//Reject lower energy peaks
		if(max_bin_n < 2){
			while(HCalE_neutron[bin]->GetBinContent(max_bin_n+1) <= HCalE_neutron[bin]->GetBinContent(max_bin_n)){
                        max_bin_n++;
                        }
                HCalE_neutron[bin]->GetXaxis()->SetRange(max_bin_n+1, HCalE_neutron[bin]->GetNbinsX());
                max_bin_n = HCalE_neutron[bin]->GetMaximumBin();
                max_bin_center_n = HCalE_neutron[bin]->GetXaxis()->GetBinCenter(max_bin_n);
                max_count_n = HCalE_neutron[bin]->GetMaximum();
                bin_width_n = HCalE_neutron[bin]->GetBinWidth(max_bin_n);
        	stdDev_n = HCalE_neutron[bin]->GetStdDev();
            	} //end of if condition
	//Fit the histogram
	double lower_bin_n = Emin + max_bin_n*bin_width_n - stdDev_n;
        double upper_bin_n = Emin + max_bin_n*bin_width_n + stdDev_n;
	neutron_gaus_fit->SetLineWidth(4);
        neutron_gaus_fit->SetParameters(max_count_n,max_bin_center_n);
        neutron_gaus_fit->SetRange(lower_bin_n,upper_bin_n);
        HCalE_neutron[bin]->Fit(neutron_gaus_fit,"RMQ");

        HCalE_neutron[bin]->Draw();
        Emean_n[bin] = neutron_gaus_fit->GetParameter(1);
        Esig_n[bin]= neutron_gaus_fit->GetParameter(2);
        HCalE_neutron[bin]->SetTitle(Form("Loop:%d Np:%f Nuc:%d Mean:%f Sigma:%f",bin,p,n,Emean_n[bin],Esig_n[bin]));
        bin_p_neu[bin]=p;
	threshold_neutron.push_back(neutron_gaus_fit->GetParameter(1)/neutron_thresh_fac);
        }// end of bin for loop
  }else{//We should never get here
  cout << "Error: Run loop has a third nucleon" << endl;
  break;
 }
 
  //Loop a second time to apply thresholds

  //Loop indices
  nevent = 0;
  //cout << "Total events " << nevents << endl; 
  while(C->GetEntry(nevent++)){
  //Quick print out to report progress
  	if(nevent%100000 == 0){
  	cout << "Processing second loop " << nuc << " MC data, event " << nevent << " / " << nevents << endl;
  	std::cout.flush();
  	}
  //We got a proton that should have been a neutron
  if(n==0 && (((int) (mc_fnucl))==0) ){
  continue;
  }
  //We got a nuetron that should have been a proton
  if(n==1 && (((int) (mc_fnucl))==1)){
  continue;
  }
  //Calculate using MC nucleon momenta
  TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
  TVector3 hcal_origin = hcaldist *hcal_zaxis;

  TVector3 vertex(mc_vx,mc_vy,mc_vz); //Target Vertex vector in Hall Coordinates
  double theta_nucleon_exp = acos(mc_pz/mc_p); //scattered nucleon theta
  double phi_nucleon_exp = atan2(mc_py,mc_px); //scattered phi nucleon  
  //projected q-vector
  TVector3 pN_hat(sin(theta_nucleon_exp)*cos(phi_nucleon_exp),sin(theta_nucleon_exp)*sin(phi_nucleon_exp) ,cos(theta_nucleon_exp));

  //Intersection of q-vector at plane of HCal
  double s_intersect = (hcal_origin-vertex).Dot(hcal_zaxis)/(pN_hat.Dot(hcal_zaxis));

  TVector3 hcal_intersect = vertex + s_intersect*pN_hat;

  double xhcal_exp = (hcal_intersect - hcal_origin).Dot(hcal_xaxis);
  double yhcal_exp = (hcal_intersect - hcal_origin).Dot(hcal_yaxis);

  double dx = xhcal - xhcal_exp;
  double dy = yhcal - yhcal_exp;

  //Base level acceptance cut 
  if( offhcal ) continue;
  
  if(n==0){//proton
  double bin_width = hEdepvP_p->GetXaxis()->GetBinWidth(1);
  int bin_i = int((mc_p-pmin)/bin_width);
  
  //Check for energy above threshold
  	if(ehcal > threshold_proton[bin_i]){
	HCalE_cut_proton[bin_i]->Fill(ehcal);
	hEdepvP_p_Ecut->Fill(mc_p,ehcal);
	}

  }//end proton conditional

  if(n==1){//neutron
  double bin_width = hEdepvP_n->GetXaxis()->GetBinWidth(1);
  int bin_i = int((mc_p-pmin)/bin_width);

   //Check for energy above threshold
    if(ehcal > threshold_neutron[bin_i]){
        HCalE_cut_neutron[bin_i]->Fill(ehcal);
        hEdepvP_n_Ecut->Fill(mc_p,ehcal);
        }

  }//end neutron conditional

  } // end of event loop
  
C->Reset();
}//end nucleon loop



 //Draw relevant graphs

 //Convert to arrays because for some reason TGraphErrors does not like vectors 
 double * Emean_p_arr = Emean_p.data(); 
 double * binerror_p_arr = binerror_p.data();
 double * Esig_p_arr = Esig_p.data();
 double * bin_p_pro_arr = bin_p_pro.data();
 c0->cd();
 hEdepvP_p->Draw("colz");

 c0->Update();

 auto gr_p = new TGraphErrors(num_bin,bin_p_pro_arr,Emean_p_arr,binerror_p_arr,Esig_p_arr);  
 gr_p->SetTitle("Proton");
 gr_p->SetMarkerColor(kRed);
 gr_p->SetMarkerStyle(33);
 gr_p->SetMarkerSize(2);
 gr_p->SetLineColor(kRed);
 gr_p->SetLineWidth(2);
 gr_p->Write();
 multi_graph->Add(gr_p);

 gr_p->Draw("P");

 c0->Modified();

 c0->Write();


 //Convert to arrays because for some reason TGraphErrors does not like vectors
 double * Emean_n_arr = Emean_n.data();
 double * binerror_n_arr = binerror_n.data();
 double * Esig_n_arr = Esig_n.data();
 double * bin_p_neu_arr = bin_p_neu.data();
 c1->cd();
 hEdepvP_n->Draw("colz");

 c1->Update();

 auto gr_n = new TGraphErrors(num_bin,bin_p_neu_arr,Emean_n_arr,binerror_n_arr,Esig_n_arr);
 gr_n->SetTitle("Neutron");
 gr_n->SetMarkerColor(kCyan);
 gr_n->SetMarkerStyle(34);
 gr_n->SetMarkerSize(2);
 gr_n->SetLineColor(kCyan);
 gr_n->SetLineWidth(2);
 gr_n->Write();
 multi_graph->Add(gr_n);

 gr_n->Draw("P");
 
 c1->Modified();

 c1->Write();

 multi_graph->SetTitle("HCal E vs Nucleon p");
 multi_graph->GetXaxis()->SetTitle("Nucleon p (GeV)");
 multi_graph->GetYaxis()->SetTitle("E_{hcal}");

double * bin_p_pro_array;
double * bin_p_neu_array;

auto m_graph = new TMultiGraph();

TCanvas *c2 = new TCanvas("c2","HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak})",1600,1200);
TCanvas *c3 = new TCanvas("c3","HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak})",1600,1200);
TCanvas *c4 = new TCanvas("c4","HCal dx Sigma vs proton p (MC)",1600,1200);
TCanvas *c5 = new TCanvas("c5","HCal dx Sigma vs neutron p (MC)",1600,1200);
TCanvas *c6 = new TCanvas("c6","HCal X Res vs Nucleon p (MC)",1600,1200);
TCanvas *c7 = new TCanvas("c7","HCal dy Sigma vs proton p (MC)",1600,1200);
TCanvas *c8 = new TCanvas("c8","HCal dy Sigma vs neutron p (MC)",1600,1200);
TCanvas *c9 = new TCanvas("c9","HCal Y Res vs Nucleon p (MC)",1600,1200);
TCanvas *c10 = new TCanvas("c10","HCal Spatial Resolution (4x4 cluster)",1600,1200);

 //Now we will calculation efficiencies
  for(int b = 0; b < num_bin; b++){
  double p = b*p_step+pmin;
  bin_p_pro[b]=p;
  bin_p_neu[b]=p;
  //cut is numerator. The ones above threshold. Regular is just the energy without any threshold cut. Or all of them
  double proton_eff = ((HCalE_cut_proton[b]->Integral())/(HCalE_proton[b]->Integral()))*100; 
  double neutron_eff = ((HCalE_cut_neutron[b]->Integral())/(HCalE_neutron[b]->Integral()))*100;

  HCalEff_proton[b] = proton_eff;
  HCalEff_neutron[b] = neutron_eff;

  HCalEff_np_ratio[b] = HCalEff_neutron[b]/HCalEff_proton[b];

  }//End for loop to determine efficiency
 
  //Convert to array because graph doesnt like vector
  bin_p_pro_array = bin_p_pro.data();
  bin_p_neu_array = bin_p_neu.data();
  double * HCalEff_neutron_arr = HCalEff_neutron.data();
  double * HCalEff_proton_arr = HCalEff_proton.data();
  c2->cd();
  //Draw graphs
  auto graph_p = new TGraph(num_bin,bin_p_pro_array,HCalEff_proton_arr);
  graph_p->SetTitle("Proton");
  graph_p->SetMarkerColor(kRed);
  graph_p->SetMarkerStyle(20);
  graph_p->SetMarkerSize(1);
  graph_p->SetLineColor(kRed);
  graph_p->SetLineWidth(0);
  m_graph->Add(graph_p);

  auto graph_n = new TGraph(num_bin,bin_p_neu_array,HCalEff_neutron_arr);
  graph_n->SetTitle("Neutron");
  graph_n->SetMarkerColor(kCyan);
  graph_n->SetMarkerStyle(21);
  graph_n->SetMarkerSize(1);
  graph_n->SetLineColor(kCyan);
  graph_n->SetLineWidth(0);
  m_graph->Add(graph_n);

  m_graph->SetTitle(Form("HCAL Efficiency (E_{T}=1/%0.0d E_{Peak}) (4x4 cluster)",proton_thresh_fac));
  m_graph->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");
  m_graph->GetYaxis()->SetTitle("Efficiency (%)");
  m_graph->Draw("AP");

  m_graph->GetYaxis()->SetRangeUser(80.,105.);

  c2->Modified();

  c2->BuildLegend();

  c2->Write();

 c3->SetGrid();

 c3->cd();
  
 double * HCalEff_np_ratio_arr = HCalEff_np_ratio.data();

 auto gr_r = new TGraph(num_bin,bin_p_pro_array,HCalEff_np_ratio_arr);
 gr_r->SetTitle("HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak})");
 gr_r->SetMarkerColor(kMagenta);
 gr_r->SetMarkerStyle(20);
 gr_r->SetMarkerSize(1);
 gr_r->SetLineColor(kMagenta);
 gr_r->SetLineWidth(0);
 gr_r->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");
 gr_r->GetYaxis()->SetTitle("Efficiency Ratio (N/P)");
 gr_r->Draw("AP");
 gr_r->GetYaxis()->SetRangeUser(0.9,1.05);

 c3->Modified();

 c3->Write();

 c4->SetGrid();
 c4->cd();


 vector<double> dbin_p (num_bin);
 vector<double> dxsig_p (num_bin);
 vector<double> dxsig_n (num_bin);
 vector<double> dysig_p (num_bin);
 vector<double> dysig_n (num_bin);

 auto dx_mg = new TMultiGraph();
 auto dy_mg = new TMultiGraph();
 auto all_mg = new TMultiGraph();

//Loop over the nucleons for slices, n=0 proton, n=1 nuetron
  for(int n=0; n<2; n++){
  double p0, p1;

  double fit_l, fit_h;

  vector<TH1D*> pbindx_slice_proton (num_bin);
  vector<TH1D*> pbindy_slice_proton (num_bin);

  vector<TH1D*> pbindx_slice_neutron (num_bin);
  vector<TH1D*> pbindy_slice_neutron (num_bin);
  
	for(int b=0; b<num_bin; b++){
	//Get expected mean from fit to hcale vs nucleon p
	double p =  b*p_step+pmin;
	double fitp1_exp=0;
	double fitp2_exp=0.10;
	
		if(n==0){//proton
		pbindx_slice_proton[b] = hdxvp_p->ProjectionY(Form("pbindxslice_%d_proton",b+1),b+1,b+1);
		pbindy_slice_proton[b] = hdyvp_p->ProjectionY(Form("pbindyslice_%d_proton",b+1),b+1,b+1);
		}else if(n==1){//neutron
		pbindx_slice_neutron[b] = hdxvp_n->ProjectionY(Form("pbindxslice_%d_neutron",b+1),b+1,b+1);
                pbindy_slice_neutron[b] = hdyvp_n->ProjectionY(Form("pbindyslice_%d_neutron",b+1),b+1,b+1);
		}
	fit_l= fitp1_exp - 3*fitp2_exp;
	fit_h= fitp1_exp + 3*fitp2_exp;
	
	TF1 *gausfitdx = new TF1("gausfitdx",gaussian_fit,fit_l,fit_h,3);
	gausfitdx->SetParameter(0,800);
	gausfitdx->SetParameter(1,fitp1_exp);
	gausfitdx->SetParLimits(1,fit_l,fit_h);
	gausfitdx->SetParameter(2,fitp2_exp);
	gausfitdx->SetParLimits(2,0,0.25);

	TF1 *gausfitdy = new TF1("gausfitdy",gaussian_fit,fit_l,fit_h,3);
	gausfitdy->SetParameter(0,800);
	gausfitdy->SetParameter(1,fitp1_exp);
	gausfitdy->SetParLimits(1,fit_l,fit_h);
	gausfitdy->SetParameter(2,fitp2_exp);
	gausfitdy->SetParLimits(2,0,0.25);

	dbin_p[b]=p;

		if(n==0){//proton
		pbindx_slice_proton[b]->Fit("gausfitdx","qRBM");
		pbindx_slice_proton[b]->Draw();
		pbindy_slice_proton[b]->Fit("gausfitdy","qRBM");
		pbindy_slice_proton[b]->Draw();

		dxsig_p[b] = gausfitdx->GetParameter(2); //convert to cm
		dysig_p[b] = gausfitdy->GetParameter(2);
		pbindx_slice_proton[b]->SetTitle(Form("dxslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dxsig_p[b]));
		pbindy_slice_proton[b]->SetTitle(Form("dyslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dysig_p[b]));
		}else if(n==1){//neutron
		pbindx_slice_neutron[b]->Fit("gausfitdx","qRBM");
                pbindx_slice_neutron[b]->Draw();
                pbindy_slice_neutron[b]->Fit("gausfitdy","qRBM");
                pbindy_slice_neutron[b]->Draw();

                dxsig_n[b] = gausfitdx->GetParameter(2); //convert to cm
                dysig_n[b] = gausfitdy->GetParameter(2);
                pbindx_slice_neutron[b]->SetTitle(Form("dxslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dxsig_n[b]));
                pbindy_slice_neutron[b]->SetTitle(Form("dyslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dysig_n[b]));
		}

	}//end bin loop
  }//end nucleon for loop

  double * dbin_p_arr = dbin_p.data();
  double * dxsig_p_arr = dxsig_p.data();

  hdxvp_p->Draw("colz");
  c4->Update();

  //Draw graphs
  auto dxgr_p = new TGraph(num_bin,dbin_p_arr,dxsig_p_arr);
  dxgr_p->SetTitle("Proton X (RMS)");
  dxgr_p->SetMarkerColor(kRed);
  dxgr_p->SetMarkerStyle(21);
  dxgr_p->SetMarkerSize(1);
  dxgr_p->SetLineColor(kRed);
  dxgr_p->SetLineWidth(0);
  dxgr_p->Write();
  dx_mg->Add(dxgr_p);
  all_mg->Add(dxgr_p);

  dxgr_p->Draw("P");

  c4->Modified();
  c4->Write();

  c5->SetGrid();
  c5->cd();

  double * dxsig_n_arr = dxsig_n.data();

  hdxvp_n->Draw("colz");
  c5->Update();

  auto dxgr_n = new TGraph(num_bin,dbin_p_arr,dxsig_n_arr);
  dxgr_n->SetTitle("Neutron X (RMS)");
  dxgr_n->SetMarkerColor(kCyan);
  dxgr_n->SetMarkerStyle(20);
  dxgr_n->SetMarkerSize(1);
  dxgr_n->SetLineColor(kCyan);
  dxgr_n->SetLineWidth(0);
  dxgr_n->Write();
  dx_mg->Add(dxgr_n);
  all_mg->Add(dxgr_n);

  dxgr_n->Draw("P");

  c5->Modified();
  c5->Write();

  c6->SetGrid();
  c6->cd();

  dx_mg->SetTitle("HCal X Res vs Nucleon p (MC)");
  dx_mg->GetXaxis()->SetTitle("Nucleon p (GeV)");
  dx_mg->GetYaxis()->SetTitle("dx sigma (m)");
  dx_mg->Draw("AP");

  c6->BuildLegend();
  c6->Write();

  c7->SetGrid();
  c7->cd();

  double * dysig_p_arr = dysig_p.data();

  hdyvp_p->Draw("colz");
  c7->Update();

  auto dygr_p = new TGraph(num_bin,dbin_p_arr,dysig_p_arr);
  dygr_p->SetTitle("Proton Y (RMS)");
  dygr_p->SetMarkerColor(kRed);
  dygr_p->SetMarkerStyle(25);
  dygr_p->SetMarkerSize(1);
  dygr_p->SetLineColor(kRed);
  dygr_p->SetLineWidth(0);
  dygr_p->Write();
  dy_mg->Add(dygr_p);
  all_mg->Add(dygr_p);

  dygr_p->Draw("P");

  c7->Modified();
  c7->Write();

  c8->SetGrid();
  c8->cd();

  double * dysig_n_arr = dysig_n.data();

  hdyvp_n->Draw("colz");
  c8->Update();

  auto dygr_n = new TGraph(num_bin,dbin_p_arr,dysig_n_arr);
  dygr_n->SetTitle("Neutron Y (RMS)");
  dygr_n->SetMarkerColor(kCyan);
  dygr_n->SetMarkerStyle(24);
  dygr_n->SetMarkerSize(1);
  dygr_n->SetLineColor(kCyan);
  dygr_n->SetLineWidth(0);
  dygr_n->Write();
  dy_mg->Add(dygr_n);
  all_mg->Add(dygr_n);

  dygr_n->Draw("P");

  c8->Modified();
  c8->Write();

  c9->SetGrid();
  c9->cd();

  dy_mg->SetTitle("HCal Y Res vs Nucleon p (MC)");
  dy_mg->GetXaxis()->SetTitle("Nucleon p (GeV)");
  dy_mg->GetYaxis()->SetTitle("dy sigma (m)");
  dy_mg->Draw("AP");

  c9->BuildLegend();
  c9->Write();

  c10->SetGrid();
  c10->cd();

  all_mg->SetTitle("HCal Spatial Resolution (4x4 cluster)");
  all_mg->GetXaxis()->SetTitle("Nucleon momentum (GeV/c)");
  all_mg->GetYaxis()->SetTitle("X and Y Resolution (RMS) (m)");
  all_mg->Draw("AP");

  c10->Modified();

  c10->BuildLegend();

  c10->Write();



file_out->Write();

TString plotname = outfile_name;
plotname.ReplaceAll(".root",".pdf");

c0->Print(plotname + "(");
c1->Print(plotname);
c2->Print(plotname);
c3->Print(plotname);
c4->Print(plotname);
c5->Print(plotname);
c6->Print(plotname);
c7->Print(plotname);
c8->Print(plotname);
c9->Print(plotname);
c10->Print(plotname + ")");

watch->Stop();
// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;
//cout << "Passed the clock" << endl;
}  

