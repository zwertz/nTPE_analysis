//author: Ezekiel Wertz
//short script to to test different scale field values for MC. Created for pass 2 data-simulation peak matching. I did sbs8.

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

//Let's define some global parameters. Mainly what is needed for the config file parser and presumably the rest of the script
TString MC_file, Exp, kin, kinematic_file_name, targ;
TCut globalcut = "";
int SBS_field, MAXNTRACKS, useAlshield;
double W2_low, W2_high, dxO_n,dyO_n,dxsig_n,dysig_n,dxO_p,dyO_p,dxsig_p,dysig_p,dx_low,dx_high,dy_low,dy_high,dxsig_n_fac,dxsig_p_fac,dysig_fac,dxmax,sf;

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
                        if(key == "MC_file"){
                        MC_file = val;
                        //cout << "MC File: " << MC_file << endl;   
                         }
                        else{
                        //We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
                        cout << "Found a key that this script can't handle. Fix that!" << endl;
                        return;
                        }
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
			else if(key == "sf"){
			sf = val.Atof();
			//cout << "Scale Field " << sf << endl;
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
}
void sf_dxdy(const char *setup_file_name){

// Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//Constructor for the data structure we will use to make plots
TChain *C = new TChain("T");

//parse function to get in the information that The One Config file has and is manipulated
parseMainConfig(setup_file_name);

//add the file to the TChain

C->Add(MC_file);
//cout << MC_file << endl;
long Nentries = C->GetEntries();

	if(Nentries < 1){
	cout << "ERROR:File is empty or missing!" << endl;
	return;
	}
//Handles information for kinematics
kinematic_obj my_kin(kinematic_file_name,kin);
//make the outputfile name for later
TString outputfile_name = makeOutputFileName_MC(Exp,kin,SBS_field,sf,targ);



//get useful kinematic information and store it
double Ebeam = my_kin.getBeamEnergy();
double hcaldist = my_kin.getHCalDist();
double hcaltheta = my_kin.getHCalAngle_Rad();
double bbtheta = my_kin.getBBAngle_Rad();
double hcal_offset = getHCalOffset(kin);
double hcalheight = 0;
my_kin.printKinInfo();

//varibles for info from tree branches
double BBtr_px[MAXNTRACKS], BBtr_py[MAXNTRACKS], BBtr_pz[MAXNTRACKS], BBtr_p[MAXNTRACKS];
double vz[MAXNTRACKS];
double BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
double xhcal,yhcal,ehcal,kineW2;
double ntrack, nhits,nucleon;

 // Cut on global parameters from setup config
 TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

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
C->SetBranchStatus( "bb.ps.e", 1 );
C->SetBranchStatus( "bb.ps.x", 1 );
C->SetBranchStatus( "bb.ps.y", 1 );
C->SetBranchStatus( "bb.sh.e", 1 );
C->SetBranchStatus( "bb.sh.x", 1 );
C->SetBranchStatus( "bb.sh.y", 1 );
C->SetBranchStatus("bb.gem.track.nhits", 1);
C->SetBranchStatus("e.kine.W2",1);
C->SetBranchStatus("MC.mc_nucl",1);
C->SetBranchStatus( "sbs.hcal.nclus", 1 );

C->SetBranchAddress("bb.tr.vz",vz);
C->SetBranchAddress("bb.tr.px",BBtr_px);
C->SetBranchAddress("bb.tr.py",BBtr_py);
C->SetBranchAddress("bb.tr.pz",BBtr_pz);
C->SetBranchAddress("bb.tr.p",BBtr_p);
C->SetBranchAddress("sbs.hcal.x",&xhcal);
C->SetBranchAddress("sbs.hcal.y",&yhcal);
C->SetBranchAddress("sbs.hcal.e",&ehcal);
C->SetBranchAddress( "bb.ps.e", &BBps_e );
C->SetBranchAddress( "bb.ps.x", &BBps_x );
C->SetBranchAddress( "bb.ps.y", &BBps_y );
C->SetBranchAddress( "bb.sh.e", &BBsh_e );
C->SetBranchAddress( "bb.sh.x", &BBsh_x );
C->SetBranchAddress( "bb.sh.y", &BBsh_y );
C->SetBranchAddress( "e.kine.W2", &kineW2 );
C->SetBranchAddress("bb.tr.n",&ntrack);
C->SetBranchAddress("bb.gem.track.nhits",&nhits);
C->SetBranchAddress("MC.mc_nucl", &nucleon);

//output file construction
TFile *fout = new TFile(outputfile_name,"RECREATE");

//Let's make some useful histograms
TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal X  vs Y, no cuts;HCal Y  (m); HCal X  (m)", 300, -2.0, 2.0, 500, -2.5, 2.5 );
TH2D *hxy_fidcut = new TH2D("hxy_fidcut","HCal X vs Y, fiducial cut;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -3.0, 3.0);
TH2D *hxy_expect_nocut = new TH2D("hxy_expect_nocut","HCal X Expect vs Y Expect, no cuts;HCal Y Expect (m); HCal X Expect (m)", 300, -2.5, 2.5, 500, -2.5, 2.5 );
TH2D *hxy_expect_fidcut = new TH2D("hxy_expect_fidcut","HCal X Expect vs Y Expect, fid cuts;HCal Y Expect (m); HCal X Expect (m)", 300, -2.5, 2.5, 500, -2.5, 2.5 );
TH1D *hX = new TH1D( "X", "HCal X (m); m", 250, dx_low, dx_high );
TH1D *hX_expect = new TH1D( "X_expect", "HCal X Expect (m); m", 250, dx_low, dx_high+1 );

TH1D *hdx = new TH1D( "dx", "HCal dx (m); m", 250, dx_low, dx_high );
TH1D *hdx_cut = new TH1D( "dx_cut", "HCal dx (m),  cuts; m", 250, dx_low, dx_high );
TH1D *hdx_fidcut = new TH1D( "dx_fidcut","HCal dx fiducial cut; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
TH1D *hdy = new TH1D( "dy", "HCal dy (m); m", 250, dy_low, dy_high );
TH1D *hdy_fidcut = new TH1D( "dy_fidcut","HCal dy fiducial cut; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );
TH1D *hY = new TH1D( "Y", "HCal Y (m); m", 250, dy_low, dy_high );
TH1D *hY_expect = new TH1D( "Y_expect", "HCal Y Expect (m); m", 250, dy_low, dy_high );
TH1D *h_W2 = new TH1D( "W2", "W2 (GeV) No Cuts; GeV", 250, -1.0, 6.0 );
TH1D *hQ2 = new TH1D("Q2","Q2",250,0.5,6.0);

TH1D *hdx_p = new TH1D( "dx_p","HCal dx; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
TH1D *hdx_n = new TH1D( "dx_n","HCal dx; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
TH1D *hdx_pcut = new TH1D( "dx_pcut","HCal dx just pcut; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );
TH1D *hdx_ncut = new TH1D( "dx_ncut","HCal dx just ncut; x_{HCAL}-x_{expect} (m)", 250, dx_low, dx_high );

TH1D *h_vert_z = new TH1D( "vert_z", "Vertex Position z-direction (m); m", 200, -0.2, 0.2 );
double est_peak_loc = -0.01342*sf-0.004486;

cout << "Peak loc: "<< est_peak_loc << endl;
//loop over the events in the TChain
long nevent = 0;
int treenum = 0, currenttreenum = 0,currentrun =0;
while(C->GetEntry(nevent++)){
currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
    //  cout << "My leaves" << endl;
    treenum = currenttreenum;
    GlobalCut->UpdateFormulaLeaves();
   }
bool failedglobal = GlobalCut->EvalInstance(0) == 0;

double p_ep = BBtr_p[0]; // Obtain the magnitude of scattered electron momentum
double etheta = acos( BBtr_pz[0]/p_ep); //Use the uncorrected track momentum to reconstruct e' thetha
double ephi = atan2( BBtr_py[0], BBtr_px[0] );
TVector3 vertex( 0, 0, vz[0] ); //z location of vertex in hall coordinates
TLorentzVector Pbeam(0,0,Ebeam,Ebeam);//Mass of e negligable
TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
TLorentzVector Ptarg( 0, 0, 0, M_p ); //Just use proton?
TLorentzVector q = Pbeam - kprime; //Standard q-vector
TLorentzVector PgammaN = Ptarg+q; // (-px, -py, Ebeam-pz,Mp+Ebeam-p)
double pelastic = Ebeam /(1.0+(Ebeam/M_p)*(1.0-cos(etheta)));
double E_ep = sqrt(pow(M_e,2)+pow(BBtr_p[0],2)); //Obtain the scattered electron energy
double nu = Ebeam-E_ep; // Obtain energy transfer
double pp = sqrt(pow(nu,2)+ 2*M_p*nu);
double phinucleon = ephi + TMath::Pi(); //assume coplanarity
double thetanucleon = acos( (Ebeam - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics

TVector3 pNhat( sin(thetanucleon)*cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));

//Define HCal coordinate systemof event while loop
TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
TVector3 hcal_xaxis(0,-1,0);
TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
TVector3 hcal_origin = hcaldist *hcal_zaxis+ hcalheight *hcal_xaxis;
TVector3 hcalpos = hcal_origin + xhcal * hcal_xaxis + yhcal * hcal_yaxis;
//Define interesection points for hadron vector
double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/pNhat.Dot( hcal_zaxis );
TVector3 hcal_intersect = vertex + sintersect * pNhat;
//Define the expected position of hadron on HCal from BB track. This assumes a neutron hypothesis since the expect values for the HCal positions do not account for the deflection of the protons due to the magnetic field.
double xhcal_expect = (hcal_intersect-hcal_origin).Dot( hcal_xaxis );
double yhcal_expect = (hcal_intersect-hcal_origin).Dot( hcal_yaxis );
double W2 = kineW2;//Get the invariant mass transfer W and four-momentum of the scattered nucleon
double Q2 = 2*Ebeam*E_ep*( 1-cos(etheta) );


//define dx,dy
double dx = xhcal - xhcal_expect;
double dy = yhcal - yhcal_expect;

hX->Fill(xhcal);
hX_expect->Fill(xhcal_expect);
hxy_expect_nocut->Fill(yhcal_expect,xhcal_expect);
hY->Fill(yhcal);
hY_expect->Fill(yhcal_expect);
h_W2->Fill(W2);
hQ2->Fill(Q2);
h_vert_z->Fill(vz[0]);

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
//if( offhcal ) continue;


//Define Fiducial cut
double HCal_fid_top = HCal_top - dxsig_p;
double HCal_fid_bot = HCal_bot + dxsig_n;
double HCal_fid_left = HCal_left + dysig_p; //currently dysig is same for both p and n
double HCal_fid_right = HCal_right - dysig_p;
//Actual fiducial cut. This assumes a neutron and proton hypothesis since the expect values for the HCal positions do not account for the deflection of the protons due to the magnetic field.

bool neutron_hyp_fid = ((xhcal_expect - est_peak_loc) >= HCal_fid_top) && xhcal_expect <= HCal_bot ;
bool proton_hyp_fid = ((xhcal_expect + est_peak_loc) <= HCal_fid_bot) && xhcal_expect >= HCal_top ;

bool inFiducial = (neutron_hyp_fid || proton_hyp_fid);

hdy->Fill(dy);
hdx->Fill(dx);
hxy_nocut->Fill(yhcal,xhcal);

//neutrons
if(nucleon==0){
hdx_n->Fill(dx);
}
//protons
if(nucleon==1){
hdx_p->Fill(dx);
}

bool bad_dy = (abs(dy-dyO_n) > (dysig_fac*dysig_n)) || (abs(dy-dyO_p)>(dysig_fac*dysig_p));

if(bad_dy){
continue;
}


if(!failedglobal && !offhcal && inFiducial){
hdy_fidcut->Fill(dy);
hdx_fidcut->Fill(dx);
hxy_fidcut->Fill(yhcal,xhcal);
hxy_expect_fidcut->Fill(yhcal_expect,xhcal_expect);
	//neutrons
	if(nucleon==0){
	hdx_ncut->Fill(dx);
	}
	//protons
	if(nucleon==1){
	hdx_pcut->Fill(dx);
	}
}



hdx_cut->Fill(dx);

}//end of event loop

fout->Write();

watch->Stop();

// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;
}
