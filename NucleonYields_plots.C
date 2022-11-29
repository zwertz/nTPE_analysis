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


void NucleonYields_plots( double runnum, const char *data_file_name, const char *kinematic_file_name, double sbsfieldscale=1.0){
//Convert from double to char and then to TString. This is just for book keeping
 stringstream ss;
 ss << runnum;
 std::string runnum_temp;
 ss >> runnum_temp;
 const char *runnum_char = runnum_temp.c_str();
 TString runnum_string(runnum_char);
 //cout << runnum_char << endl;
 //The next couple blocks of code use a iostream to read-in data from a file to properly constructed the location of the root file data. So it figures out the path for you and all you need to remember is the run number. 
 ifstream datafile(data_file_name);
 TString currentLine;
 TString datRun,pass,kinematic,target;
 bool gotRun = false;
 while(currentLine.ReadLine(datafile)){
 if(currentLine.BeginsWith("#")){
  //Treat # as comments and ignore them
  //cout << "Cond 1" << endl;
    continue;
  }
 else if(currentLine.BeginsWith(runnum_string)){
 //We found the right run number. Store the info
  TObjArray *tokens = currentLine.Tokenize(" ");
  //Assuming ordering is runnum, pass, SBSKin, Target, sbsfieldscale
  datRun =  ((TObjString*) (*tokens)[0])->GetString();
  pass = ((TObjString*) (*tokens)[1])->GetString();
  kinematic =  ((TObjString*) (*tokens)[2])->GetString();
  target =  ((TObjString*) (*tokens)[3])->GetString();

  cout << "The run " << datRun << " The Pass " << pass << " The Kin " << kinematic << " The Target " << target << endl;
  gotRun = true;
  //cout << gotRun << endl;
 }

 else{
  //Where are still searching but it's not a comment
  //  //cout << "Cond 3" << endl; 
  continue;
  }
 }

  if ((datafile.eof()) && !gotRun){
  //Conditional that we checked the entire data file and did not find the runnum
  cout << "Did not find run number: " << runnum_char << " in the data file! Quitting, figure it out!" << endl;
  return;
 }
 //All of this was to have a modular input directory. So let's make it
 string input_directory = "/work/halla/sbs/sbs-gmn";
  const char *input_directory_char = input_directory.c_str();
  const char *datRun_char = datRun.Data();
  const char *pass_char = pass.Data();
  const char *kin_char = kinematic.Data();
  const char *tar_char = target.Data();
 TString inputfile = Form("%s/%s/%s/%s/rootfiles/e1209019_fullreplay_%s_*.root",input_directory_char,pass_char,kin_char,tar_char,datRun_char);
 //cout << "File Location " << inputfile << endl;



//Read-in loop and conditionals to find the correct kinematic information so we don't have to remember it every time we want to call the script.
double Ebeam,bbtheta,sbstheta,sbsdist,hcaldist;

ifstream kinfile(kinematic_file_name);
 TString datLine;
 TString datKin;
 bool gotKin = false;
 while(datLine.ReadLine(kinfile)){
 if(currentLine.BeginsWith("#")){
  //Treat # as comments and ignore them
  //cout << "Cond 1" << endl;
   continue;
  }
 else if(datLine.BeginsWith(kinematic)){
 //We found the right kinematic. Store the info
 TObjArray *myobjs = datLine.Tokenize(" ");
  //Assuming ordering is kinematic, Beam Energy, BB Angle, SBS angle, SBS dist, HCal dist
  datKin = ((TObjString*) (*myobjs)[0])->GetString();
  Ebeam = (((TObjString*) (*myobjs)[1])->GetString()).Atof();
  bbtheta = (((TObjString*) (*myobjs)[2])->GetString()).Atof();
  sbstheta = (((TObjString*) (*myobjs)[3])->GetString()).Atof();
  sbsdist = (((TObjString*) (*myobjs)[4])->GetString()).Atof();
  hcaldist = (((TObjString*) (*myobjs)[5])->GetString()).Atof();
  cout << "Kinematic " << datKin << "  Beam Energy " << Ebeam << " BB Angle  " << bbtheta << " SBS Angle  " << sbstheta << " SBS Dist  " << sbsdist << " HCal Dist  " << hcaldist  << endl; 
  gotKin = true;
   //cout << gotRun << endl;
   }
 else{
  //Where are still searching but it's not a comment
  //cout << "Cond 3" << endl; 
  continue;
  }
 }
  
 if ((kinfile.eof()) && !gotKin){
  //Conditional that we checked the entire kinematic file and did not find the kinematic info
 cout << "Did not find kinematic: " << datKin << " in the kinematic file! Quitting, figure it out!" << endl;
  return;
 }






//To convert from degrees to radians
  sbstheta *= TMath::Pi()/180.0;
  bbtheta *= TMath::Pi()/180.0;

//Looks like this is assuming it is taking one root file
  TChain *C = new TChain("T");

  C->Add(inputfile);
 // C->Add(rootfilename); 
//need to understand what these variables mean
  TCut globalcut = "bb.ps.e>0.15&&abs(bb.tr.vz)<0.27&&sbs.hcal.nclus>0&&bb.tr.n==1";
  
  TEventList *elist = new TEventList("elist","");
  
  C->Draw(">>elist",globalcut);

  int MAXNTRACKS=10;
  
  //variables we need are BigBite track px,py,pz and sbs hcal x, y, e

  double ntrack;

  double epx[MAXNTRACKS];
  double epy[MAXNTRACKS];
  double epz[MAXNTRACKS];
  double ep[MAXNTRACKS];

  double vx[MAXNTRACKS];
  double vy[MAXNTRACKS];
  double vz[MAXNTRACKS];

  double xhcal,yhcal,ehcal;

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

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.tr.px",epx);
  C->SetBranchAddress("bb.tr.py",epy);
  C->SetBranchAddress("bb.tr.pz",epz);
  C->SetBranchAddress("bb.tr.p",ep);
  C->SetBranchAddress("sbs.hcal.x",&xhcal);
  C->SetBranchAddress("sbs.hcal.y",&yhcal);
  C->SetBranchAddress("sbs.hcal.e",&ehcal);

  TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
  TLorentzVector Ptarg(0,0,0,0.5*(0.938272+0.939565));

  long nevent=0;

  double W2min = 0.88-0.4;
  double W2max = 0.88+0.4;
//Unique name structure for output fules
 
 TString outfile = Form("elastic_output/elastic_Zeke_%s_%s_%s_%s.root",pass_char,kin_char,tar_char,datRun_char);
 
  TFile *fout = new TFile(outfile,"RECREATE");

  TH2D *hdxdy_all = new TH2D("hdxdy_all","All events;#Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wcut = new TH2D("hdxdy_Wcut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wanticut = new TH2D("hdxdy_Wanticut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);

  TH1D *hW2_all = new TH1D("hW2_all","All events; W^{2} (GeV^{2});", 250, -1, 4 );

  TH1D *hEHCAL_all = new TH1D("hEHCAL_all","All events; HCAL energy sum (GeV);",250,0,0.5);
  TH1D *hEHCAL_Wcut = new TH1D("hEHCAL_Wcut","|W^{2}-0.88|<0.4;HCAL energy sum (GeV);",250,0,0.5);

  TVector3 hcal_origin( -hcaldist*sin(sbstheta), 0, hcaldist*cos(sbstheta) );

  TVector3 hcal_zaxis = hcal_origin.Unit();
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
  
  while( C->GetEntry(elist->GetEntry(nevent++) ) ){
    if( ntrack == 1.0 ){
      TLorentzVector kprime( epx[0], epy[0], epz[0], ep[0] );
      TLorentzVector q = Pbeam - kprime;

      TVector3 qdir = q.Vect().Unit();

      TVector3 vertex(0,0,vz[0]);

      double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qdir.Dot( hcal_zaxis );

      TVector3 hcal_intersect = vertex + sintersect * qdir; 

      double xhcal_expect = hcal_intersect.Dot( hcal_xaxis );
      double yhcal_expect = hcal_intersect.Dot( hcal_yaxis );

      hdxdy_all->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
      double W2recon = (Ptarg + q).M2();

      hW2_all->Fill( W2recon );

      hEHCAL_all->Fill( ehcal );

      if( W2recon > W2min && W2recon < W2max ){
	hdxdy_Wcut->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
	hEHCAL_Wcut->Fill( ehcal );
      } else {
	hdxdy_Wanticut->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
      }
    }
  }

  elist->Delete(); 
  fout->Write();

} 