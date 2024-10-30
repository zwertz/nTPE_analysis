//Author: Ezekiel Wertz
//10/23/2024
//Purpose: A script that will take MC information and consider position dependent HCal efficiency. Currently it is only compatible with single magnetic field settings.
//This is based off of scripts intially developed by A. Puckett and P.Datta

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
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/fits.C"
#include "../../src/data_object.C"
#include "../../src/cuts.C"
#include "../../src/physics.C"
#include "../../src/parse_config.C"
#include "../../src/plots.C"
#include "../../src/calc_FFs_RCS_obj.C"
#include "../../src/fit_histogram.C"

//Main
void HCal_p_Eff_uniformity_MC(const char *setup_file_name){
//Set this default to true so that way fits to histogram should be more correct. This effects statistical error
TH1::SetDefaultSumw2(kTRUE);

//Define a clock to check macro processing time
  TStopwatch *watch = new TStopwatch();
  watch->Start( kTRUE );

  //gStyle->SetErrorX(0);

  //parse object to get in the information that The One Config file has and is manipulated
  parse_config mainConfig(setup_file_name);
  //Need to update print statement for like HCal Efficiency. Revisit once script mostly determined
  //mainConfig.printDataYields();
  
  //store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
  TCut globalcut = mainConfig.getGlobalCut();
  TString exp = mainConfig.getExp();
  TString kin = mainConfig.getKin();
  TString pass = mainConfig.getPass();
  int sbs_field = mainConfig.getSBSField();
  TString target = mainConfig.getTarg();
  TString kinematic_file = mainConfig.getKinFileName();
  TString MC_input_file_name = mainConfig.getMCFileName();
  double ehcal_min = mainConfig.getHCaleMin(); 
  double dxO_p = mainConfig.get_dxOp();
  double dyO_p = mainConfig.get_dyOp();
  double dxsig_p = mainConfig.get_dxsigp();
  double dysig_p = mainConfig.get_dysigp();
  double W2fitmax = mainConfig.getW2FitMax();
  double binfac = mainConfig.getBinFac();
  double fidx_min = mainConfig.get_fidxmin();
  double fidx_max = mainConfig.get_fidxmax();
  double fidy_min = mainConfig.get_fidymin();
  double fidy_max = mainConfig.get_fidymax();
  double W2_low = mainConfig.getW2Low();
  double W2_high = mainConfig.getW2High();
  double spot_sig = mainConfig.get_spotsig();
  double fitx_low = mainConfig.get_fitxlow();
  double fitx_high = mainConfig.get_fitxhigh();
  double fity_low = mainConfig.get_fitylow();
  double fity_high = mainConfig.get_fityhigh();


  //setup hcal active area with bounds that match database depending on pass
  vector<double> hcalaa = cuts::hcal_Position_MC();

  TFile *mc_file = new TFile(MC_input_file_name.Data());

  //setup output file
  TString outfile = utility::makeOutputFileNameHCalEffUniformity(exp,pass,kin,sbs_field,target);
  TFile *fout = new TFile(outfile,"RECREATE");

  //Histograms we will fill
  TH1D *h_W2_globcut = new TH1D( "hW2_globcut", "W2 (GeV^{2}) global cut and fiducial cut; GeV^{2}", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_hcalcut = new TH1D( "hW2_hcalcut", "W2 (GeV^{2}) global+fiducial+HCal cut; GeV^{2}", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_anticut = new TH1D( "hW2_anticut", "W2 (GeV^{2}) global+fiducial+!HCal cut; GeV^{2}", binfac*W2fitmax, 0.0, W2fitmax );

  TH1D *hx_expect_all = new TH1D("hx_expect_all","HCal X Expect, global cut;HCal X Expect (m)", 147, -3.5, 3.0 );
  TH1D *hx_expect_hcalcut = new TH1D("hx_expect_hcalcut","HCal X Expect, global+HCal cut;HCal X Expect (m)", 147, -3.5, 3.0 );
  TH1D *hx_expect_anticut = new TH1D("hx_expect_anticut","HCal X Expect, global+ !HCal cut;HCal X Expect (m)", 147, -3.5, 3.0 );

  TH1D *hy_expect_all = new TH1D("hy_expect_all","HCal Y Expect, global cut;HCal Y Expect (m)", 67, -1.25, 1.25 );
  TH1D *hy_expect_hcalcut = new TH1D("hy_expect_hcalcut","HCal Y Expect, global+HCal cut;HCal Y Expect (m)", 67, -1.25, 1.25 );
  TH1D *hy_expect_anticut = new TH1D("hy_expect_anticut","HCal Y Expect, global+ !HCal cut;HCal Y Expect (m)", 67, -1.25, 1.25 );

  TH2D *hxy_expect_all = new TH2D("hxy_expect_all","HCal X Expect vs Y Expect, global cut;HCal Y Expect (m); HCal X Expect (m)", 67, -1.25, 1.25, 147, -3.5, 2.0 );
  TH2D *hxy_expect_hcalcut = new TH2D("hxy_expect_hcalcut","HCal X Expect vs Y Expect, global cut+HCal cut;HCal Y Expect (m); HCal X Expect (m)", 67, -1.25, 1.25, 147, -3.5, 2.0 );
  TH2D *hxy_expect_anticut = new TH2D("hxy_expect_anticut","HCal X Expect vs Y Expect, global cut+!HCal cut;HCal Y Expect (m); HCal X Expect (m)", 67, -1.25, 1.25, 147, -3.5, 2.0 );

  TH1D *hx_HCal_all = new TH1D("hx_HCal_all","HCal X observed, global cut; HCal X  (m)", 147, -3.5, 2.0 );
  TH1D *hy_HCal_all = new TH1D("hy_HCal_all","HCal Y observed, global cut; HCal Y  (m)", 67, -1.25, 1.25 );
  TH2D *hxy_HCal_all = new TH2D("hxy_HCal_all", "HCal X vs Y observed, global cut;HCal Y  (m),HCal X  (m)",67, -1.25, 1.25,147, -3.5, 2.0);

  TH1D *hx_HCal_cut = new TH1D("hx_HCal_cut","HCal X observed, global cut+HCal cut; HCal X  (m)", 147, -3.5, 2.0 );
  TH1D *hy_HCal_cut = new TH1D("hy_HCal_cut","HCal Y observed, global cut+HCal cut; HCal Y  (m)", 67, -1.25, 1.25 );
  TH2D *hxy_HCal_cut = new TH2D("hxy_HCal_cut", "HCal X vs Y observed, global cut+HCal cut;HCal Y  (m),HCal X  (m)",67, -1.25, 1.25,147, -3.5, 2.0);

  TH1D *hx_HCal_anticut = new TH1D("hx_HCal_anticut","HCal X observed, global cut+!HCal cut; HCal X  (m)", 147, -3.5, 2.0 );
  TH1D *hy_HCal_anticut = new TH1D("hy_HCal_anticut","HCal Y observed, global cut+!HCal cut; HCal Y  (m)", 67, -1.25, 1.25 );
  TH2D *hxy_HCal_anticut = new TH2D("hxy_HCal_anticut", "HCal X vs Y observed, global cut+!HCal cut;HCal Y  (m),HCal X  (m)",67, -1.25, 1.25,147, -3.5, 2.0);

  //Disabling row column info for now. Not sure if essentially. Would need to implement info first in parser to enable here
  
  TH1D *hrowHCal_all = new TH1D("hrowHCal_all","HCal Row, global cut; HCal row",24,-0.5,23.5);
  TH1D *hcolHCal_all = new TH1D("hcolHCal_all","HCal Col, global cut; HCal col",12,-0.5,11.5);
  TH2D *hrowcolHCal_all = new TH2D("hrowcolHCal_all","HCal Row vs Col, global cut;HCal col;HCal row",12,-0.5,11.5,24,-0.5,23.5);

  TH1D *hrowHCal_cut = new TH1D("hrowHCal_cut","HCal Row, global+HCal cut; HCal row",24,-0.5,23.5);
  TH1D *hcolHCal_cut = new TH1D("hcolHCal_cut","HCal Col, global+HCAl cut; HCal col",12,-0.5,11.5);
  TH2D *hrowcolHCal_cut = new TH2D("hrowcolHCal_cut","HCal Row vs Col, global+HCal cut;HCal col;HCal row",12,-0.5,11.5,24,-0.5,23.5);

  TH1D *hrowHCal_anticut = new TH1D("hrowHCal_anticut","HCal Row, global+!HCal cut; HCal row",24,-0.5,23.5);
  TH1D *hcolHCal_anticut = new TH1D("hcolHCal_anticut","HCal Col, global+!HCAl cut; HCal col",12,-0.5,11.5);
  TH2D *hrowcolHCal_anticut = new TH2D("hrowcolHCal_anticut","HCal Row vs Col, global+!HCal cut;HCal col;HCal row",12,-0.5,11.5,24,-0.5,23.5);

  TH1D *hdx_all = new TH1D("hdx_all","HCal dx, global cuts; dx (m)",250,-3,2);
  TH1D *hdx_cut = new TH1D("hdx_cut","HCal dx, global+HCal cuts; dx (m)",250,-3,2);
  TH1D *hdx_anticut = new TH1D("hdx_anticut","HCal dx, global+!HCal cuts; dx (m)",250,-3,2);

  TH1D *hdy_all = new TH1D("hdy_all","HCal dy, global cuts; dy (m)",250,-1.25,1.25);
  TH1D *hdy_cut = new TH1D("hdy_cut","HCal dy, global+HCal cuts; dy (m)",250,-1.25,1.25);
  TH1D *hdy_anticut = new TH1D("hdy_anticut","HCal dy, global+!HCal cuts; dy (m)",250,-1.25,1.25);

  TH2D *hdxdy_all = new TH2D("hdxdy_all","HCal dx vs dy, global cuts; dy (m); dx (m)",250,-1.25,1.25,250,-2,1.25);
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","HCal dx vs dy, global+HCal cuts; dy (m); dx (m)",250,-1.25,1.25,250,-2,1.25);
  TH2D *hdxdy_anticut = new TH2D("hdxdy_anticut","HCal dx vs dy, global+!HCal cuts; dy (m); dx (m)",250,-1.25,1.25,250,-2,1.25);

  TH1D *hehcal_all = new TH1D("hehcal_all","E_{HCal}, global cuts; E_{HCal} (GeV)",500,0,2.5);
  TH1D *hehcal_cut = new TH1D("hehcal_cut","E_{HCal}, global+HCal cuts; E_{HCal} (GeV)",500,0,2.5);
  TH1D *hehcal_anticut = new TH1D("hehcal_anticut","E_{HCal}, global+!HCal cuts; E_{HCal} (GeV)",500,0,2.5);

  //Intialize the TChain for data
  TChain *C = new TChain("Parse");
  C->Add(MC_input_file_name.Data());

  //Branch variables
  double ehcal, xhcal, yhcal, dx, dy, proton_deflection, xexp, yexp, W2, BBps_e, BBtr_vz, BBtr_p, BBgem_nhits, rowblkHCAL, colblkHCAL, final_mc_weight, BBgem_chi2ndf;
  int BBsh_nclus;

  //Get the branches we need from Parse
  C->SetBranchStatus("ehcal",1);
  C->SetBranchStatus("xhcal",1);
  C->SetBranchStatus("yhcal",1);
  C->SetBranchStatus("xexp",1);
  C->SetBranchStatus("yexp",1);
  C->SetBranchStatus("dx",1);
  C->SetBranchStatus("dy",1);
  C->SetBranchStatus("proton_deflection",1);
  C->SetBranchStatus("W2",1);
  C->SetBranchStatus("BBps_e",1);
  C->SetBranchStatus("BBtr_vz",1);
  C->SetBranchStatus("BBtr_vz",1);
  C->SetBranchStatus("BBgem_nhits",1);
  C->SetBranchStatus("BBsh_nclus",1);
  C->SetBranchStatus("rowblkHCAL",1);
  C->SetBranchStatus("colblkHCAL",1);
  C->SetBranchStatus("Final_MC_weight",1);
  C->SetBranchStatus("BBgem_chi2ndf",1);

  C->SetBranchAddress("ehcal",&ehcal);
  C->SetBranchAddress("xhcal",&xhcal);
  C->SetBranchAddress("yhcal",&yhcal);
  C->SetBranchAddress("xexp",&xexp);
  C->SetBranchAddress("yexp",&yexp);
  C->SetBranchAddress("dx",&dx);
  C->SetBranchAddress("dy",&dy);
  C->SetBranchAddress("proton_deflection",&proton_deflection);
  C->SetBranchAddress("W2",&W2);
  C->SetBranchAddress("BBps_e",&BBps_e);
  C->SetBranchAddress("BBtr_vz",&BBtr_vz);
  C->SetBranchAddress("BBtr_p",&BBtr_p);
  C->SetBranchAddress("BBgem_nhits",&BBgem_nhits);
  C->SetBranchAddress("BBsh_nclus",&BBsh_nclus);
  C->SetBranchAddress("rowblkHCAL",&rowblkHCAL);
  C->SetBranchAddress("colblkHCAL",&colblkHCAL);
  C->SetBranchAddress("Final_MC_weight",&final_mc_weight);
  C->SetBranchAddress("BBgem_chi2ndf",&BBgem_chi2ndf);

  //setup global cut formula
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  //Accounting event variables
  long nevent = 0, nentries = C->GetEntries();

  //ttree formula variables
  int treenum = 0, currenttreenum = 0;

  //event loop
  while(C->GetEntry(nevent++)){

  	//single loop global cut
        currenttreenum = C->GetTreeNumber();
        if( nevent == 1 || currenttreenum != treenum ){
        treenum = currenttreenum;
        GlobalCut->UpdateFormulaLeaves();
        }
        //Is true if failed global cut
        bool failglobal = cuts::failedGlobal(GlobalCut);
	//Establish HCal Cut
	bool passHCalE = cuts::passHCalE(ehcal,ehcal_min);
	bool hcalaa_ON = cuts::hcalaa_ON(xhcal,yhcal,hcalaa);
	bool HCal_spot = cuts::passHCal_Spot(dx,dy,dxO_p,dyO_p,dxsig_p,dysig_p,spot_sig);
	bool pass_HCalCut = passHCalE && hcalaa_ON && HCal_spot;

	//Fiducial region, we want to separate the direction so not using the standard function calls
	bool fidx_cut = (xexp - proton_deflection) >= fidx_min && (xexp - proton_deflection) <= fidx_max;
	bool fidy_cut = yexp >= fidy_min && yexp <= fidy_max;

	//W2 cut
	bool goodW2 = cuts::goodW2(W2,W2_low,W2_high);

	//Make sure everything passes the global cut. The global cut only effects the e-arm
	//not fail = pass
	if(!failglobal){
		 //Check the W2 plot
		//Make sure we pass the fiducial cut
		if(fidx_cut && fidy_cut){
		h_W2_globcut->Fill(W2,final_mc_weight);
			//Check the HCal cut
			if(pass_HCalCut){
			h_W2_hcalcut->Fill(W2,final_mc_weight);
			}else{
			h_W2_anticut->Fill(W2,final_mc_weight);
			}
		}//end fiducial conditional for W2

		//Now enforce everything after this point has a W2 cut
		if(goodW2){
		
			//enforce fiducial cut in the y direction. Going to fill x direction plots
			if(fidy_cut){
			//Fill histograms for the denominator
			hx_expect_all->Fill(xexp - proton_deflection,final_mc_weight);
			hx_HCal_all->Fill(xhcal,final_mc_weight);
			hrowHCal_all->Fill(rowblkHCAL,final_mc_weight);

				//enforce the HCal cut now
				//Fill histograms for numerator
				if(pass_HCalCut){
				hx_expect_hcalcut->Fill(xexp - proton_deflection,final_mc_weight);
                        	hx_HCal_cut->Fill(xhcal,final_mc_weight);
                        	hrowHCal_cut->Fill(rowblkHCAL,final_mc_weight);
				}else{
				hx_expect_anticut->Fill(xexp - proton_deflection,final_mc_weight);
                                hx_HCal_anticut->Fill(xhcal,final_mc_weight);
                                hrowHCal_anticut->Fill(rowblkHCAL,final_mc_weight);
				}	
			}//end fidy cut


			//enforce fiducial cut in the x direction. Going to fill y direction plots
                        if(fidx_cut){
                        //Fill histograms for the denominator
                        hy_expect_all->Fill(yexp,final_mc_weight);
                        hy_HCal_all->Fill(yhcal,final_mc_weight);
                        hcolHCal_all->Fill(colblkHCAL,final_mc_weight);

                                //enforce the HCal cut now
                                //Fill histograms for numerator
                                if(pass_HCalCut){
                                hy_expect_hcalcut->Fill(yexp,final_mc_weight);
                                hy_HCal_cut->Fill(yhcal,final_mc_weight);
                                hcolHCal_cut->Fill(colblkHCAL,final_mc_weight);
                                }else{
                                hy_expect_anticut->Fill(yexp,final_mc_weight);
                                hy_HCal_anticut->Fill(yhcal,final_mc_weight);
                                hcolHCal_anticut->Fill(colblkHCAL,final_mc_weight);
                                }
                        }//end fidx cut


			//check the xy expect distribution under HCal cut but not fid
			hxy_expect_all->Fill(yexp,xexp - proton_deflection,final_mc_weight);
			if(pass_HCalCut){
			hxy_expect_hcalcut->Fill(yexp,xexp - proton_deflection,final_mc_weight);
			}else{
			hxy_expect_anticut->Fill(yexp,xexp - proton_deflection,final_mc_weight);
			}

			//For pretty much everything else to fill enforce the fiducial cut and then check which ones are passing HCal cut
			if(fidx_cut && fidy_cut){
				hdx_all->Fill(dx,final_mc_weight);
				hdy_all->Fill(dy,final_mc_weight);
				hdxdy_all->Fill(dy,dx,final_mc_weight);
				hxy_HCal_all->Fill(yhcal,xhcal,final_mc_weight);
				hehcal_all->Fill(ehcal,final_mc_weight);
				hrowcolHCal_all->Fill(colblkHCAL,rowblkHCAL,final_mc_weight);


				if(pass_HCalCut){
					hdx_cut->Fill(dx,final_mc_weight);
                               		hdy_cut->Fill(dy,final_mc_weight);
                                	hdxdy_cut->Fill(dy,dx,final_mc_weight);
                                	hxy_HCal_cut->Fill(yhcal,xhcal,final_mc_weight);
                                	hehcal_cut->Fill(ehcal,final_mc_weight);
					hrowcolHCal_cut->Fill(colblkHCAL,rowblkHCAL,final_mc_weight);
				}else{
					hdx_anticut->Fill(dx,final_mc_weight);
                                        hdy_anticut->Fill(dy,final_mc_weight);
                                        hdxdy_anticut->Fill(dy,dx,final_mc_weight);
                                        hxy_HCal_anticut->Fill(yhcal,xhcal,final_mc_weight);
                                        hehcal_anticut->Fill(ehcal,final_mc_weight);
					hrowcolHCal_anticut->Fill(colblkHCAL,rowblkHCAL,final_mc_weight);
				}
			}//end total fiducial conditional

		}//end W2 conditional

	}// end global conditional

  }//end of event loop

  //Make histogram for position dependent efficiency in the x-direction
  TH1D *heff_vs_xexpect = new TH1D( *hx_expect_hcalcut );
  heff_vs_xexpect->SetName( "heff_vs_xexpect" );
  heff_vs_xexpect->SetTitle( "heff_vs_xexpect" );
  heff_vs_xexpect->Divide( hx_expect_hcalcut, hx_expect_all );
  //calculate the corresponding bin by bin error for the x-direction efficiency
  for( int i=1; i<=heff_vs_xexpect->GetNbinsX(); i++ ){
    double eff = heff_vs_xexpect->GetBinContent(i);
    double N = std::max(1.0,hx_expect_all->GetBinContent(i));
    double err = sqrt(eff*(1.0-eff)/N); //binomial error
    heff_vs_xexpect->SetBinError(i,err);
  }

  //Make histogram for position dependent efficiency in the y-direction
  TH1D *heff_vs_yexpect = new TH1D( *hy_expect_hcalcut );
  heff_vs_yexpect->SetName( "heff_vs_yexpect" );
  heff_vs_yexpect->SetTitle( "heff_vs_yexpect" ); 
  heff_vs_yexpect->Divide( hy_expect_hcalcut, hy_expect_all );
  //calculate the corresponding bin by bin error for the y-direction efficiency
  for( int j=1; j<=heff_vs_yexpect->GetNbinsX(); j++ ){
    double eff = heff_vs_yexpect->GetBinContent(j);
    double N = std::max(1.0,hy_expect_all->GetBinContent(j));
    double err = sqrt(eff*(1.0-eff)/N); //binomial error
    heff_vs_yexpect->SetBinError(j,err);
  }

  TH2D *heff_vs_xyexpect = new TH2D( *hxy_expect_hcalcut );
  heff_vs_xyexpect->SetName("heff_vs_xyexpect" );
  heff_vs_xyexpect->SetTitle("heff_vs_xyexpect" );
  heff_vs_xyexpect->Divide( hxy_expect_hcalcut, hxy_expect_all );

  for( int i=1; i<=heff_vs_xyexpect->GetNbinsX(); i++ ){
    for( int j=1; j<=heff_vs_xyexpect->GetNbinsY(); j++ ){

      int bin = heff_vs_xyexpect->GetBin(i,j);
      
      double eff = heff_vs_xyexpect->GetBinContent(bin);
      double N = std::max(1.0,hxy_expect_all->GetBinContent(bin));
      double err = sqrt(eff*(1.0-eff)/N);
      heff_vs_xyexpect->SetBinError(bin,err);
    }
  }

  TH1D *heff_vs_x = new TH1D( *hx_HCal_cut );
  heff_vs_x->SetName("heff_vs_x");
  heff_vs_x->Divide( hx_HCal_cut, hx_HCal_all );

  TH1D *heff_vs_y = new TH1D( *hy_HCal_cut );
  heff_vs_y->SetName("heff_vs_y");
  heff_vs_y->Divide( hy_HCal_cut, hy_HCal_all );

  TH2D *heff_vs_xy = new TH2D( *hxy_HCal_cut );
  heff_vs_xy->SetName("heff_vs_xy" );
  heff_vs_xy->Divide( hxy_HCal_cut, hxy_HCal_all );

  TH1D *heff_vs_row = new TH1D( *hrowHCal_cut );
  heff_vs_row->SetName("heff_vs_row");
  heff_vs_row->Divide( hrowHCal_cut, hrowHCal_all );

  TH1D *heff_vs_col = new TH1D( *hcolHCal_cut );
  heff_vs_col->SetName("heff_vs_col");
  heff_vs_col->Divide( hcolHCal_cut, hcolHCal_all );

  TH2D *heff_vs_rowcol = new TH2D( *hrowcolHCal_cut );
  heff_vs_rowcol->SetName( "heff_vs_rowcol" );
  heff_vs_rowcol->Divide( hrowcolHCal_cut, hrowcolHCal_all );


  TH1D *heff_vs_W2 = new TH1D( *h_W2_hcalcut );
  heff_vs_W2->SetName("heff_vs_W2");
  heff_vs_W2->Divide( h_W2_hcalcut, h_W2_globcut );

  for( int i=1; i<=heff_vs_W2->GetNbinsX(); i++ ){
    double eff = heff_vs_W2->GetBinContent(i);
    double N = std::max(1.0,h_W2_globcut->GetBinContent(i));
    double err = sqrt(eff*(1.0-eff)/N);
    heff_vs_W2->SetBinError(i,err);
  }
  
  //Make the canvases which hold the fitted histogram
  TCanvas* c0 = plots::plotHCalEff(heff_vs_xexpect,"c0","heff_vs_xexpect_clone","heff_vs_xexpect_wfit",fitx_low,fitx_high);

  TCanvas* c2 = plots::plotHCalEff(heff_vs_yexpect,"c2","heff_vs_yexpect_clone","heff_vs_yexpect_wfit",fity_low,fity_high);
  
  TCanvas* c3 = plots::plot_Comp(hx_expect_all,hx_expect_hcalcut,"c3","hx_expect_all_clone","hx_expect_cut_clone");
  
  TCanvas* c4 = plots::plot_Comp(hy_expect_all,hy_expect_hcalcut,"c4","hy_expect_all_clone","hy_expect_cut_clone");

  TCanvas* c5 = plots::plot_Comp(h_W2_globcut,h_W2_hcalcut,"c5","hW2_globcut_clone","hW2_hcalcut_clone");

  //Write the info to file
  fout->Write();

  //Write stuff to a pdf
  TString plotname = outfile;
  plotname.ReplaceAll(".root",".pdf");
  TString start = Form("%s%s",plotname.Data(),"(");
  //middle is the same as the name
  TString end = Form("%s%s",plotname.Data(),")");

  c0->Print(start.Data(),"pdf");
  c2->Print(plotname.Data(),"pdf");
  c3->Print(plotname.Data(),"pdf");
  c4->Print(end.Data(),"pdf");

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;
}//end main
