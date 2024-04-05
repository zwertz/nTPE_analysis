#ifndef PLOTS_H
#define PLOTS_H

//Author Ezekiel Wertz
//A location to hold functions that make canvases or plots, beyond just standard output to root tree


namespace plots{

//A function to setup the TLines we will need to make canvases to check the HCal acceptance area and fiducial cuts
vector<TLine*> setupLines(vector<double> hcal_info, Width_t wide,Color_t datCol);

//horizontal line at X from Yi to Yf
TLine* setupLine_Horz(double Yi, double Yf, double X, Width_t wide,Color_t datCol,Style_t style);

//make the canvas to put acceptance region check on
TCanvas* plotAcceptance_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TH2D *hxy_globcut,TH2D *hxy_acceptancecut);

//make the canvas to put fiducial region check on
TCanvas* plotFid_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TLine *LineFidPro,TH2D *hxy_expect_glob_W2_cut,TH2D *hxy_expect_cut,TH2D *hxy_expect_failedfid);

//make the canvas to put fiducial neutron and proton hypothesis check on
TCanvas* plotFid_Hypothesis_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TLine *LineFidPro,TH2D *hxy_expect_fidcutn,TH2D *hxy_expect_fidcutp);

}//end namespace
#endif
