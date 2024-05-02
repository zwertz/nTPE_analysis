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

// Utility function to shift every bin of a TH1D along the x-axis
TH1D* shiftHistogramX(TH1D* origHist, double shiftValue);

//Makes residual histogram from the two given histograms also shows error.
TH1D* makeResidualWithError(TString name,TH1D* hist1, TH1D* hist2,bool match_bins, bool match_x_range);

//Makes residual histogram from a given histogram and fit, also shows error.
TH1D* makeResidualWithError(TString name,TH1D* hist, TH1* fit);

//Make a histogram with a subtracted backgroud based on a fit.
TH1D* subtractBG(TH1D* hist, TF1* bgFit);

//Make a canvas that displays the Data and MC dx plot being compared along with background from a 4th order polynomial. Display relevant yield and ratio information
TCanvas* plotDataMCFitsResiduals(TH1D* hdx_data, TH1D* hdx_mc_p, TH1D* hdx_mc_n, TF1* bg,const char *name, const char *fitName, const char* fitType, const vector<pair<double,double>> params, pair<double,double> qual,double hcalfit_low, double hcalfit_high,bool shiftfit);

//Make a canvas that displays the Data and MC dx plot being compared. With background already subtracked. Display relevant yield and ratio information
TCanvas* plotDataMCFitsResiduals_NoBG(TH1D* hdx_data, TH1D* hdx_mc_p, TH1D* hdx_mc_n,const char *name,const char *fitName, const char* fitType, const vector<pair<double,double>> params, pair<double,double> qual,double hcalfit_low, double hcalfit_high,bool shiftfit);

}//end namespace
#endif
