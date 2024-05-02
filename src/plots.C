#include "../include/plots.h"

//Author Ezekiel Wertz
//Implementations for plots related functions


namespace plots{

//A function to setup the TLines we will need to make canvases to check the HCal acceptance area and fiducial cuts
vector<TLine*> setupLines(vector<double> hcal_info, Width_t wide,Color_t datCol){

  vector<TLine*> demLines;

 //Define boundary values based on provided HCal info
 double Xi = hcal_info[0];
 double Xf = hcal_info[1];
 double Yi = hcal_info[2];
 double Yf = hcal_info[3];

 //Define lines for boundaries
 TLine *LineXi = new TLine(Yi,Xi,Yf,Xi); //horizontal line at Xi from Yi to Yf
 LineXi->SetLineWidth(wide);
 LineXi->SetLineColor(datCol);

 TLine *LineXf = new TLine(Yi,Xf,Yf,Xf); //horizontal line at Xf from Yi to Yf
 LineXf->SetLineWidth(wide);
 LineXf->SetLineColor(datCol);

 TLine *LineYi = new TLine(Yi,Xi,Yi,Xf); //vertical line at Yi from Xi to Xf
 LineYi->SetLineWidth(wide);
 LineYi->SetLineColor(datCol);
 
 TLine *LineYf = new TLine(Yf,Xi,Yf,Xf); //vertical line at Yf from Xi to Xf
 LineYf->SetLineWidth(wide);
 LineYf->SetLineColor(datCol);

 demLines.push_back(LineXi);
 demLines.push_back(LineXf);
 demLines.push_back(LineYi);
 demLines.push_back(LineYf); 

 return demLines;
}

//horizontal line at X from Yi to Yf
TLine* setupLine_Horz(double Yi, double Yf, double X, Width_t wide,Color_t datCol,Style_t style){
TLine *myLine = new TLine(Yi,X,Yf,X); 
myLine->SetLineWidth(wide);
myLine->SetLineColor(datCol);
myLine->SetLineStyle(style);

return myLine;
}

//make the canvas to put acceptance region check on
TCanvas* plotAcceptance_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TH2D *hxy_globcut,TH2D* hxy_acceptancecut){

TCanvas* myCan = new TCanvas(name,name,1600,1200); 
myCan->Divide(2,1);

//first plot
myCan->cd(1);
hxy_globcut->Draw("colz");
Lines_pos[0]->Draw("same");
Lines_pos[1]->Draw("same");
Lines_pos[2]->Draw("same");
Lines_pos[3]->Draw("same");

Lines_aa[0]->Draw("same");
Lines_aa[1]->Draw("same");
Lines_aa[2]->Draw("same");
Lines_aa[3]->Draw("same");

Lines_Fid[0]->Draw("same");
Lines_Fid[1]->Draw("same");
Lines_Fid[2]->Draw("same");
Lines_Fid[3]->Draw("same");

//make legend
auto legend1 = new TLegend(0.17,0.8,0.52,0.9);
legend1->SetTextSize(0.03);
legend1->AddEntry(Lines_pos[0],"HCal Boundary","l");
legend1->AddEntry(Lines_aa[0],"Acceptance Region","l");
legend1->AddEntry(Lines_Fid[0],"Fiducial Region","l");
legend1->Draw();

//second plot
myCan->cd(2);
hxy_acceptancecut->Draw("colz");
Lines_pos[0]->Draw("same");
Lines_pos[1]->Draw("same");
Lines_pos[2]->Draw("same");
Lines_pos[3]->Draw("same");

Lines_aa[0]->Draw("same");
Lines_aa[1]->Draw("same");
Lines_aa[2]->Draw("same");
Lines_aa[3]->Draw("same");

Lines_Fid[0]->Draw("same");
Lines_Fid[1]->Draw("same");
Lines_Fid[2]->Draw("same");
Lines_Fid[3]->Draw("same");

//make legend
auto legend2 = new TLegend(0.17,0.8,0.52,0.9);
legend2->SetTextSize(0.03);
legend2->AddEntry(Lines_pos[0],"HCal Boundary","l");
legend2->AddEntry(Lines_aa[0],"Acceptance Region","l");
legend2->AddEntry(Lines_Fid[0],"Fiducial Region","l");
legend2->Draw();

return myCan;
}

//make the canvas to put fiducial region check on
TCanvas* plotFid_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TLine *LineFidPro,TH2D *hxy_expect_glob_W2_cut,TH2D *hxy_expect_cut,TH2D *hxy_expect_failedfid){

TCanvas* myCan = new TCanvas(name,name,1600,1200);
myCan->Divide(3,1);

//first plot
myCan->cd(1);
hxy_expect_glob_W2_cut->Draw("colz");
Lines_pos[0]->Draw("same");
Lines_pos[1]->Draw("same");
Lines_pos[2]->Draw("same");
Lines_pos[3]->Draw("same");

Lines_aa[0]->Draw("same");
Lines_aa[1]->Draw("same");
Lines_aa[2]->Draw("same");
Lines_aa[3]->Draw("same");

Lines_Fid[0]->Draw("same");
Lines_Fid[1]->Draw("same");
Lines_Fid[2]->Draw("same");
Lines_Fid[3]->Draw("same");

//make legend
auto legend1 = new TLegend(0.17,0.8,0.52,0.9);
legend1->SetTextSize(0.03);
legend1->AddEntry(Lines_pos[0],"HCal Boundary","l");
legend1->AddEntry(Lines_aa[0],"Acceptance Region","l");
legend1->AddEntry(Lines_Fid[0],"Fiducial Region","l");
legend1->Draw();


//second plot
myCan->cd(2);
hxy_expect_cut->Draw("colz");
Lines_pos[0]->Draw("same");
Lines_pos[1]->Draw("same");
Lines_pos[2]->Draw("same");
Lines_pos[3]->Draw("same");

Lines_aa[0]->Draw("same");
Lines_aa[1]->Draw("same");
Lines_aa[2]->Draw("same");
Lines_aa[3]->Draw("same");

Lines_Fid[0]->Draw("same");
Lines_Fid[1]->Draw("same");
Lines_Fid[2]->Draw("same");
Lines_Fid[3]->Draw("same");

LineFidPro->Draw("same");

//make legend
auto legend2 = new TLegend(0.17,0.8,0.52,0.9);
legend2->SetTextSize(0.03);
legend2->AddEntry(Lines_pos[0],"HCal Boundary","l");
legend2->AddEntry(Lines_aa[0],"Acceptance Region","l");
legend2->AddEntry(Lines_Fid[0],"Fiducial Region","l");
legend2->Draw();

//third plot
myCan->cd(3);
hxy_expect_failedfid->Draw("colz");

Lines_pos[0]->Draw("same");
Lines_pos[1]->Draw("same");
Lines_pos[2]->Draw("same");
Lines_pos[3]->Draw("same");

Lines_aa[0]->Draw("same");
Lines_aa[1]->Draw("same");
Lines_aa[2]->Draw("same");
Lines_aa[3]->Draw("same");

Lines_Fid[0]->Draw("same");
Lines_Fid[1]->Draw("same");
Lines_Fid[2]->Draw("same");
Lines_Fid[3]->Draw("same");

LineFidPro->Draw("same");

//make legend
auto legend3 = new TLegend(0.17,0.8,0.52,0.9);
legend3->SetTextSize(0.03);
legend3->AddEntry(Lines_pos[0],"HCal Boundary","l");
legend3->AddEntry(Lines_aa[0],"Acceptance Region","l");
legend3->AddEntry(Lines_Fid[0],"Fiducial Region","l");
legend3->Draw();


return myCan;
}

//make the canvas to put fiducial neutron and proton hypothesis check on
TCanvas* plotFid_Hypothesis_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TLine *LineFidPro,TH2D *hxy_expect_fidcutn,TH2D *hxy_expect_fidcutp){

TCanvas* myCan = new TCanvas(name,name,1600,1200);
myCan->Divide(2,1);

//first plot
myCan->cd(1);
hxy_expect_fidcutn->Draw("colz");

Lines_pos[0]->Draw("same");
Lines_pos[1]->Draw("same");
Lines_pos[2]->Draw("same");
Lines_pos[3]->Draw("same");

Lines_aa[0]->Draw("same");
Lines_aa[1]->Draw("same");
Lines_aa[2]->Draw("same");
Lines_aa[3]->Draw("same");

Lines_Fid[0]->Draw("same");
Lines_Fid[1]->Draw("same");
Lines_Fid[2]->Draw("same");
Lines_Fid[3]->Draw("same");

LineFidPro->Draw("same");

//make legend
auto legend1 = new TLegend(0.17,0.8,0.52,0.9);
legend1->SetTextSize(0.03);
legend1->AddEntry(Lines_pos[0],"HCal Boundary","l");
legend1->AddEntry(Lines_aa[0],"Acceptance Region","l");
legend1->AddEntry(Lines_Fid[0],"Fiducial Region","l");
legend1->Draw();


//second plot
myCan->cd(2);
hxy_expect_fidcutp->Draw("colz");
Lines_pos[0]->Draw("same");
Lines_pos[1]->Draw("same");
Lines_pos[2]->Draw("same");
Lines_pos[3]->Draw("same");

Lines_aa[0]->Draw("same");
Lines_aa[1]->Draw("same");
Lines_aa[2]->Draw("same");
Lines_aa[3]->Draw("same");

Lines_Fid[0]->Draw("same");
Lines_Fid[1]->Draw("same");
Lines_Fid[2]->Draw("same");
Lines_Fid[3]->Draw("same");

//make legend
auto legend2 = new TLegend(0.17,0.8,0.52,0.9);
legend2->SetTextSize(0.03);
legend2->AddEntry(Lines_pos[0],"HCal Boundary","l");
legend2->AddEntry(Lines_aa[0],"Acceptance Region","l");
legend2->AddEntry(Lines_Fid[0],"Fiducial Region","l");
legend2->Draw();


return myCan;
}

// Utility function to shift every bin of a TH1D along the x-axis
TH1D* shiftHistogramX(TH1D* origHist, double shiftValue){

//total number of entries in the histogram
double totEntries = origHist->GetEntries();

//create a histogram that is indentical to the original
TH1D* shiftHist = (TH1D*) origHist->Clone("shiftHist");

//clear the contents of the cloned histogram
shiftHist->Reset();

	//apply the shift to each bin, bins start at 1
	for(int i=1; i<=origHist->GetNbinsX(); ++i){
	//get information about the bin in the orig hist
	double origBinCenter = origHist->GetBinCenter(i);
	double origBinContent = origHist->GetBinContent(i);
	double origBinError = origHist->GetBinError(i);
	
	//find the bin in the shifted histogram that corresponds to the new bin center
	int newBinCenter = shiftHist->FindBin(origBinCenter + shiftValue);

	//add the content and the error to the new bin
	//If multiple old bins shift into the same new bin, their contents and errors are added
	double newBinContent = shiftHist->GetBinContent(newBinCenter) + origBinContent;
	double newBinError = sqrt(pow(shiftHist->GetBinError(newBinCenter),2) + pow(origBinError,2) );

	//update the shifted histogram with the content and error info
	shiftHist->SetBinContent(newBinCenter,newBinContent);
	shiftHist->SetBinError(newBinCenter,newBinError);

	}

//set the total number of entries
shiftHist->SetEntries(totEntries);

return shiftHist;
}

//Makes residual histogram from the two given histograms also shows error.
TH1D* makeResidualWithError(TString name,TH1D* hist1, TH1D* hist2,bool match_bins, bool match_x_range){

//Get some info about histo1
TString hist1_name = hist1->GetName();
int hist1_Nbins = hist1->GetNbinsX();
double hist1_minX = hist1->GetXaxis()->GetXmin();
double hist1_maxX = hist1->GetXaxis()->GetXmax();

//Get some info about histo2
TString hist2_name = hist2->GetName();
int hist2_Nbins = hist2->GetNbinsX();
double hist2_minX = hist2->GetXaxis()->GetXmin();
double hist2_maxX = hist2->GetXaxis()->GetXmax();

//check histograms for binning and range if required
	if(match_bins && (hist1_Nbins != hist2_Nbins)){
	cout << "Histograms need to have matching bin numbers for accurate residual calculation." << endl;
	return nullptr;
	}
	if(match_x_range && ((hist1_minX != hist2_minX)||(hist1_maxX != hist2_maxX))){
	cout << "Histograms need to have matching x ranges for accurate residual calcuation." << endl;
	}
//create the histogram for the residual
TH1D* resid_hist = new TH1D(Form("Residual Histogram with Errors - %s: %s and %s",name.Data(),hist1_name.Data(),hist2_name.Data()),Form("Residual Histogram with Errors - %s: %s and %s",name.Data(),hist1_name.Data(),hist2_name.Data()),hist1_Nbins,hist1_minX,hist1_maxX);

	//go bin by bin starting from 1 and populate the residual histogram
	for(int bin = 1;bin<=hist1_Nbins; bin++ ){
	//hist1 info
	double value1 = hist1->GetBinContent(bin);
	double error1 = hist1->GetBinError(bin);

	//hist2 info
	double value2 = hist2->GetBinContent(bin);
        double error2 = hist2->GetBinError(bin);

	//calculate residual
	double residual = value1-value2;
	double res_error = sqrt(pow(error1,2) + pow(error2,2));

	resid_hist->SetBinContent(bin,residual);
	resid_hist->SetBinError(bin,res_error);

	}

return resid_hist;
}

//Makes residual histogram from a given histogram and fit, also shows error.
TH1D* makeResidualWithError(TString name,TH1D* hist, TF1* fit){
double totentries = hist->GetEntries();

int hist_Nbins = hist->GetNbinsX();
double hist_minX = hist->GetXaxis()->GetXmin();
double hist_maxX = hist->GetXaxis()->GetXmax();
TString hist_name = hist->GetName();
TString fit_name = fit->GetName();

//create the histogram for the residual
TH1D* resid_hist = new TH1D(Form("Residual Histogram with Errors - %s: %s and %s",name.Data(),hist_name.Data(),fit_name.Data()),Form("Residual Histogram with Errors - %s: %s and %s",name.Data(),hist_name.Data(),fit_name.Data()),hist_Nbins,hist_minX,hist_maxX);

	//go bin by bin starting from 1 and populate the residual histogram
	for(int bin = 1;bin<=hist_Nbins; bin++ ){
	double hist_val = hist->GetBinContent(bin);
	double hist_error = hist->GetBinError(bin);
	double fit_val = fit->Eval(hist->GetXaxis()->GetBinCenter(bin));
	//Might need to actually add this in quadrature, currently just keeping old hist bin error

	double new_val = hist_val - fit_val;
	resid_hist->SetBinContent(bin,new_val);
	resid_hist->SetBinError(bin,hist_error);
	}	
resid_hist->SetEntries(totentries);
return resid_hist;
}

//Make a histogram with a subtracted backgroud based on a fit.
TH1D* subtractBG(TH1D* hist, TF1* bgFit){
double totentries = hist->GetEntries();

int hist_Nbins = hist->GetNbinsX();
double hist_minX = hist->GetXaxis()->GetXmin();
double hist_maxX = hist->GetXaxis()->GetXmax();

TH1D* hist_sub = new TH1D("hist_sub","",hist_Nbins,hist_minX,hist_maxX);

	//go by each bin and subtract off the backgroun value from the hist and populate the new hist
	for(int bin =1; bin <= hist_Nbins; ++bin){
	double bgValue = bgFit->Eval(hist->GetXaxis()->GetBinCenter(bin));
	double newValue = hist->GetBinContent(bin) - bgValue;
	//Might need to actually add this in quadrature, currently just keeping old hist bin error
	double histError = hist->GetBinError(bin);
	

	hist_sub->SetBinContent(bin, newValue);
	hist_sub->SetBinError(bin,histError);
	}

hist_sub->SetEntries(totentries);
return hist_sub;
}

//Make a canvas that displays the Data and MC dx plot being compared along with background from a 4th order polynomial. Display relevant yield and ratio information
TCanvas* plotDataMCFitsResiduals(TH1D* hdx_data, TH1D* hdx_mc_p, TH1D* hdx_mc_n, TF1* bg,const char *name,const char *fitName, const char* fitType, const vector<pair<double,double>> params, pair<double,double> qual,double hcalfit_low, double hcalfit_high,bool shiftfit){

//Recreate the fit
TF1* fit = new TF1("fit",fitType,hcalfit_low,hcalfit_high,params.size());
	//file the fit with the parameters
	for(int i=0; i < params.size();++i){
	fit->SetParameter(i,params[i].first);
	fit->SetParError(i,params[i].second);
	}
fit->SetChisquare(qual.first);
fit->SetNDF(qual.second);

//Create a canvas we will use to store information
TCanvas *myTotCan = new TCanvas(name,Form("Data with Fits and Residuals %s",fitName),1600,1200);
//create two pads on the canvas
TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);

//Draw the pads and set the margins
pad1->SetBottomMargin(0.1); // Upper and lower plot are not joined
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);
pad1->Draw();
pad2->Draw();

//Draw the histogram with fits on the first pad
pad1->cd();
//data dx
hdx_data->SetLineColor(kBlack);
hdx_data->SetLineWidth(1);
hdx_data->SetTitle(Form("dx, %s;m",fitName));
hdx_data->SetStats(0);
hdx_data->Draw("hist");
//fit the data dx
fit->SetLineColor(kMagenta);
fit->SetLineWidth(2);
fit->Draw("same");
//proton MC dx
	if(shiftfit){
	hdx_mc_p = plots::shiftHistogramX(hdx_mc_p,params[2].first);
	}
hdx_mc_p->Scale(params[0].first);
hdx_mc_p->SetLineColor(kRed);
hdx_mc_p->SetFillColorAlpha(kRed-9,0.5);
hdx_mc_p->SetFillStyle(3001);
hdx_mc_p->SetLineWidth(2);
hdx_mc_p->Draw("hist same");
//neutron MC dx
	if(shiftfit){
	hdx_mc_n = plots::shiftHistogramX(hdx_mc_n,params[3].first);
	}
hdx_mc_n->Scale(params[1].first);
hdx_mc_n->SetLineColor(kBlue);
hdx_mc_n->SetFillColorAlpha(kBlue-9,0.5);
hdx_mc_n->SetFillStyle(3001);
hdx_mc_n->SetLineWidth(2);
hdx_mc_n->Draw("hist same");

//background function
bg->SetLineColor(kGreen);
bg->SetFillColorAlpha(kGreen-9,0.5);
bg->SetFillStyle(1001);
bg->Draw("same");

//yield info for protons
double p_sum_error;
double p_sum = hdx_mc_p->IntegralAndError(0,hdx_mc_p->GetNbinsX()+1,p_sum_error,"");

double p_scale = params[0].first;
double p_scale_error = params[0].second;

//yield info for neutrons
double n_sum_error;
double n_sum = hdx_mc_n->IntegralAndError(0,hdx_mc_n->GetNbinsX()+1,n_sum_error,"");

double n_scale = params[1].first;
double n_scale_error = params[1].second;

//get the ratios and the errors
double np_sum_ratio = n_sum/p_sum;
double np_scale_ratio = n_scale/p_scale;

double np_sum_ratio_error = np_sum_ratio * sqrt(pow(n_sum_error/n_sum,2) + pow(p_sum_error/p_sum,2));
double np_scale_ratio_error = np_scale_ratio * sqrt(pow(n_scale_error/n_scale,2) + pow(p_scale_error/p_scale,2));

//Make a legend to hold all the information we care about
TLegend* legend = new TLegend(0.68, 0.5, 0.9, 0.9);
legend->AddEntry(hdx_data,"Data","l");
legend->AddEntry(fit, "Total Fit","l");
legend->AddEntry(bg,"Background","l");
legend->AddEntry(hdx_mc_p,"Proton SIMC MC","l");
legend->AddEntry(hdx_mc_n,"Neutron SIMC MC","l");
legend->AddEntry((TObject*)0,Form("data N events : %0.0f",hdx_data->GetEntries()),"");
legend->AddEntry((TObject*)0,Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_scale_ratio,np_scale_ratio_error),"");
legend->AddEntry((TObject*)0,Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio,np_sum_ratio_error),"");
legend->AddEntry((TObject*)0,Form("dx shift pars, n/p : %0.3f / %0.3f ",params[3].first,params[2].first),"");
legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()),"");
legend->Draw("same");

//Residual plot
//Add the two MC histograms together for direct comparison
TH1D* sumHist = new TH1D("sumHist1","dx for both MC protons and neutrons",hdx_mc_p->GetNbinsX(),hdx_mc_p->GetXaxis()->GetXmin(),hdx_mc_p->GetXaxis()->GetXmax());

sumHist->Add(hdx_mc_p,hdx_mc_n);

//background subtraction to directly compare
	for(int bin=1; bin<= hdx_data->GetNbinsX(); ++bin){
	double bg_value = bg->Eval(hdx_data->GetXaxis()->GetBinCenter(bin));
	hdx_data->SetBinContent(bin,hdx_data->GetBinContent(bin) - bg_value); 
	hdx_data->SetBinError(bin,hdx_data->GetBinError(bin));
	}
//Make the residual plot
TH1D* residual = plots::makeResidualWithError("dx",hdx_data,sumHist,true, false);
TH1D* residual_fit = plots::makeResidualWithError("dx_fit",hdx_data,fit);

residual->SetTitle("");
residual->GetXaxis()->SetRangeUser(hcalfit_low,hcalfit_high);
residual->SetLineColor(kBlue+2);
residual->SetLineWidth(2);
residual->SetStats(0);
residual->GetYaxis()->SetTitle("Residuals");
residual->GetYaxis()->SetTitleSize(0.07);
residual->GetYaxis()->SetTitleOffset(0.5);
residual->GetYaxis()->SetLabelSize(0.07);
residual_fit->SetLineColor(kOrange-3);
residual_fit->SetLineWidth(2);
residual_fit->SetStats(0);

//Make a legend to hold all the information we care about
TLegend* legend_res = new TLegend(0.8, 0.75, 0.9, 1.0);
legend_res->AddEntry(residual,"Data-MC combine","l");
legend_res->AddEntry(residual_fit,"Data-Total Fit","l");

//Draw the residual on the second pad
pad2->cd();
residual->Draw("hist");
residual_fit->Draw("hist same");
legend_res->Draw("same");
//Draw a line at y=0 just for better reading of the residual
TF1* zeroLine = new TF1("zeroLine","0",residual->GetXaxis()->GetXmin(),residual->GetXaxis()->GetXmax());
zeroLine->SetLineColor(kBlack);
zeroLine->SetLineStyle(2);
zeroLine->Draw("same");

return myTotCan;
}


//Make a canvas that displays the Data and MC dx plot being compared. With background already subtracked. Display relevant yield and ratio information
TCanvas* plotDataMCFitsResiduals_NoBG(TH1D* hdx_data, TH1D* hdx_mc_p, TH1D* hdx_mc_n,const char *name,const char *fitName, const char* fitType, const vector<pair<double,double>> params, pair<double,double> qual,double hcalfit_low, double hcalfit_high,bool shiftfit){
//Recreate the fit
TF1* fit = new TF1("fit",fitType,hcalfit_low,hcalfit_high,params.size());
        //file the fit with the parameters
	for(int i=0; i < params.size();++i){
        fit->SetParameter(i,params[i].first);
        fit->SetParError(i,params[i].second);
        }
fit->SetChisquare(qual.first);
fit->SetNDF(qual.second);

//Create a canvas we will use to store information
TCanvas *myTotCan = new TCanvas(name,Form("Data with Fits and Residuals %s",fitName),1600,1200);
//create two pads on the canvas
TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);

//Draw the pads and set the margins
pad1->SetBottomMargin(0.1); // Upper and lower plot are not joined
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);
pad1->Draw();
pad2->Draw();

//Draw the histogram with fits on the first pad
pad1->cd();
//data dx
hdx_data->SetLineColor(kBlack);
hdx_data->SetLineWidth(1);
hdx_data->SetTitle(Form("dx, %s;m",fitName));
hdx_data->SetStats(0);
hdx_data->Draw("hist");
//fit the data dx
fit->SetLineColor(kMagenta);
fit->SetLineWidth(2);
fit->Draw("same");
//proton MC dx
        if(shiftfit){
        hdx_mc_p = plots::shiftHistogramX(hdx_mc_p,params[2].first);
        }
hdx_mc_p->Scale(params[0].first);
hdx_mc_p->SetLineColor(kRed);
hdx_mc_p->SetFillColorAlpha(kRed-9,0.5);
hdx_mc_p->SetFillStyle(3001);
hdx_mc_p->SetLineWidth(2);
hdx_mc_p->Draw("hist same");
//neutron MC dx
        if(shiftfit){
        hdx_mc_n = plots::shiftHistogramX(hdx_mc_n,params[3].first);
        }
hdx_mc_n->Scale(params[1].first);
hdx_mc_n->SetLineColor(kBlue);
hdx_mc_n->SetFillColorAlpha(kBlue-9,0.5);
hdx_mc_n->SetFillStyle(3001);
hdx_mc_n->SetLineWidth(2);
hdx_mc_n->Draw("hist same");

//yield info for protons
double p_sum_error;
double p_sum = hdx_mc_p->IntegralAndError(0,hdx_mc_p->GetNbinsX()+1,p_sum_error,"");

double p_scale = params[0].first;
double p_scale_error = params[0].second;

//yield info for neutrons
double n_sum_error;
double n_sum = hdx_mc_n->IntegralAndError(0,hdx_mc_n->GetNbinsX()+1,n_sum_error,"");

double n_scale = params[1].first;
double n_scale_error = params[1].second;

//get the ratios and the errors
double np_sum_ratio = n_sum/p_sum;
double np_scale_ratio = n_scale/p_scale;

double np_sum_ratio_error = np_sum_ratio * sqrt(pow(n_sum_error/n_sum,2) + pow(p_sum_error/p_sum,2));
double np_scale_ratio_error = np_scale_ratio * sqrt(pow(n_scale_error/n_scale,2) + pow(p_scale_error/p_scale,2));

//Make a legend to hold all the information we care about
TLegend* legend = new TLegend(0.68, 0.5, 0.9, 0.9);
legend->AddEntry(hdx_data,"Data","l");
legend->AddEntry(fit, "Total Fit","l");
legend->AddEntry(hdx_mc_p,"Proton SIMC MC","l");
legend->AddEntry(hdx_mc_n,"Neutron SIMC MC","l");
legend->AddEntry((TObject*)0,Form("data N events : %0.0f",hdx_data->GetEntries()),"");
legend->AddEntry((TObject*)0,Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_scale_ratio,np_scale_ratio_error),"");
legend->AddEntry((TObject*)0,Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio,np_sum_ratio_error),"");
legend->AddEntry((TObject*)0,Form("dx shift pars, n/p : %0.3f / %0.3f ",params[3].first,params[2].first),"");
legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()),"");
legend->Draw("same");

//Residual plot
//Add the two MC histograms together for direct comparison
TH1D* sumHist = new TH1D("sumHist2","dx for both MC protons and neutrons",hdx_mc_p->GetNbinsX(),hdx_mc_p->GetXaxis()->GetXmin(),hdx_mc_p->GetXaxis()->GetXmax());

sumHist->Add(hdx_mc_p,hdx_mc_n);
TH1D* residual = plots::makeResidualWithError("dx",hdx_data,sumHist,true, false);
TH1D* residual_fit = plots::makeResidualWithError("dx_fit",hdx_data,fit);

residual->SetTitle("");
residual->GetXaxis()->SetRangeUser(hcalfit_low,hcalfit_high);
residual->SetLineColor(kBlue+2);
residual->SetLineWidth(2);
residual->SetStats(0);
residual->GetYaxis()->SetTitle("Residuals");
residual->GetYaxis()->SetTitleSize(0.07);
residual->GetYaxis()->SetTitleOffset(0.5);
residual->GetYaxis()->SetLabelSize(0.07);
residual_fit->SetLineColor(kOrange-3);
residual_fit->SetLineWidth(2);
residual_fit->SetStats(0);

//Make a legend to hold all the information we care about
TLegend* legend_res = new TLegend(0.8, 0.75, 0.9, 1.0);
legend_res->AddEntry(residual,"Data-MC combine","l");
legend_res->AddEntry(residual_fit,"Data-Total Fit","l");

//Draw the residual on the second pad
pad2->cd();
residual->Draw("hist");
residual_fit->Draw("hist same");
legend_res->Draw("same");

//Draw a line at y=0 just for better reading of the residual
TF1* zeroLine = new TF1("zeroLine","0",residual->GetXaxis()->GetXmin(),residual->GetXaxis()->GetXmax());
zeroLine->SetLineColor(kBlack);
zeroLine->SetLineStyle(2);
zeroLine->Draw("same");

return myTotCan;

}



}//end namespace
