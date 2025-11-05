#include "../include/plots.h"

//Author Ezekiel Wertz
//A location to hold functions that make canvases or plots, beyond just standard output to root tree. These plots tend to be repeatedly of interest throughout the analysis. Or were used more than once and therefore implemented as functions. Here each function has some basic description.


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

//vertical line at Y from Xi to Xf
TLine* setupLine_Vert(double Xi, double Xf, double Y, Width_t wide,Color_t datCol,Style_t style){
TLine *myLine = new TLine(Y,Xi,Y,Xf);
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
hxy_expect_glob_W2_cut->SetStats(0);
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
hxy_expect_cut->SetStats(0);
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
hxy_expect_failedfid->SetStats(0);
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
hxy_expect_fidcutn->SetStats(0);
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
hxy_expect_fidcutp->SetStats(0);
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
TH1D* subtractBG(TH1D* hist, TF1* bgFit, TFitResultPtr fit_ptr){
double totentries = hist->GetEntries();

int hist_Nbins = hist->GetNbinsX();
double hist_minX = hist->GetXaxis()->GetXmin();
double hist_maxX = hist->GetXaxis()->GetXmax();
TString hist_name = hist->GetName();
TString fit_name = bgFit->GetName();

TH1D* hist_sub = new TH1D(Form("hist_sub: %s and %s",hist_name.Data(),fit_name.Data()),"",hist_Nbins,hist_minX,hist_maxX);

	//go by each bin and subtract off the backgroun value from the hist and populate the new hist
	for(int bin =1; bin <= hist_Nbins; ++bin){
	double bgValue = bgFit->Eval(hist->GetXaxis()->GetBinCenter(bin));
	double newValue = hist->GetBinContent(bin) - bgValue;
	//Add the errors from the histogram and the fit error in quadrature. See fit error function for details
	double histError = hist->GetBinError(bin);
	//This only work for a poly order 4
	//double bgError = fits::FitErrorFunc(hist->GetXaxis()->GetBinCenter(bin),fit_ptr);

	//cout << "Error when subtracting BG Hist:" << histError << " Fit: " << bgError << endl;

	//double newError = sqrt(pow(histError,2) + pow(bgError,2));

	hist_sub->SetBinContent(bin, newValue);
	//hist_sub->SetBinError(bin,newError);
	hist_sub->SetBinError(bin,histError);
	}

hist_sub->SetEntries(totentries);
return hist_sub;
}

//Make a histogram, subtract hist2 from hist2
TH1D* subtractHist(TH1D* hist1,TH1D* hist2){
double totentries= hist1->GetEntries() - hist2->GetEntries();
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

TH1D* sub_hist = new TH1D(Form("Histogram with Errors, %s - %s",hist1_name.Data(),hist2_name.Data()),Form("Histogram with Errors, %s - %s",hist1_name.Data(),hist2_name.Data()),hist1_Nbins,hist1_minX,hist1_maxX);

        //go bin by bin starting from 1 and populate the residual histogram
	for(int bin = 1;bin<=hist1_Nbins; bin++ ){
        //hist1 info
	double value1 = hist1->GetBinContent(bin);
        double error1 = hist1->GetBinError(bin);

        //hist2 info
	double value2 = hist2->GetBinContent(bin);
        double error2 = hist2->GetBinError(bin);

        //calculate subtraction
	double residual = value1-value2;
        double res_error = sqrt(pow(error1,2) + pow(error2,2));

        sub_hist->SetBinContent(bin,residual);
        sub_hist->SetBinError(bin,res_error);

        }
sub_hist->SetEntries(totentries);
return sub_hist;
}

//Make a canvas that displays the Data and MC dx plot being compared along with background from a 4th order polynomial. Display relevant yield and ratio information

TCanvas* plotDataMCFitsResiduals(fit_histogram* myFitHisto,TF1* bg,const char *name, bool shiftfit, std::ofstream& report, calc_FFs_RCS_obj daFFInfo){

//Initialize some information from the fit_histogram object
TH1D* hdx_data = myFitHisto->get_hist_data();
TH1D* hdx_mc_p = myFitHisto->get_hist_p();
TH1D* hdx_mc_n = myFitHisto->get_hist_n();
string fitName = myFitHisto->get_fitName();
string fitType = myFitHisto->get_fitType();
vector<pair<double,double>> params = myFitHisto->get_fitParamsErrs();
double hcalfit_low = myFitHisto->get_xMin();
double hcalfit_high = myFitHisto->get_xMax();
double n_shift = myFitHisto->get_shift_n();
double p_shift = myFitHisto->get_shift_p();

//Recreate the fit
TF1* fit;

	//This conditional works but we will need to modify it if we need any new types of functions
	if(fitType == "fitFull_polyBG"){
	fit = new TF1("fit",myFitHisto,&fit_histogram::fitFull_polyBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFull_polyBG");
	}else if(fitType =="fitFullShift_polyBG"){
	fit = new TF1("fit",myFitHisto,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFullShift_polyBG");
	}else if(fitType =="fitFullShift_gaussBG"){
	fit = new TF1("fit",myFitHisto,&fit_histogram::fitFullShift_gaussBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFullShift_gaussBG");
	}else if(fitType =="fitFullShift_InterpolateBG"){
	fit = new TF1("fit",myFitHisto,&fit_histogram::fitFullShift_InterpolateBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFullShift_InterpolateBG");
	}else{
	cout << "The plot function you are trying to implement " << fitType << " is no good! Figure it out now!" << endl;
	fit = new TF1("fit",myFitHisto,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFullShift_polyBG");
	}
	//file the fit with the parameters
	for(int i=0; i < params.size();++i){
	fit->SetParameter(i,params[i].first);
	fit->SetParError(i,params[i].second);
	}
fit->SetChisquare(myFitHisto->get_ChiSq());
fit->SetNDF(myFitHisto->get_NDF());

//Make clones of the mc histograms
//TH1D *hdx_mc_p_before = (TH1D*)(hdx_mc_p->Clone("hdx_mc_p_before"));
TH1D *hdx_mc_p_after = (TH1D*)(hdx_mc_p->Clone("hdx_mc_p_after"));

//TH1D *hdx_mc_n_before = (TH1D*)(hdx_mc_n->Clone("hdx_mc_n_before"));
TH1D *hdx_mc_n_after = (TH1D*)(hdx_mc_n->Clone("hdx_mc_n_after"));

//Create a canvas we will use to store information
TCanvas *myTotCan = new TCanvas(name,Form("Data with Fits and Residuals %s",fitName.c_str()),1600,1200);
//create two pads on the canvas
TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);
pad1->Draw();
pad2->Draw();

//Draw the pads and set the margins
pad1->SetBottomMargin(0.1); // Upper and lower plot are not joined
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);

//Draw the histogram with fits on the first pad
pad1->cd();
//data dx
hdx_data->SetLineColor(kBlack);
hdx_data->SetLineWidth(1);
hdx_data->SetTitle(Form("dx, %s;m",fitName.c_str()));
hdx_data->SetStats(0);
hdx_data->Draw("hist");
//fit the data dx
fit->SetLineColor(kMagenta);
fit->SetLineWidth(2);
fit->Draw("same");
//proton MC dx
	if(shiftfit){
	//hdx_mc_p_before = plots::shiftHistogramX(hdx_mc_p_before,p_shift);
        hdx_mc_p_after = plots::shiftHistogramX(hdx_mc_p_after,p_shift);
	}

//Yield info and integrals come before scaling to get sigma_n/sigma_p from MC. Need to further scale.

//yield info for protons
//double p_sum_mc_error;
//double p_sum_mc = hdx_mc_p_before->IntegralAndError(0,hdx_mc_p_before->GetNbinsX()+1,p_sum_mc_error,"");

//proton scale info
double p_scale = myFitHisto->get_scale_p();
double p_scale_error = myFitHisto->get_scale_p_err();

//neutron scale info
//double n_scale = params[1].first;
//double n_scale_error = params[1].second;
//Get Rsf from the fit parameters
double R_sf = myFitHisto->get_Rsf();
double R_sf_error = myFitHisto->get_Rsf_err();

//scale info for neutrons
double n_scale = myFitHisto->get_scale_n();
double n_scale_error = myFitHisto->get_scale_n_err();

hdx_mc_p_after->Scale(p_scale);
hdx_mc_p_after->SetLineColor(kRed);
hdx_mc_p_after->SetFillColorAlpha(kRed-9,0.5);
hdx_mc_p_after->SetFillStyle(1001);
hdx_mc_p_after->SetLineWidth(2);
hdx_mc_p_after->Draw("hist same");
//neutron MC dx
	if(shiftfit){
	//hdx_mc_n_before = plots::shiftHistogramX(hdx_mc_n_before,n_shift);
        hdx_mc_n_after = plots::shiftHistogramX(hdx_mc_n_after,n_shift);
	}
//yield info for neutrons
//double n_sum_mc_error;
//double n_sum_mc = hdx_mc_n_before->IntegralAndError(0,hdx_mc_n_before->GetNbinsX()+1,n_sum_mc_error,"");
	
hdx_mc_n_after->Scale(n_scale);
hdx_mc_n_after->SetLineColor(kBlue);
hdx_mc_n_after->SetFillColorAlpha(kBlue-9,0.5);
hdx_mc_n_after->SetFillStyle(1001);
hdx_mc_n_after->SetLineWidth(2);
hdx_mc_n_after->Draw("hist same");


//background function
bg->SetLineColor(kGreen);
bg->SetFillColorAlpha(kGreen-9,0.5);
bg->SetFillStyle(1001);
bg->Draw("same");

//Just MC
//Instead calculate the reduced cross-section for protons and neutrons using the FF parameterizations in SIMC. Here we are using a class developed by me to do this purpose
double n_red_cross_section = daFFInfo.get_reduced_cross_section_n();
double n_red_cross_section_err = daFFInfo.get_reduced_cross_section_n_err();

double p_red_cross_section = daFFInfo.get_reduced_cross_section_p();
double p_red_cross_section_err = daFFInfo.get_reduced_cross_section_p_err();
//Then do the n/p ratio. This is what we call Rsim
double np_red_cross_section_ratio = daFFInfo.get_reduced_cross_section_ratio();
double np_red_cross_section_ratio_err = daFFInfo.get_reduced_cross_section_ratio_err();

//Scale the MC information by the n and p scale factor or Rsf
double n_red_cross_section_exp = n_red_cross_section*n_scale;
double n_red_cross_section_exp_err = n_red_cross_section_exp * sqrt(pow(n_scale_error/n_scale,2) + pow(n_red_cross_section_err/n_red_cross_section,2));

double p_red_cross_section_exp = p_red_cross_section*p_scale;
double p_red_cross_section_exp_err = p_red_cross_section_exp * sqrt(pow(p_scale_error/p_scale,2) + pow(p_red_cross_section_err/p_red_cross_section,2));

//Yield info and integrals come after scaling to get sigma_n/sigma_p for experiment. Don't need further scale.

//The yield method is being turned off. As it includes detection efficiency and and nuclear and radiative effects baked in
//yield info for protons
//double p_sum_exp_error;
//double p_sum_exp = hdx_mc_p_after->IntegralAndError(0,hdx_mc_p_after->GetNbinsX()+1,p_sum_exp_error,"");

//yield info for neutrons
//double n_sum_exp_error;
//double n_sum_exp = hdx_mc_n_after->IntegralAndError(0,hdx_mc_n_after->GetNbinsX()+1,n_sum_exp_error,"");

//double np_sum_mc_ratio = n_sum_mc/p_sum_mc;
//double np_sum_exp_ratio = n_sum_exp/p_sum_exp;

//get the exp np ratio
double np_red_cross_section_ratio_exp = R_sf*np_red_cross_section_ratio;
//get the np_scale ratio, we call it Rsf
//double np_scale_ratio = n_scale/p_scale;

//Get the corresponding errors
//double np_sum_mc_ratio_error = np_sum_mc_ratio * sqrt(pow(n_sum_mc_error/n_sum_mc,2) + pow(p_sum_mc_error/p_sum_mc,2));
//double np_sum_exp_ratio_error = np_sum_exp_ratio * sqrt(pow(n_sum_exp_error/n_sum_exp,2) + pow(p_sum_exp_error/p_sum_exp,2));
double np_red_cross_section_ratio_exp_err = np_red_cross_section_ratio_exp * sqrt(pow(R_sf_error/R_sf,2) + pow(np_red_cross_section_ratio_err/np_red_cross_section_ratio,2));
//double np_scale_ratio_error = np_scale_ratio * sqrt(pow(n_scale_error/n_scale,2) + pow(p_scale_error/p_scale,2));

//Temporary Legend change for nice presentation. Change back later
//Make a legend to hold all the information we care about
TLegend* legend = new TLegend(0.57, 0.5, 0.9, 0.9);
legend->AddEntry(hdx_data,"Data","l");
legend->AddEntry(fit, "Total Fit","l");
legend->AddEntry(bg,"Background","f l");
legend->AddEntry(hdx_mc_p_after,"Proton SIMC MC","f l");
legend->AddEntry(hdx_mc_n_after,"Neutron SIMC MC","f l");
legend->AddEntry((TObject*)0,Form("data N events : %0.0f",hdx_data->GetEntries()),"");
legend->AddEntry((TObject*)0,Form("n/p scale ratio R_{sf} : %0.5f #pm %0.5f",R_sf,R_sf_error),"");
//legend->AddEntry((TObject*)0,Form("MC n yield: %1.0f #pm %1.0f",n_sum_mc,n_sum_mc_error),"");
//legend->AddEntry((TObject*)0,Form("MC p yield: %1.0f #pm %1.0f",p_sum_mc,p_sum_mc_error),"");
//legend->AddEntry((TObject*)0,Form("MC n/p yield ratio: %0.3f #pm %0.3f",np_sum_mc_ratio,np_sum_mc_ratio_error),"");
//legend->AddEntry((TObject*)0,Form("Exp n yield: %1.0f #pm %1.0f",n_sum_exp,n_sum_exp_error),"");
//legend->AddEntry((TObject*)0,Form("Exp p yield: %1.0f #pm %1.0f",p_sum_exp,p_sum_exp_error),"");
//legend->AddEntry((TObject*)0,Form("Exp n/p yield ratio: %0.3f #pm %0.3f",np_sum_exp_ratio,np_sum_exp_ratio_error),"");
//legend->AddEntry((TObject*)0,Form("MC n red cross-section: %0.5f #pm %0.1f",n_red_cross_section,n_red_cross_section_err),"");
//legend->AddEntry((TObject*)0,Form("MC p red cross-section: %0.5f #pm %0.1f",p_red_cross_section,p_red_cross_section_err),"");
//legend->AddEntry((TObject*)0,Form("MC n/p RCS ratio: %0.3f #pm %0.3f",np_red_cross_section_ratio,np_red_cross_section_ratio_err),"");
//legend->AddEntry((TObject*)0,Form("Exp n red cross-section: %0.5f #pm %f",n_red_cross_section_exp,n_red_cross_section_exp_err),"");
//legend->AddEntry((TObject*)0,Form("Exp p red cross-section: %0.5f #pm %f",p_red_cross_section_exp,p_red_cross_section_exp_err),"");
//legend->AddEntry((TObject*)0,Form("Exp n/p RCS ratio: %0.4f #pm %0.4f",np_red_cross_section_ratio_exp,np_red_cross_section_ratio_exp_err),"");
//legend->AddEntry((TObject*)0,Form("dx shift pars, n/p : %0.3f / %0.3f ",n_shift,p_shift),"");
legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()),"");
legend->SetTextSize(0.025);
legend->Draw("same");

//Make the residual plot
//TH1D* residual = plots::makeResidualWithError("dx",hdx_data_bgsub,sumHist,true, false);
TH1D* residual_fit = plots::makeResidualWithError("dx_fit",hdx_data,fit);
TString residual_fit_name = residual_fit->GetName();
residual_fit->SetTitle("");
residual_fit->GetXaxis()->SetRangeUser(hcalfit_low,hcalfit_high);
double maxRange1 = residual_fit->GetMaximum();

double maxRange = 1.5*maxRange1;
residual_fit->GetYaxis()->SetRangeUser(-maxRange,maxRange);
//residual->SetLineColor(kBlue+2);
//residual->SetLineWidth(2);
residual_fit->SetStats(0);
residual_fit->GetYaxis()->SetTitle("Residuals");
residual_fit->GetYaxis()->SetTitleSize(0.07);
residual_fit->GetYaxis()->SetTitleOffset(0.5);
residual_fit->GetYaxis()->SetLabelSize(0.07);
residual_fit->SetLineColor(kOrange-3);
residual_fit->SetLineWidth(1);

//Make a legend to hold all the information we care about
TLegend* legend_res = new TLegend(0.75, 0.8, 0.9, 1.0);
//legend_res->AddEntry(residual,"Data-MC combine","l");
legend_res->AddEntry(residual_fit,"Data-Total Fit","l");
//Draw the residual on the second pad
pad2->cd();
//residual->Draw("hist");
gStyle->SetErrorX(0);
residual_fit->Draw("E1");
legend_res->SetTextSize(0.075);
//legend_res->Draw("same");
//Draw a line at y=0 just for better reading of the residual
TF1* zeroLine = new TF1("zeroLine","0",residual_fit->GetXaxis()->GetXmin(),residual_fit->GetXaxis()->GetXmax());
zeroLine->SetLineColor(kBlack);
zeroLine->SetLineStyle(2);
zeroLine->Draw("same");

report << "===================================================================" << endl;
report << "Report Information for " << fitName <<"." << endl;
report << "Neutron scale factor: " << n_scale << " +/- " << n_scale_error << endl;
report << "Proton scale factor: " << p_scale << " +/- " << p_scale_error << endl;
report << "n/p scale ratio, R_sf: " << R_sf << " +/- " << R_sf_error << endl;
report << "MC neutron reduced cross-section (before scale): " << n_red_cross_section << " +/- " << n_red_cross_section_err << endl;
report << "MC proton reduced cross-section (before scale): " << p_red_cross_section << " +/- " << p_red_cross_section_err << endl;
report << "MC n/p reduced cross-section ratio (before scale): " << np_red_cross_section_ratio << " +/- " << np_red_cross_section_ratio_err << endl;
//report << "MC neutron yield (before scale): " << n_sum_mc << " +/- " << n_sum_mc_error << endl;
//report << "MC proton  yield (before scale): " << p_sum_mc << " +/- " << p_sum_mc_error << endl;
//report << "MC n/p yield ratio (before scale): " << np_sum_mc_ratio << " +/- " << np_sum_mc_ratio_error << endl;
report << "Experiment neutron reduced cross-section (after scale): " << n_red_cross_section_exp << " +/- " << n_red_cross_section_exp_err << endl;
report << "Experiment proton reduced cross-section (after scale): " << p_red_cross_section_exp << " +/- " << p_red_cross_section_exp_err << endl;
report << "Experiment n/p RCS ratio (after scale): " << np_red_cross_section_ratio_exp << " +/- " << np_red_cross_section_ratio_exp_err << endl;
//report << "Experiment neutron yield (after scale): " << n_sum_exp << " +/- " << n_sum_exp_error << endl;
//report << "Experiment proton yield (after scale): " << p_sum_exp << " +/- " << p_sum_exp_error << endl;
//report << "Experiment n/p yield ratio (after scale): " << np_sum_exp_ratio << " +/- " << np_sum_exp_ratio_error << endl;
report << "neutron dx shift parameter: " << n_shift << endl;
report << "proton dx shift parameter: " << p_shift << endl;
report << "Chi^2/ndf: " << fit->GetChisquare() << "/" << fit->GetNDF() << endl;
report << "===================================================================" << endl << endl;


return myTotCan;
}


//Make a canvas that displays the Data and MC dx plot being compared. With background already subtracked. Display relevant yield and ratio information
TCanvas* plotDataMCFitsResiduals_NoBG(fit_histogram* myFitHisto,const char *name, bool shiftfit, std::ofstream& report, calc_FFs_RCS_obj daFFInfo){

//Initialize some information from the fit_histogram object
TH1D* hdx_data = myFitHisto->get_hist_data();
TH1D* hdx_mc_p = myFitHisto->get_hist_p();
TH1D* hdx_mc_n = myFitHisto->get_hist_n();
string fitName = myFitHisto->get_fitName();
string fitType = myFitHisto->get_fitType();
vector<pair<double,double>> params = myFitHisto->get_fitParamsErrs();
double hcalfit_low = myFitHisto->get_xMin();
double hcalfit_high = myFitHisto->get_xMax();
double n_shift = myFitHisto->get_shift_n();
double p_shift = myFitHisto->get_shift_p();

//Recreate the fit
TF1* fit;

        //This conditional works but we will need to modify it if we need any new types of functions
        if(fitType == "fitFullNoBG"){
        fit = new TF1("fit",myFitHisto,&fit_histogram::fitFullNoBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFullNoBG");
        }else if(fitType =="fitFullShiftNoBG"){
        fit = new TF1("fit",myFitHisto,&fit_histogram::fitFullShiftNoBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFullShiftNoBG");
	}else{
        cout << "The plot function you are trying to implement " << fitType << " is no good! Figure it out now!" << endl;
        fit = new TF1("fit",myFitHisto,&fit_histogram::fitFullShiftNoBG,hcalfit_low,hcalfit_high,params.size(),"fit_histogram","fitFullShiftNoBG");
        }
        //file the fit with the parameters
        for(int i=0; i < params.size();++i){
        fit->SetParameter(i,params[i].first);
        fit->SetParError(i,params[i].second);
        }
fit->SetChisquare(myFitHisto->get_ChiSq());
fit->SetNDF(myFitHisto->get_NDF());


//Make clones of the mc histograms
//TH1D *hdx_mc_p_before = (TH1D*)(hdx_mc_p->Clone("hdx_mc_p_before"));
TH1D *hdx_mc_p_after = (TH1D*)(hdx_mc_p->Clone("hdx_mc_p_after"));

//TH1D *hdx_mc_n_before = (TH1D*)(hdx_mc_n->Clone("hdx_mc_n_before"));
TH1D *hdx_mc_n_after = (TH1D*)(hdx_mc_n->Clone("hdx_mc_n_after"));

//Create a canvas we will use to store information
TCanvas *myTotCan = new TCanvas(name,Form("Data with Fits and Residuals %s",fitName.c_str()),1600,1200);
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
hdx_data->SetTitle(Form("dx, %s;m",fitName.c_str()));
hdx_data->SetStats(0);
hdx_data->Draw("hist");
//fit the data dx
fit->SetLineColor(kMagenta);
fit->SetLineWidth(2);
fit->Draw("same");
//proton MC dx
        if(shiftfit){
	//hdx_mc_p_before = plots::shiftHistogramX(hdx_mc_p_before,p_shift);
        hdx_mc_p_after = plots::shiftHistogramX(hdx_mc_p_after,p_shift);
	}

//Yield info and integrals come before scaling to get sigma_n/sigma_p from MC. Need to further scale.

//yield info for protons
//double p_sum_mc_error;
//double p_sum_mc = hdx_mc_p_before->IntegralAndError(0,hdx_mc_p_before->GetNbinsX()+1,p_sum_mc_error,"");

//proton scale info
double p_scale = myFitHisto->get_scale_p();
double p_scale_error = myFitHisto->get_scale_p_err();

//neutron scale info
//double n_scale = params[1].first;
//double n_scale_error = params[1].second;
//Get Rsf from the fit parameters
double R_sf = myFitHisto->get_Rsf();
double R_sf_error = myFitHisto->get_Rsf_err();

//scale info for neutrons
double n_scale = myFitHisto->get_scale_n();
double n_scale_error = myFitHisto->get_scale_n_err();

hdx_mc_p_after->Scale(p_scale);
hdx_mc_p_after->SetLineColor(kRed);
hdx_mc_p_after->SetFillColorAlpha(kRed-9,0.5);
hdx_mc_p_after->SetFillStyle(1001);
hdx_mc_p_after->SetLineWidth(2);
hdx_mc_p_after->Draw("hist same");
//neutron MC dx
        if(shiftfit){
        //hdx_mc_n_before = plots::shiftHistogramX(hdx_mc_n_before,n_shift);
	hdx_mc_n_after = plots::shiftHistogramX(hdx_mc_n_after,n_shift);
        }
//yield info for neutrons
//double n_sum_mc_error;
//double n_sum_mc = hdx_mc_n_before->IntegralAndError(0,hdx_mc_n_before->GetNbinsX()+1,n_sum_mc_error,"");
	
hdx_mc_n_after->Scale(n_scale);
hdx_mc_n_after->SetLineColor(kBlue);
hdx_mc_n_after->SetFillColorAlpha(kBlue-9,0.5);
hdx_mc_n_after->SetFillStyle(1001);
hdx_mc_n_after->SetLineWidth(2);
hdx_mc_n_after->Draw("hist same");

//The yield method is being turned off. As it includes detection efficiency and and nuclear and radiative effects baked in
//Yield info and integrals come after scaling to get sigma_n/sigma_p for experiment. Don't need further scale.

//yield info for protons
//double p_sum_exp_error;
//double p_sum_exp = hdx_mc_p_after->IntegralAndError(0,hdx_mc_p_after->GetNbinsX()+1,p_sum_exp_error,"");

//yield info for neutrons
//double n_sum_exp_error;
//double n_sum_exp = hdx_mc_n_after->IntegralAndError(0,hdx_mc_n_after->GetNbinsX()+1,n_sum_exp_error,"");

//Just MC
//Instead calculate the reduced cross-section for protons and neutrons using the FF parameterizations in SIMC. Here we are using a class developed by me to do this purpose
double n_red_cross_section = daFFInfo.get_reduced_cross_section_n();
double n_red_cross_section_err = daFFInfo.get_reduced_cross_section_n_err();

double p_red_cross_section = daFFInfo.get_reduced_cross_section_p();
double p_red_cross_section_err = daFFInfo.get_reduced_cross_section_p_err();
//Then do the n/p ratio. This is what we call Rsim
double np_red_cross_section_ratio = daFFInfo.get_reduced_cross_section_ratio();
double np_red_cross_section_ratio_err = daFFInfo.get_reduced_cross_section_ratio_err();

//Scale the MC information by the n and p scale factor or Rsf
double n_red_cross_section_exp = n_red_cross_section*n_scale;
double n_red_cross_section_exp_err = n_red_cross_section_exp * sqrt(pow(n_scale_error/n_scale,2) + pow(n_red_cross_section_err/n_red_cross_section,2));

double p_red_cross_section_exp = p_red_cross_section*p_scale;
double p_red_cross_section_exp_err = p_red_cross_section_exp * sqrt(pow(p_scale_error/p_scale,2) + pow(p_red_cross_section_err/p_red_cross_section,2));

//get the ratios and the errors
//double np_sum_mc_ratio = n_sum_mc/p_sum_mc;
//double np_sum_exp_ratio = n_sum_exp/p_sum_exp;

//get the exp np ratio
double np_red_cross_section_ratio_exp = R_sf*np_red_cross_section_ratio;
//get the np_scale ratio, we call it Rsf
//double np_scale_ratio = n_scale/p_scale;

//double np_sum_mc_ratio_error = np_sum_mc_ratio * sqrt(pow(n_sum_mc_error/n_sum_mc,2) + pow(p_sum_mc_error/p_sum_mc,2));
//double np_sum_exp_ratio_error = np_sum_exp_ratio * sqrt(pow(n_sum_exp_error/n_sum_exp,2) + pow(p_sum_exp_error/p_sum_exp,2));
//Get the corresponding errors
double np_red_cross_section_ratio_exp_err = np_red_cross_section_ratio_exp * sqrt(pow(R_sf_error/R_sf,2) + pow(np_red_cross_section_ratio_err/np_red_cross_section_ratio,2));
//double np_scale_ratio_error = np_scale_ratio * sqrt(pow(n_scale_error/n_scale,2) + pow(p_scale_error/p_scale,2));

//Make a legend to hold all the information we care about
TLegend* legend = new TLegend(0.57, 0.5, 0.9, 0.9);
legend->AddEntry(hdx_data,"Data","l");
legend->AddEntry(fit, "Total Fit","l");
legend->AddEntry(hdx_mc_p_after,"Proton SIMC MC","f l");
legend->AddEntry(hdx_mc_n_after,"Neutron SIMC MC","f l");
legend->AddEntry((TObject*)0,Form("data N events : %0.0f",hdx_data->GetEntries()),"");
legend->AddEntry((TObject*)0,Form("n/p scale ratio R_{sf} : %0.4f #pm %0.4f",R_sf,R_sf_error),"");
//legend->AddEntry((TObject*)0,Form("MC n yield: %1.0f #pm %1.0f",n_sum_mc,n_sum_mc_error),"");
//legend->AddEntry((TObject*)0,Form("MC p yield: %1.0f #pm %1.0f",p_sum_mc,p_sum_mc_error),"");
//legend->AddEntry((TObject*)0,Form("MC n/p yield ratio: %0.3f #pm %0.3f",np_sum_mc_ratio,np_sum_mc_ratio_error),"");
//legend->AddEntry((TObject*)0,Form("Exp n yield: %1.0f #pm %1.0f",n_sum_exp,n_sum_exp_error),"");
//legend->AddEntry((TObject*)0,Form("Exp p yield: %1.0f #pm %1.0f",p_sum_exp,p_sum_exp_error),"");
//legend->AddEntry((TObject*)0,Form("Exp n/p yield ratio: %0.3f #pm %0.3f",np_sum_exp_ratio,np_sum_exp_ratio_error),"");
legend->AddEntry((TObject*)0,Form("MC n red cross-section: %0.5f #pm %0.1f",n_red_cross_section,n_red_cross_section_err),"");
legend->AddEntry((TObject*)0,Form("MC p red cross-section: %0.5f #pm %0.1f",p_red_cross_section,p_red_cross_section_err),"");
legend->AddEntry((TObject*)0,Form("MC n/p RCS ratio: %0.3f #pm %0.3f",np_red_cross_section_ratio,np_red_cross_section_ratio_err),"");
legend->AddEntry((TObject*)0,Form("Exp n red cross-section: %0.5f #pm %0.6f",n_red_cross_section_exp,n_red_cross_section_exp_err),"");
legend->AddEntry((TObject*)0,Form("Exp p red cross-section: %0.5f #pm %0.6f",p_red_cross_section_exp,p_red_cross_section_exp_err),"");
legend->AddEntry((TObject*)0,Form("Exp n/p RCS ratio: %0.4f #pm %0.4f",np_red_cross_section_ratio_exp,np_red_cross_section_ratio_exp_err),"");
legend->AddEntry((TObject*)0,Form("dx shift pars, n/p : %0.3f / %0.3f ",params[3].first,params[2].first),"");
legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()),"");
legend->SetTextSize(0.025);
legend->Draw("same");

//Residual plot
//Add the two MC histograms together for direct comparison
//TH1D* sumHist = new TH1D("sumHist2","dx for both MC protons and neutrons",hdx_mc_p->GetNbinsX(),hdx_mc_p->GetXaxis()->GetXmin(),hdx_mc_p->GetXaxis()->GetXmax());

//sumHist->Add(hdx_mc_p,hdx_mc_n);
TH1D* residual_fit = plots::makeResidualWithError("dx_fit",hdx_data,fit);

residual_fit->SetTitle("");
residual_fit->GetXaxis()->SetRangeUser(hcalfit_low,hcalfit_high);
double maxRange1 = residual_fit->GetMaximum();
double maxRange = 1.5*maxRange1;
residual_fit->GetYaxis()->SetRangeUser(-maxRange,maxRange);

residual_fit->GetYaxis()->SetTitle("Residuals");
residual_fit->GetYaxis()->SetTitleSize(0.07);
residual_fit->GetYaxis()->SetTitleOffset(0.5);
residual_fit->GetYaxis()->SetLabelSize(0.07);
residual_fit->SetLineColor(kOrange-3);
residual_fit->SetLineWidth(1);
residual_fit->SetStats(0);

//Make a legend to hold all the information we care about
TLegend* legend_res = new TLegend(0.75, 0.8, 0.9, 1.0);
legend_res->AddEntry(residual_fit,"Data-Total Fit","l");

//Draw the residual on the second pad
pad2->cd();
residual_fit->Draw("E1 ");
legend_res->SetTextSize(0.075);
legend_res->Draw("same");

//Draw a line at y=0 just for better reading of the residual
TF1* zeroLine = new TF1("zeroLine","0",residual_fit->GetXaxis()->GetXmin(),residual_fit->GetXaxis()->GetXmax());
zeroLine->SetLineColor(kBlack);
zeroLine->SetLineStyle(2);
zeroLine->Draw("same");

return myTotCan;

}
//Function to plot the background histogram with corresponding fit function. Then a residual between the fit and histogram.
TCanvas* plotBGResiduals(TH1D* hdx_data, TH1D* hdx_mc_p, TH1D* hdx_mc_n, TF1* bg,const char *name,const char *fitName, const char* fitType, const vector<pair<double,double>> params,double hcalfit_low, double hcalfit_high,bool shiftfit){

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

//proton MC dx
	if(shiftfit){
        hdx_mc_p = plots::shiftHistogramX(hdx_mc_p,params[2].first);
        }
hdx_mc_p->Scale(params[0].first);
//neutron MC dx
	if(shiftfit){
        hdx_mc_n = plots::shiftHistogramX(hdx_mc_n,params[3].first);
        }
hdx_mc_n->Scale(params[1].first*params[0].first);

//Add the two MC histograms together for direct comparison
TH1D* sumHist = new TH1D(Form("sumHist %s %s",fitType, fitName),"dx for both MC protons and neutrons",hdx_mc_p->GetNbinsX(),hdx_mc_p->GetXaxis()->GetXmin(),hdx_mc_p->GetXaxis()->GetXmax());
sumHist->Add(hdx_mc_p,hdx_mc_n);
//Subtract the MC info from data to get hist that looks like background
TH1D* hdx_bg = plots::subtractHist(hdx_data,sumHist);

//Draw the histogram with fits on the first pad
pad1->cd();
hdx_bg->SetLineColor(kBlack);
hdx_bg->SetLineWidth(1);
hdx_bg->SetTitle(Form("dx BG, %s;m",fitName));
hdx_bg->SetStats(0);
hdx_bg->Draw("hist");

bg->SetLineColor(kGreen);
bg->SetLineWidth(2);
bg->Draw("same");

//Make a legend to hold all the information we care about
TLegend* legend = new TLegend(0.8, 0.8, 0.9, 0.9);
legend->AddEntry(hdx_bg,"BG Hist","l");
legend->AddEntry(bg, "BG Fit","l");
legend->SetTextSize(0.04);
legend->Draw("same");

//Draw the residual on the second pad
pad2->cd();

TH1D* residual_fit = plots::makeResidualWithError("dx_bg_resid",hdx_bg,bg);
residual_fit->SetTitle("");
residual_fit->GetXaxis()->SetRangeUser(hcalfit_low,hcalfit_high);
double maxRange1 = residual_fit->GetMaximum();
double maxRange = 1.75*maxRange1;
residual_fit->GetYaxis()->SetRangeUser(-maxRange,maxRange);

residual_fit->GetYaxis()->SetTitle("Residuals");
residual_fit->GetYaxis()->SetTitleSize(0.07);
residual_fit->GetYaxis()->SetTitleOffset(0.5);
residual_fit->GetYaxis()->SetLabelSize(0.07);
residual_fit->SetLineColor(kBlue+2);
residual_fit->SetLineWidth(1);
residual_fit->SetStats(0);

TLegend* legend_res = new TLegend(0.75, 0.8, 0.9, 1.0);
legend_res->AddEntry(residual_fit,"BGHist-BG Fit","l");

residual_fit->Draw("E1 ");
legend_res->SetTextSize(0.075);
legend_res->Draw("same");

//Draw a line at y=0 just for better reading of the residual
TF1* zeroLine = new TF1("zeroLine","0",residual_fit->GetXaxis()->GetXmin(),residual_fit->GetXaxis()->GetXmax());
zeroLine->SetLineColor(kBlack);
zeroLine->SetLineStyle(2);
zeroLine->Draw("same");


return myTotCan;
}


//Function to plot the fit over the HCal position dependent efficiency histogram. Should work for both x and y directions
TCanvas* plotHCalEff(TH1D* h_eff_plot,const char *can_name, const char *name, const char *fitName, double fit_low, double fit_high, vector<TLine*> Lines_Fid, TString target, TString spot_choice,int sbs_field){

	gStyle->SetGridStyle(1);
	//Create a canvas we will use to store information
	TCanvas *myTotCan = new TCanvas(can_name,Form("HCal Efficiency: %s",fitName),1600,1200);

	//Make clones of the efficiency histograms so we can fit them and evaluate the average relative efficiency
	TH1D *heff_plot_clone = (TH1D*) (h_eff_plot->Clone(name));
	heff_plot_clone->GetYaxis()->SetTitle("Relative Rate");
	//sets grid divisions
	string mystring(name);

	if(mystring.find("xexpect") != string::npos ){
	 heff_plot_clone->GetXaxis()->SetNdivisions(25);
        heff_plot_clone->GetYaxis()->SetNdivisions(10);
	}else if(mystring.find("yexpect") != string::npos){
	 heff_plot_clone->GetXaxis()->SetNdivisions(12);
        heff_plot_clone->GetYaxis()->SetNdivisions(10);
	}else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

	if(target == "LD2" && sbs_field > 0){
		if(spot_choice == "proton"){
		heff_plot_clone->GetYaxis()->SetRangeUser(0.0,0.45);
		}else if(spot_choice == "neutron"){
		heff_plot_clone->GetYaxis()->SetRangeUser(0.0,0.25);
		}
	}

	//Make the fits for the histogram
	TF1 *eff_fit = new TF1(fitName,"pol0",fit_low,fit_high);
	eff_fit->SetNpx(1000);

	//Fit the histo
	heff_plot_clone->Fit(eff_fit,"QR");
	
	//Draw the histogram and fit
        heff_plot_clone->SetMarkerStyle(20);
	//enables grid
	myTotCan->SetGrid();
	heff_plot_clone->Draw("E");
	eff_fit->Draw("same");
	
	
        if(mystring.find("xexpect") != string::npos ){
        Lines_Fid[0]->Draw("same");
        Lines_Fid[1]->Draw("same");

        }else if(mystring.find("yexpect") != string::npos){
        Lines_Fid[2]->Draw("same");
        Lines_Fid[3]->Draw("same");

        }else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

	//Make a legend to hold all the information we care about
	TLegend* legend = new TLegend(0.65, 0.2, 0.9, 0.3);	
	legend->AddEntry((TObject*)0,Form("p0: %0.4f #pm %0.4f",eff_fit->GetParameter(0),eff_fit->GetParError(0)),"");
	legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf: %0.3f/%d",eff_fit->GetChisquare(),eff_fit->GetNDF()),"");
	legend->SetFillStyle(0);
	legend->Draw("same");
	

	return myTotCan;
}

TCanvas* plotHCalEff(TH1D* h_eff_plot,const char *can_name, const char *name, const char *fitName, double fit_low, double fit_high, vector<TLine*> Lines_Fid){
gStyle->SetGridStyle(1);
        //Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("HCal Efficiency: %s",fitName),1600,1200);

        //Make clones of the efficiency histograms so we can fit them and evaluate the average relative efficiency
        TH1D *heff_plot_clone = (TH1D*) (h_eff_plot->Clone(name));
        heff_plot_clone->GetYaxis()->SetTitle("Relative Rate");
        //sets grid divisions
        string mystring(name);

        if(mystring.find("xexpect") != string::npos ){
         heff_plot_clone->GetXaxis()->SetNdivisions(25);
        heff_plot_clone->GetYaxis()->SetNdivisions(10);
        }else if(mystring.find("yexpect") != string::npos){
         heff_plot_clone->GetXaxis()->SetNdivisions(12);
        heff_plot_clone->GetYaxis()->SetNdivisions(10);
        }else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

        //Make the fits for the histogram
        TF1 *eff_fit = new TF1(fitName,"pol0",fit_low,fit_high);
        eff_fit->SetNpx(1000);

        //Fit the histo
        heff_plot_clone->Fit(eff_fit,"QR");

        //Draw the histogram and fit
        heff_plot_clone->SetMarkerStyle(20);
//enables grid
        myTotCan->SetGrid();
        heff_plot_clone->Draw("E");
        eff_fit->Draw("same");


        if(mystring.find("xexpect") != string::npos ){
        Lines_Fid[0]->Draw("same");
        Lines_Fid[1]->Draw("same");

        }else if(mystring.find("yexpect") != string::npos){
        Lines_Fid[2]->Draw("same");
        Lines_Fid[3]->Draw("same");

        }else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

        //Make a legend to hold all the information we care about
        TLegend* legend = new TLegend(0.65, 0.2, 0.9, 0.3);
        legend->AddEntry((TObject*)0,Form("p0: %0.4f #pm %0.4f",eff_fit->GetParameter(0),eff_fit->GetParError(0)),"");
        legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf: %0.3f/%d",eff_fit->GetChisquare(),eff_fit->GetNDF()),"");
        legend->SetFillStyle(0);
        legend->Draw("same");


        return myTotCan;
}


//A function to compare a histogram under different sets of cuts on the same scale.
TCanvas* plot_Comp(TH1D* plot_nocut,TH1D* plot_cut, const char *can_name, const char *name_nocut, const char *name_cut){

	//Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("Plot Comparison: %s and %s",plot_nocut->GetName(),plot_cut->GetName()),1600,1200);

	//Make clones of the histograms
	TH1D* plot_nocut_clone = (TH1D*) (plot_nocut->Clone(name_nocut));
	TH1D* plot_cut_clone = (TH1D*) (plot_cut->Clone(name_cut));
	plot_nocut_clone->SetMarkerStyle(20);
	plot_cut_clone->SetMarkerStyle(20);
	plot_nocut_clone->SetLineColor(kBlack);
	plot_nocut_clone->SetMarkerColor(kBlack);
	plot_cut_clone->SetLineColor(kRed);
	plot_cut_clone->SetMarkerColor(kRed);
	
	plot_nocut_clone->Draw("P");
	plot_cut_clone->Draw("same");

	//Make a legend to hold all the information we care about
	TLegend* legend = new TLegend(0.35, 0.35, 0.55, 0.45);
	legend->AddEntry(plot_nocut_clone,"Only E-arm cuts","f l");
	legend->AddEntry(plot_cut_clone," E-arm + H-arm cuts","f l");
	legend->SetFillStyle(0);
	legend->Draw("same");

	return myTotCan;
}


TCanvas* plot_Comp_1DEff(TH1D* plot_nocut,TH1D* plot_cut, const char *can_name, const char *name_nocut, const char *name_cut,const char *other_name,vector<TLine*> Lines_Fid, TString target){
	gStyle->SetGridStyle(1);
	//Create a canvas we will use to store information
	TCanvas *myTotCan = new TCanvas(can_name,Form("Plot Comparison: %s and %s",plot_nocut->GetName(),plot_cut->GetName()),1600,1200);

        //Make clones of the histograms
        TH1D* plot_nocut_clone = (TH1D*) (plot_nocut->Clone(name_nocut));
        TH1D* plot_cut_clone = (TH1D*) (plot_cut->Clone(name_cut));
        plot_nocut_clone->SetMarkerStyle(20);
        plot_cut_clone->SetMarkerStyle(20);
        plot_nocut_clone->SetLineColor(kBlack);
        plot_nocut_clone->SetMarkerColor(kBlack);
        plot_cut_clone->SetLineColor(kRed);
        plot_cut_clone->SetMarkerColor(kRed);
	//sets grid divisions
        string mystring(other_name);

        if(mystring.find("xexpect") != string::npos ){
        plot_nocut_clone->GetXaxis()->SetNdivisions(25);
        plot_nocut_clone->GetYaxis()->SetNdivisions(10);
        }else if(mystring.find("yexpect") != string::npos){
        plot_nocut_clone->GetXaxis()->SetNdivisions(12);
        plot_nocut_clone->GetYaxis()->SetNdivisions(10);
        }else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

	if(target == "LD2"){
                plot_nocut_clone->GetYaxis()->SetRangeUser(0.0,0.45);
		plot_cut_clone->GetYaxis()->SetRangeUser(0.0,0.45);
        }

	//enables grid
        myTotCan->SetGrid();
        plot_nocut_clone->Draw("P");
        plot_cut_clone->Draw("same");

        if(mystring.find("xexpect") != string::npos ){
        Lines_Fid[0]->Draw("same");
        Lines_Fid[1]->Draw("same");

        }else if(mystring.find("yexpect") != string::npos){
        Lines_Fid[2]->Draw("same");
        Lines_Fid[3]->Draw("same");

        }else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

        //Make a legend to hold all the information we care about
        TLegend* legend = new TLegend(0.7, 0.35, 0.9, 0.45);
        legend->AddEntry(plot_nocut_clone,name_nocut,"f l");
        legend->AddEntry(plot_cut_clone,name_cut,"f l");
        legend->SetFillStyle(0);
        legend->Draw("same");

        return myTotCan;


}

//A function to plot the 2D efficiency map on a canvas
TCanvas* plot_HCalEffMap(TH2D* eff_map,const char *can_name, const char *name, TString target,TString spot_choice){

	//Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("%s",name),1600,1200);

	//Make a clone of the 2D efficiency histogram to be able to manipulate it
	TH2D* eff_map_clone = (TH2D*) (eff_map->Clone(name));
	gStyle->SetPalette(55);
	gStyle->SetNumberContours(256);
	eff_map_clone->SetStats(0);

	/*if(target == "LD2"){
                if(spot_choice == "proton"){
                eff_map_clone->GetZaxis()->SetRangeUser(0.0,0.45);
                }else if(spot_choice == "neutron"){
                eff_map_clone->GetZaxis()->SetRangeUser(0.0,0.2);
		}
        }*/
	eff_map_clone->Draw("colz");

	return myTotCan;
}

TCanvas* plot_HCalEffMap_1D(TH1D* eff_vs_expect,const char *can_name, const char *name, const char *label,vector<TLine*> Lines_Fid){
	gStyle->SetGridStyle(1);
	//Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("%s",label),1600,1200);
	
	//Make a clone of the 1D efficiency histogram to be able to manipulate it
	TH1D* eff_vs_expect_clone = (TH1D*) (eff_vs_expect->Clone(name));
	eff_vs_expect_clone->SetTitle(label);
	eff_vs_expect_clone->SetStats(0);
	double expect_max = eff_vs_expect_clone->GetMaximum() + 0.6;

	eff_vs_expect_clone->GetYaxis()->SetRangeUser(0.0,expect_max);
	//eff_vs_expect_clone->GetYaxis()->SetRangeUser(0.0,1.5);
	//eff_vs_expect_clone->GetYaxis()->SetTitle("Relative Rate Ratio");
	
	string mystring(name);
        if(mystring.find("xexpect") != string::npos ){
        eff_vs_expect_clone->GetXaxis()->SetNdivisions(25);
        eff_vs_expect_clone->GetYaxis()->SetNdivisions(10);
        }else if(mystring.find("yexpect") != string::npos){
        eff_vs_expect_clone->GetXaxis()->SetNdivisions(12);
        eff_vs_expect_clone->GetYaxis()->SetNdivisions(10);0
        }else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }	
	
	myTotCan->SetGrid();
	eff_vs_expect_clone->Draw();

	if(mystring.find("xexpect") != string::npos ){
	Lines_Fid[0]->Draw("same");
        Lines_Fid[1]->Draw("same");
	
	}else if(mystring.find("yexpect") != string::npos){
	Lines_Fid[2]->Draw("same");
        Lines_Fid[3]->Draw("same");

	}else{
	cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
	}

	return myTotCan;
}

//A function to plot the 2D efficiency map on a canvas
TCanvas* plot_HCalEffMap_Comp(TH2D* eff_map,const char *can_name, const char *name){

        //Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("%s",name),1600,1200);

        //Make a clone of the 2D efficiency histogram to be able to manipulate it
        TH2D* eff_map_clone = (TH2D*) (eff_map->Clone(name));
        gStyle->SetPalette(55);
        gStyle->SetNumberContours(256);
	eff_map_clone->SetTitle(name);
	eff_map_clone->SetStats(0);
        eff_map_clone->GetZaxis()->SetRangeUser(0.0,2.0);
	eff_map_clone->Draw("colz");

        return myTotCan;
}

//A function to plot the 2D efficiency map on a canvas with some HCal info overlaid
TCanvas* plot_HCalEffMap_overlay(TH2D* eff_map,const char *can_name, const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_Fid){

        //Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("%s",name),1600,1200);

        //Make a clone of the 2D efficiency histogram to be able to manipulate it
        TH2D* eff_map_clone = (TH2D*) (eff_map->Clone(name));
        gStyle->SetPalette(55);
        gStyle->SetNumberContours(256);
        
	eff_map_clone->SetStats(0);
        eff_map_clone->Draw("colz");
        
	Lines_pos[0]->Draw("same");
	Lines_pos[1]->Draw("same");
	Lines_pos[2]->Draw("same");
	Lines_pos[3]->Draw("same");

	Lines_Fid[0]->Draw("same");
	Lines_Fid[1]->Draw("same");
	Lines_Fid[2]->Draw("same");
	Lines_Fid[3]->Draw("same");

	//make legend
	auto legend1 = new TLegend(0.1,0.8,0.3,0.9);
	legend1->SetTextSize(0.03);
	legend1->AddEntry(Lines_pos[0],"HCal Boundary","l");
	legend1->AddEntry(Lines_Fid[0],"Fiducial Region","l");
	legend1->Draw();

	return myTotCan;
}

//A function to plot the 2D efficiency map on a canvas with some HCal info overlaid
TCanvas* plot_HCalEffMap_overlay_Comp(TH2D* eff_map,const char *can_name, const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_Fid){

        //Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("%s",name),1600,1200);

        //Make a clone of the 2D efficiency histogram to be able to manipulate it
        TH2D* eff_map_clone = (TH2D*) (eff_map->Clone(name));
        gStyle->SetPalette(55);
        gStyle->SetNumberContours(256);
	eff_map_clone->SetTitle(name);
        eff_map_clone->SetStats(0);
        eff_map_clone->GetZaxis()->SetRangeUser(0.0,2.0);
	eff_map_clone->Draw("colz");

        Lines_pos[0]->Draw("same");
        Lines_pos[1]->Draw("same");
        Lines_pos[2]->Draw("same");
        Lines_pos[3]->Draw("same");

        Lines_Fid[0]->Draw("same");
        Lines_Fid[1]->Draw("same");
        Lines_Fid[2]->Draw("same");
        Lines_Fid[3]->Draw("same");

        //make legend
        auto legend1 = new TLegend(0.1,0.8,0.3,0.9);
        legend1->SetTextSize(0.03);
        legend1->AddEntry(Lines_pos[0],"HCal Boundary","l");
        legend1->AddEntry(Lines_Fid[0],"Fiducial Region","l");
        legend1->Draw();

        return myTotCan;
}

TCanvas* plot_Comp_1DEff_Scale(TH1D* plot_nocut,TH1D* plot_cut, const char *can_name, const char *name_nocut, const char *name_cut,const char *other_name,vector<TLine*> Lines_Fid, TString target){
        gStyle->SetGridStyle(1);
        //Create a canvas we will use to store information
        TCanvas *myTotCan = new TCanvas(can_name,Form("Plot Comparison: %s and %s",plot_nocut->GetName(),plot_cut->GetName()),1600,1200);

        //Make clones of the histograms
        TH1D* plot_nocut_clone = (TH1D*) (plot_nocut->Clone(name_nocut));
        TH1D* plot_cut_clone = (TH1D*) (plot_cut->Clone(name_cut));
        plot_nocut_clone->SetMarkerStyle(20);
        plot_cut_clone->SetMarkerStyle(20);
        plot_nocut_clone->SetLineColor(kBlack);
        plot_nocut_clone->SetMarkerColor(kBlack);
        plot_cut_clone->SetLineColor(kRed);
        plot_cut_clone->SetMarkerColor(kRed);
        //sets grid divisions
        string mystring(other_name);

	double plot_nocut_clone_int;
        double plot_cut_clone_int;
        if(mystring.find("xexpect") != string::npos ){
        plot_nocut_clone->GetXaxis()->SetNdivisions(25);
        plot_nocut_clone->GetYaxis()->SetNdivisions(10);
	int bin_low_1 = plot_nocut_clone->FindBin(-0.5);
	//int bin_low_1 = plot_nocut_clone->FindBin(-1.5);
	int bin_high_1 = plot_nocut_clone->FindBin(0.5);
	//int bin_high_1 = plot_nocut_clone->FindBin(0.0);
	int bin_low_2 = plot_cut_clone->FindBin(-0.5);
        //int bin_low_2 = plot_cut_clone->FindBin(-1.5);
       	int bin_high_2 = plot_cut_clone->FindBin(0.5);
	//int bin_high_2 = plot_cut_clone->FindBin(0.0);
	plot_nocut_clone_int = plot_nocut_clone->Integral(bin_low_1,bin_high_1);
        plot_cut_clone_int = plot_cut_clone->Integral(bin_low_2,bin_high_2);
        }else if(mystring.find("yexpect") != string::npos){
        plot_nocut_clone->GetXaxis()->SetNdivisions(12);
        plot_nocut_clone->GetYaxis()->SetNdivisions(10);
        int bin_low_1 = plot_nocut_clone->FindBin(-0.3);
        int bin_high_1 = plot_nocut_clone->FindBin(0.3);
        int bin_low_2 = plot_cut_clone->FindBin(-0.3);
        int bin_high_2 = plot_cut_clone->FindBin(0.3);
        plot_nocut_clone_int = plot_nocut_clone->Integral(bin_low_1,bin_high_1);
        plot_cut_clone_int = plot_cut_clone->Integral(bin_low_2,bin_high_2);
	}else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

        //Scale the red histogram up to the black to compare the two histogram
	double my_scale_fac = plot_nocut_clone_int/plot_cut_clone_int;
	
	plot_cut_clone->Scale(my_scale_fac);

	if(target == "LD2"){
                plot_nocut_clone->GetYaxis()->SetRangeUser(0.0,0.45);
                plot_cut_clone->GetYaxis()->SetRangeUser(0.0,0.45);
        }

        //enables grid
        myTotCan->SetGrid();
        plot_nocut_clone->Draw("P");
        plot_cut_clone->Draw("same");

        if(mystring.find("xexpect") != string::npos ){
        Lines_Fid[0]->Draw("same");
        Lines_Fid[1]->Draw("same");

        }else if(mystring.find("yexpect") != string::npos){
	Lines_Fid[2]->Draw("same");
        Lines_Fid[3]->Draw("same");

        }else{
        cout << "We did not find the substring properly in plot_HCalEffMap_1D. Figure it out!" << endl;
        }

        //Make a legend to hold all the information we care about
        TLegend* legend = new TLegend(0.7, 0.35, 0.9, 0.45);
        legend->AddEntry(plot_nocut_clone,name_nocut,"f l");
        legend->AddEntry(plot_cut_clone,Form("%s_Scaled",name_cut),"f l");
        legend->SetFillStyle(0);
        legend->Draw("same");

        return myTotCan;


}
		

}//end namespace
