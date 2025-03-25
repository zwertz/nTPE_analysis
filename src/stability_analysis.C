//Author: Ezekiel Wertz
//This is a companion class to cutvar.h and cutvar.C. I ran into a challenge where I needed to be able to manipulate the information for all 3 histograms for data, mc p, and mc n.Building this into cutvar seemed difficult. So instead this class will take output, after some processing, from cutvar. Then do further cut stability study manipulation/evaluation here. This class should handle all information from data and MC at once. Will be modified as stability and systematic studies progress.

#include "../include/stability_analysis.h"

stability_analysis::stability_analysis(cutvar &dataVar,cutvar &MC_pVar,cutvar &MC_nVar, vector<TH1D*> daDataHisto, vector<TH1D*> daMCPHisto, vector<TH1D*> daMCNHisto,const char* fit_type){

	if(dataVar.getCutVar() != MC_pVar.getCutVar() || dataVar.getCutVar() != MC_nVar.getCutVar()){
	cout << "Something went wrong, the cutvar types were not the same between the data, MC p, and MC n cutvars! Figure it out!" << endl;
	cout << dataVar.getCutVar() << " " << MC_pVar.getCutVar() << " " << MC_nVar.getCutVar() << endl;
	return;
	}
	data_Var = new cutvar(dataVar);
	MC_p_Var = new cutvar(MC_pVar);
	MC_n_Var = new cutvar(MC_nVar);
	
	myCutVar = dataVar.getCutVar();
	
	data_Histo_size = daDataHisto.size();
	mc_p_Histo_size = daMCPHisto.size();
	mc_n_Histo_size = daMCNHisto.size();
	//Basic check that all of the vectors are the same size
	if((data_Histo_size != mc_p_Histo_size) || (data_Histo_size != mc_n_Histo_size) || (mc_p_Histo_size != mc_n_Histo_size)){
	cout << "Histogram vector sizes do not match, figure out what is going on" << data_Histo_size << " , " << mc_p_Histo_size << " , " << mc_n_Histo_size << endl;
	}

	//Initialize the class variables with the input parameters
	
	slice_histo_data = daDataHisto;
	slice_histo_mc_p = daMCPHisto;
	slice_histo_mc_n = daMCNHisto;

	fitType = fit_type;

	//call the stability calculation Rsf function to populate relevant class variable vectors
	stability_calculateRsfQuantities();
	//If applicable scale and shift the MC histograms accordingly. Do it here at the highest level during the constructor
	string mystring("Shift");
	for(int slice = 0; slice < slice_histo_mc_p.size(); slice++){
		//If it is a shift fit scale and shift. If not just scale
		if(fitType.find(mystring) != string::npos){
		slice_histo_mc_p[slice] = scaleAndShiftHisto(slice_histo_mc_p[slice],scale_p_vector[slice],shift_p_vector[slice]);
		slice_histo_mc_n[slice] = scaleAndShiftHisto(slice_histo_mc_n[slice],scale_n_vector[slice],shift_n_vector[slice]);
		}else{
		slice_histo_mc_p[slice]->Scale(scale_p_vector[slice]);
                slice_histo_mc_n[slice]->Scale(scale_n_vector[slice]); 
		}
	}

}//end constructor implemenation

//destructor
//We have dynamically allocated memory, I think. But no pointers so we need to explicitly handle that.
//Everything else should be handled implicitly when it goes out of scope
stability_analysis::~stability_analysis(){
delete data_Var;
delete MC_p_Var;
delete MC_n_Var;
}

//Corresponding getter function implementations

vector<TH1D*> stability_analysis::get_slice_histo_data(){return slice_histo_data;}

vector<TH1D*> stability_analysis::get_slice_histo_mc_p(){return slice_histo_mc_p;}

vector<TH1D*> stability_analysis::get_slice_histo_mc_n(){return slice_histo_mc_n;}

TString stability_analysis::getCutVar(){return myCutVar;}

cutvar* stability_analysis::getDataVar(){return data_Var;}

cutvar* stability_analysis::getMCpVar(){return MC_p_Var;}

cutvar* stability_analysis::getMCnVar(){return MC_n_Var;}

vector<double> stability_analysis::get_scale_p_vector(){return scale_p_vector;}

vector<double> stability_analysis::get_scale_n_vector(){return scale_n_vector;}

vector<double> stability_analysis::get_scale_p_err_vector(){return scale_p_err_vector;}

vector<double> stability_analysis::get_scale_n_err_vector(){return scale_n_err_vector;}

vector<double> stability_analysis::get_shift_p_vector(){return shift_p_vector;}

vector<double> stability_analysis::get_shift_n_vector(){return shift_n_vector;}

vector<double> stability_analysis::get_shift_p_err_vector(){return shift_p_err_vector;}

vector<double> stability_analysis::get_shift_n_err_vector(){return shift_n_err_vector;}

vector<double> stability_analysis::get_ChiSq_vector(){return ChiSq_vector;}

vector<double> stability_analysis::get_ndf_vector(){return ndf_vector;}

vector<double> stability_analysis::get_Rsf_vector(){return Rsf_vector;}

vector<double> stability_analysis::get_Rsf_err_vector(){return Rsf_err_vector;}

vector<vector<double>> stability_analysis::get_poly_result_vector_vectors(){return poly_result_vector_vectors;}

vector<vector<double>> stability_analysis::get_poly_result_err_vector_vectors(){return poly_result_err_vector_vectors;}

vector<TF1*> stability_analysis::get_poly_fit_result_vector(){return poly_fit_result_vector;}

vector<TF1*> stability_analysis::get_total_fit_result_vector(){return total_fit_result_vector;}

//A function that will consider each sliced dx histogram for a given cut variable. And calculate Rsf from all 3 histograms for each slice. Then store relevant information in class vectors. Which will be accessed with getter functions.This is called by the constructor so all of the private variables are properly initialized.
void stability_analysis::stability_calculateRsfQuantities(){
	data_Histo_size = slice_histo_data.size();
	
	if(data_Var->get_dx_hist_low() != MC_p_Var->get_dx_hist_low() || data_Var->get_dx_hist_low() != MC_n_Var->get_dx_hist_low()){
	cout << "The dx hist low for the data, MC p, or MC n does not match. Look into what is going on!" << endl;
	}else if(data_Var->get_dx_hist_high() != MC_p_Var->get_dx_hist_high() || data_Var->get_dx_hist_high() != MC_n_Var->get_dx_hist_high()){
	cout << "The dx hist high for the data, MC p, or MC n does not match. Look into what is going on!" << endl;
	}

	double xmin = data_Var->get_dx_hist_low();
	double xmax = data_Var->get_dx_hist_high();

	//loop over the histos in the vectors, since we require that all of the histo vectors have the same size in the constructor. It should not matter which we use as loop control
	for(int slice_num = 0; slice_num<data_Histo_size; slice_num++){
	//Store the histogram clones for that slice locally
	
	
	TH1D* histo_data = (TH1D*) slice_histo_data[slice_num]->Clone(Form("hist_data_%s_%i",myCutVar.Data(),slice_num));
	TH1D* histo_p = (TH1D*) slice_histo_mc_p[slice_num]->Clone(Form("hist_p_%s_%i",myCutVar.Data(),slice_num));
	TH1D* histo_n = (TH1D*) slice_histo_mc_n[slice_num]->Clone(Form("hist_n_%s_%i",myCutVar.Data(),slice_num));

	//Attempting to make this compatible with the fit histogram class, since that class's job is to handle fitting and getting Rsf and so forth
	fit_histogram *da_slice_histo = new fit_histogram(histo_data,histo_p,histo_n,Form("overall_fit_%s_%i",myCutVar.Data(),slice_num),fitType.c_str(),2,7,xmin,xmax,"QR"); 
	
	//retrieve important fit results
	double scale_p  = da_slice_histo->get_scale_p();
        double scale_p_err = da_slice_histo->get_scale_p_err();
        //Physics result R_sf is a fit parameter
	double Rsf  = da_slice_histo->get_Rsf();
        double Rsf_err = da_slice_histo->get_Rsf_err();

	double scale_n = da_slice_histo->get_scale_n();
	double scale_n_err = da_slice_histo->get_scale_n_err();

        double shift_p= da_slice_histo->get_shift_p();
        double shift_p_err= da_slice_histo->get_shift_p_err();
        double shift_n = da_slice_histo->get_shift_n();
        double shift_n_err = da_slice_histo->get_shift_n_err();

        double ChiSq = da_slice_histo->get_ChiSq();
        double ndf = da_slice_histo->get_NDF();

	vector<double> poly_result = da_slice_histo->get_bkgd_params();
	vector<double> poly_result_err = da_slice_histo->get_bkgd_param_errs();

	//save results into the vectors, which are class variables
	scale_p_vector.push_back(scale_p);
        scale_n_vector.push_back(scale_n);
        scale_p_err_vector.push_back(scale_p_err);
        scale_n_err_vector.push_back(scale_n_err);
        shift_p_vector.push_back(shift_p);
        shift_n_vector.push_back(shift_n);
        shift_p_err_vector.push_back(shift_p_err);
        shift_n_err_vector.push_back(shift_n_err);
        ChiSq_vector.push_back(ChiSq);
        ndf_vector.push_back(ndf);
        Rsf_vector.push_back(Rsf);
        Rsf_err_vector.push_back(Rsf_err);
        poly_result_vector_vectors.push_back(poly_result);
        poly_result_err_vector_vectors.push_back(poly_result_err);

	//Total fit parameters
	vector<pair<double,double>> total_params = da_slice_histo->get_fitParamsErrs();
	double hcalfit_low = da_slice_histo->get_xMin();
	double hcalfit_high = da_slice_histo->get_xMax();

	TF1* total_fit;
	//This conditional works but we will need to modify it if we need any new types of functions
        if(fitType == "fitFull_polyBG"){
        total_fit = new TF1("fit",da_slice_histo,&fit_histogram::fitFull_polyBG,hcalfit_low,hcalfit_high,total_params.size(),"fit_histogram","fitFull_polyBG");
        }else if(fitType =="fitFullShift_polyBG"){
        total_fit = new TF1("fit",da_slice_histo,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,total_params.size(),"fit_histogram","fitFullShift_polyBG");
        }else{
        cout << "The plot function you are trying to implement " << fitType << " is no good! Figure it out now!" << endl;
        total_fit = new TF1("fit",da_slice_histo,&fit_histogram::fitFullShift_polyBG,hcalfit_low,hcalfit_high,total_params.size(),"fit_histogram","fitFullShift_polyBG");
        }
        //file the fit with the parameters
        for(int i=0; i < total_params.size();++i){
        total_fit->SetParameter(i,total_params[i].first);
        total_fit->SetParError(i,total_params[i].second);
        }
	total_fit_result_vector.push_back(total_fit);

	//Create the TF1s for the poly backgrounds here. Use for later functions
	double poly_result_array[poly_result.size()];
	double poly_result_err_array[poly_result.size()];
	

	for(int j = 0; j < poly_result.size(); j++){
	poly_result_array[j] = poly_result[j];
	poly_result_err_array[j] = poly_result_err[j];
	}

	

	TF1 *myFit = new TF1(Form("poly_%i_",slice_num),Form("pol%i",da_slice_histo->get_polyorder()), xmin, xmax);
	myFit->SetParameters( poly_result_array);
	myFit->SetParErrors(poly_result_err_array);
      	myFit->SetNpx(500);
	poly_fit_result_vector.push_back(myFit);

	}//end for loop

}//end Rsf calculate function

//Function to print to terminal some stat info for Rsf
void stability_analysis::printRsfInfo(){

cout<<"---------------------------------------------------------------------------------------------"<<endl;
  double Rsf_mean = utility::calculateMean(Rsf_vector);
  double Rsf_stdev = utility::calculateStDev(Rsf_vector);
  double Rsf_weight_mean = utility::calculateWeightMean(Rsf_vector,Rsf_err_vector);
  double Rsf_weight_stdev = utility::calculateWeightStDev(Rsf_vector,Rsf_err_vector);
  double Rsf_pull = utility::calculatePull(Rsf_vector,Rsf_err_vector);
  cout<<"Rsf mean = " << Rsf_mean<<endl;
  cout<<"Rsf StDev = " << Rsf_stdev<<endl;
  cout<<"weighted Rsf mean = "<<Rsf_weight_mean<<endl;
  cout<<"weighted Rsf StDev = "<<Rsf_weight_stdev<<endl;
  cout<<"Avg. pull on Rsf = " <<Rsf_pull<<endl;
  cout<<"---------------------------------------------------------------------------------------------"<<endl;
  cout<<endl;

}

//Function that plots Rsf for the different slice as a TGraphErrors. Return the canvas
TCanvas* stability_analysis::plotRsfTGraphError(){
  //calculate some stat info for Rsf
  double Rsf_mean = utility::calculateMean(Rsf_vector);
  double Rsf_stdev = utility::calculateStDev(Rsf_vector);
  double Rsf_weight_mean = utility::calculateWeightMean(Rsf_vector,Rsf_err_vector);
  double Rsf_weight_stdev = utility::calculateWeightStDev(Rsf_vector,Rsf_err_vector);
  double Rsf_pull = utility::calculatePull(Rsf_vector,Rsf_err_vector);

  //need info from the cutvar
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();
  //convert the slice ranges to array for TGraph errors, uninitialized arrays
  double x[num_slices];
  double y[num_slices];
  double x_err[num_slices];
  double y_err[num_slices];

  	//Fill the arrays
  	//loop over the number of slices and plot the Rsf value at the central value of the range
  	for(int slice_num = 0; slice_num < num_slices; slice_num++){
	
	//Get the min and max values from the data cutvar vector pair
	double xMin = xMin_xMax_vec[slice_num].first;
	double xMax = xMin_xMax_vec[slice_num].second;
	double xCenter = (xMin + xMax)/2;
	double xWidth = xMax - xMin;

	x[slice_num] = xCenter;
	y[slice_num] = Rsf_vector[slice_num];
        x_err[slice_num] = xWidth/2;
	y_err[slice_num] = Rsf_err_vector[slice_num];

	}

  //Initialize the TGraphErrors
  TGraphErrors *Rsf_graph = new TGraphErrors(num_slices, x, y, x_err, y_err);
  utility::customizeGraph(Rsf_graph, 33, kBlue, 3,"",myCutVar.Data(),"Rsf",1.4,1.4);

  //Fit straight line for the Rsf graph
  double x0 = Rsf_graph ->GetX()[0];
  double xEnd = Rsf_graph->GetX()[Rsf_graph->GetN()-1];
  TF1* fit_Rsf_graph = new TF1("fit_Rsf_graph", "pol1",x0, xEnd);
  fit_Rsf_graph->SetLineColor(kRed);
  Rsf_graph->Fit(fit_Rsf_graph, "Q RN"); // N supresses the drawing of it automatically

  // Fit a pol0 to the graph 
  TF1* fit_pol0_Rsf_graph = new TF1("fit_pol0_Rsf_graph", "pol0",Rsf_graph ->GetX()[0], Rsf_graph->GetX()[Rsf_graph->GetN()-1]);
  fit_pol0_Rsf_graph->SetLineColor(kViolet);
  fit_pol0_Rsf_graph->SetLineWidth(2);
  fit_pol0_Rsf_graph->SetLineStyle(9);
  Rsf_graph->Fit(fit_pol0_Rsf_graph, "Q RN"); // N supresses the drawing of it automatically

  //Make the canvas
  TCanvas *graphCanvas = new TCanvas(Form("Rsf graph stability: %s",myCutVar.Data()),Form("Rsf graph stability: %s",myCutVar.Data()),1600,1200);
  graphCanvas->SetGrid();
  Rsf_graph->Draw("AP");
  fit_Rsf_graph->Draw("same");
  fit_pol0_Rsf_graph->Draw("same");
  graphCanvas->Update();

  //Get fit parameters pol1
  double constant = fit_Rsf_graph->GetParameter(0);
  double constant_Error = fit_Rsf_graph->GetParError(0);
  double slope = fit_Rsf_graph->GetParameter(1);
  double slope_Error = fit_Rsf_graph->GetParError(1);
  double chi2_pol1 = fit_Rsf_graph->GetChisquare();
  int ndf_pol1 = fit_Rsf_graph->GetNDF();
  double chi2_ndf_pol1 = chi2_pol1/ndf_pol1;

  // Get fit parameters pol0
  double constant_pol0 = fit_pol0_Rsf_graph->GetParameter(0);       // The constant value
  double constantError_pol0 = fit_pol0_Rsf_graph->GetParError(0);    // The error on the constan
  double chi2_pol0 = fit_pol0_Rsf_graph->GetChisquare();             // The chi-squared value
  int ndf_pol0= fit_pol0_Rsf_graph->GetNDF();// The number of degrees of freedom
  double chi2_ndf_pol0 = chi2_pol0 / ndf_pol0;

  //Calculate the rise of the linear fit
  double y0 = slope*x0 + constant;
  double yEnd = slope*xEnd +constant;
  double rise = yEnd - y0;

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.025);
  latex.DrawLatex(0.14, 0.87, Form("pol0: y =  %.5f #pm %.5f, Fit #chi^{2}/ndf: %.3f",  constant_pol0,constantError_pol0,chi2_ndf_pol0));
  latex.DrawLatex(0.14, 0.84, Form("pol1: y =  %.5f x + %.5f, Fit #chi^{2}/ndf: %.3f",  slope,constant,chi2_ndf_pol1));
  latex.DrawLatex(0.14, 0.16, Form("Rise across range: %.5f", rise));
  //latex.DrawLatex(0.14, 0.77, Form("Rsf Weighted Mean = %.5f #pm Weighted StDev = %.5f",Rsf_weight_mean, Rsf_weight_stdev));
  //latex.DrawLatex(0.14, 0.73, Form("Rsf pull = %.5f",Rsf_pull) );
  latex.DrawLatex(0.14, 0.13, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  graphCanvas->Update();

  return graphCanvas;
}//end plotRsfTGraphError

//Function that plots each slice dx histogram on a canvas
TCanvas* stability_analysis::plotdxSliceHistos(){

  int nHist = data_Histo_size;
  int nCols = 3;
  int nRows = (nHist + nCols - 1 )/nCols;

  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();

  TCanvas* fits_canvas = new TCanvas(Form("dx slice graph: %s",myCutVar.Data()), Form("dx slice graph: %s",myCutVar.Data()), 1600, 1200);
	if (nHist!=1){
	fits_canvas->Divide(nCols, nRows);
  	}

	for (int i = 0; i < nHist; i++) {
    	fits_canvas->cd(i + 1);
    	fits_canvas->SetGrid();
    	slice_histo_data[i]->SetStats(0);
	slice_histo_data[i]->SetTitle(Form("dx slice graph: %s, Range: %.3f to %.3f",myCutVar.Data(),xMin_xMax_vec[i].first,xMin_xMax_vec[i].second));
	slice_histo_data[i]->Draw();

    	int nEntries = slice_histo_data[i]->GetEntries();
    	double padHeight = fits_canvas->GetWh() / nRows;
   	double legendTextSize = 0.03 * (600.0 / padHeight);  // Adjust based on canvas height

    	// // Create and customize the legend
   	TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    	legend->SetTextSize(legendTextSize);  // Adjust size dynamically
    	legend->SetMargin(0.10);  // Adjust margin to reduce space (default is around 0.25)
    	legend->AddEntry("", Form("R_sf = %.4f #pm %.4f ", Rsf_vector[i],Rsf_err_vector[i]), "");
    	legend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq_vector[i] ,ndf_vector[i]), "");
    	legend->AddEntry("", Form("Entries: %i", nEntries ), "");
    	legend->Draw();

	poly_fit_result_vector[i]->SetLineColor(kGreen);
    	poly_fit_result_vector[i]->Draw("same");
    	//Alter the proton distribution
	slice_histo_mc_p[i]->SetLineColor(kRed);
	slice_histo_mc_p[i]->SetFillColorAlpha(kRed-9,0.5);
	slice_histo_mc_p[i]->SetFillStyle(1001);
	slice_histo_mc_p[i]->SetLineWidth(2);
	slice_histo_mc_p[i]->SetStats(0);
	slice_histo_mc_p[i]->Draw("same");

	//Alter the neutron distribution
    	slice_histo_mc_n[i]->SetLineColor(kBlue);
	slice_histo_mc_n[i]->SetFillColorAlpha(kBlue-9,0.5);
	slice_histo_mc_n[i]->SetFillStyle(1001);
	slice_histo_mc_n[i]->SetLineWidth(2);
	slice_histo_mc_n[i]->SetStats(0);
    	slice_histo_mc_n[i]->Draw("same");
    	
	//Draw the total fit
	total_fit_result_vector[i]->SetLineColor(kMagenta);
    	total_fit_result_vector[i]->Draw("same");

  	}// end loop over slices
  fits_canvas->Update();

  return fits_canvas;
}//end plotdxSliceHistos


//Function that plots each slice residual on a canvas. Residual is between data histogram and total fit
TCanvas* stability_analysis::plotdxSliceResidual(){

  int nHist = data_Histo_size;
  int nCols = 3;
  int nRows = (nHist + nCols - 1 )/nCols;

  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();

  TCanvas* resid_canvas = new TCanvas(Form("Residual slice graph: %s",myCutVar.Data()), Form("Residual slice graph: %s",myCutVar.Data()), 1600, 1200);
        if (nHist!=1){
        resid_canvas->Divide(nCols, nRows);
        }

        for (int i = 0; i < nHist; i++) {
        resid_canvas->cd(i + 1);
        resid_canvas->SetGrid();
	TH1D* myResidual = makeResidualError(slice_histo_data[i],total_fit_result_vector[i]);
	
	myResidual->Draw("E sames");
	}

	resid_canvas->Update();

return resid_canvas;
}

//Function that scales and shifts, if applicable the MC histograms.
TH1D* stability_analysis::scaleAndShiftHisto(TH1D* origHist, double scale, double shift){

  //create a histogram that is indentical to the original
  TH1D* scaleShiftHist = (TH1D*) origHist->Clone("scaleShiftHist");
  //loop over the bins and apply the shifting and scaling
  for(int i=1; i<=origHist->GetNbinsX(); ++i){
  //get information about the bin in the orig hist
  double origBinCenter = origHist->GetBinCenter(i);

  double shiftedBinCenter = origBinCenter - shift;
  double value = origHist->Interpolate(shiftedBinCenter);

  //Check if the shifted Histogrom bin is within the histogram range
  	if((shiftedBinCenter >= origHist->GetXaxis()->GetXmin()) && (shiftedBinCenter <= origHist->GetXaxis()->GetXmax()) ){
	scaleShiftHist->SetBinContent(i,scale * value);

	//propogate the error
	double error = origHist->GetBinError(origHist->FindBin(shiftedBinCenter));
	scaleShiftHist->SetBinError(i,scale*error);

  	}else{
	//if outside the range set bin content and error to zero
	scaleShiftHist->SetBinContent(i,0);
	scaleShiftHist->SetBinError(i,0);

	}
  }

return scaleShiftHist;
}

//Function that makes a residual between the data histogram and the total fit
TH1D* stability_analysis::makeResidualError(TH1D* hist, TF1* fit){
int hist_Nbins = hist->GetNbinsX();
TString hist_name = hist->GetName();
double hist_minX = hist->GetXaxis()->GetXmin();
double hist_maxX = hist->GetXaxis()->GetXmax();


TH1D* resid_hist = new TH1D(Form("Residual Histogram: %s and total fit",hist_name.Data()),Form("Residual Histogram: %s and total fit",hist_name.Data()),hist_Nbins,hist_minX,hist_maxX);

	for(int bin = 1; bin <= hist_Nbins; bin++){
	
	double hist_val = hist->GetBinContent(bin);
        double hist_error = hist->GetBinError(bin);
	double fit_val = fit->Eval(hist->GetXaxis()->GetBinCenter(bin));

	double new_val = hist_val - fit_val;
        resid_hist->SetBinContent(bin,new_val);
        resid_hist->SetBinError(bin,hist_error);
	}
return resid_hist;
}

//Function that shows the slice region overlayed on the 2D plot
TCanvas* stability_analysis::plot2DCutOverlay(){

  //need info from the cutvar
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();

  TCanvas* cut_2Dcanvas = new TCanvas(Form("Cut Visualization 2D: %s",myCutVar.Data()),Form("Cut Visualization 2D: %s",myCutVar.Data()),1600,1200);
  cut_2Dcanvas->Divide(3,1);
  int color_1 = kRed;

  //Data histogram
  cut_2Dcanvas->cd(1);
  //Draw the data histogram 
  TH2D* hist_2D_data = (TH2D*)data_Var->get2DdxCutHisto()->Clone();;
  hist_2D_data->SetTitle(Form("Data Cut Visualization 2D: dx vs. %s",myCutVar.Data()));
  hist_2D_data->SetStats(0);
  hist_2D_data->Draw("COLZ");
  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_2D_data->GetYaxis()->GetXmin(), xLeft, hist_2D_data->GetYaxis()->GetXmax());
    TLine* lineRight = new TLine(xRight, hist_2D_data->GetYaxis()->GetXmin(),xRight, hist_2D_data->GetYaxis()->GetXmax());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  //MC p histogram
  cut_2Dcanvas->cd(2);
  //Draw the MC p histogram 
  TH2D* hist_2D_mc_p = (TH2D*)MC_p_Var->get2DdxCutHisto()->Clone();
  hist_2D_mc_p->SetTitle(Form("MC P Cut Visualization 2D: dx vs. %s",myCutVar.Data()));
  hist_2D_mc_p->SetStats(0);
  hist_2D_mc_p->Draw("COLZ");
  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_2D_mc_p->GetYaxis()->GetXmin(), xLeft, hist_2D_mc_p->GetYaxis()->GetXmax());
    TLine* lineRight = new TLine(xRight, hist_2D_mc_p->GetYaxis()->GetXmin(),xRight, hist_2D_mc_p->GetYaxis()->GetXmax());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  //MC n histogram
  cut_2Dcanvas->cd(3);
  //Draw the MC p histogram 
  TH2D* hist_2D_mc_n = (TH2D*)MC_n_Var->get2DdxCutHisto()->Clone();
  hist_2D_mc_n->SetTitle(Form("MC N Cut Visualization 2D: dx vs. %s",myCutVar.Data()));
  hist_2D_mc_n->SetStats(0);
  hist_2D_mc_n->Draw("COLZ");
  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_2D_mc_n->GetYaxis()->GetXmin(), xLeft, hist_2D_mc_n->GetYaxis()->GetXmax());
    TLine* lineRight = new TLine(xRight, hist_2D_mc_n->GetYaxis()->GetXmin(),xRight, hist_2D_mc_n->GetYaxis()->GetXmax());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

cut_2Dcanvas->Update();
return cut_2Dcanvas;
}

TCanvas* stability_analysis::plot2DCutOverlay_CutRegion(){

  //need info from the cutvar
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();

  TCanvas* cut_2Dcanvas = new TCanvas(Form("Cut Visualization 2D in Region: %s",myCutVar.Data()),Form("Cut Visualization 2D in Region: %s",myCutVar.Data()),1600,1200);
  cut_2Dcanvas->Divide(3,1);
  double xmin_2D = xMin_xMax_vec.front().first;
  double xmax_2D = xMin_xMax_vec.back().second;
  int color_1 = kRed;

  //Data histogram
  cut_2Dcanvas->cd(1);
  //Draw the data histogram
  TH2D* hist_2D_data = (TH2D*)data_Var->get2DdxCutHisto()->Clone();
  hist_2D_data->GetXaxis()->SetRangeUser(xmin_2D,xmax_2D);
  hist_2D_data->SetTitle(Form("Data Cut Visualization 2D: dx vs. %s",myCutVar.Data()));
  hist_2D_data->SetStats(0);
  hist_2D_data->Draw("COLZ");
  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_2D_data->GetYaxis()->GetXmin(), xLeft, hist_2D_data->GetYaxis()->GetXmax());
    TLine* lineRight = new TLine(xRight, hist_2D_data->GetYaxis()->GetXmin(),xRight, hist_2D_data->GetYaxis()->GetXmax());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  //MC p histogram
  cut_2Dcanvas->cd(2);
  //Draw the MC p histogram
  TH2D* hist_2D_mc_p = (TH2D*)MC_p_Var->get2DdxCutHisto()->Clone();

  hist_2D_mc_p->GetXaxis()->SetRangeUser(xmin_2D,xmax_2D);
  hist_2D_mc_p->SetTitle(Form("MC P Cut Visualization 2D: dx vs. %s",myCutVar.Data()));
  hist_2D_mc_p->SetStats(0);
  hist_2D_mc_p->Draw("COLZ");
  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_2D_mc_p->GetYaxis()->GetXmin(), xLeft, hist_2D_mc_p->GetYaxis()->GetXmax());
    TLine* lineRight = new TLine(xRight, hist_2D_mc_p->GetYaxis()->GetXmin(),xRight, hist_2D_mc_p->GetYaxis()->GetXmax());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  //MC n histogram
  cut_2Dcanvas->cd(3);
  //Draw the MC p histogram
  TH2D* hist_2D_mc_n = (TH2D*)MC_n_Var->get2DdxCutHisto()->Clone();
  hist_2D_mc_n->GetXaxis()->SetRangeUser(xmin_2D,xmax_2D);
  hist_2D_mc_n->SetTitle(Form("MC N Cut Visualization 2D: dx vs. %s",myCutVar.Data()));
  hist_2D_mc_n->SetStats(0);
  hist_2D_mc_n->Draw("COLZ");
  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_2D_mc_n->GetYaxis()->GetXmin(), xLeft, hist_2D_mc_n->GetYaxis()->GetXmax());
    TLine* lineRight = new TLine(xRight, hist_2D_mc_n->GetYaxis()->GetXmin(),xRight, hist_2D_mc_n->GetYaxis()->GetXmax());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
   }
    cut_2Dcanvas->Update();
return cut_2Dcanvas;
}

//Function that shows the slice region overlayed on the 1D plot
TCanvas* stability_analysis::plot1DCutOverlay(){
  
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();
  int color_1 = kRed;
  TCanvas* cut_1Dcanvas = new TCanvas(Form("Cut Visualization 1D: %s",myCutVar.Data()),Form("Cut Visualization 1D: %s",myCutVar.Data()),1600,1200);
  cut_1Dcanvas->Divide(1,3);
  
  //Data
  cut_1Dcanvas->cd(1);
  TH1D* hist_1D_data = (TH1D*) (data_Var->get2DdxCutHisto())->ProjectionX()->Clone();
  int daDataEntries = hist_1D_data->GetEntries();
  hist_1D_data->SetTitle(Form("Data Cut Visualization 1D: %s",myCutVar.Data()));
  hist_1D_data->SetStats(0);
  hist_1D_data->Draw();

  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_1D_data->GetMinimum(), xLeft, hist_1D_data->GetMaximum());
    TLine* lineRight = new TLine(xRight, hist_1D_data->GetMinimum(),xRight, hist_1D_data->GetMaximum());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.065);
  latex.DrawLatex(0.7, 0.7, Form("Number Entries: %i",  daDataEntries));


  //MC P
  cut_1Dcanvas->cd(2);
  TH1D* hist_1D_mc_p = (TH1D*) (MC_p_Var->get2DdxCutHisto())->ProjectionX()->Clone();
  hist_1D_mc_p->SetTitle(Form("MC P Cut Visualization 1D: %s",myCutVar.Data()));
  hist_1D_mc_p->SetStats(0);
  hist_1D_mc_p->Draw();

  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_1D_mc_p->GetMinimum(), xLeft, hist_1D_mc_p->GetMaximum());
    TLine* lineRight = new TLine(xRight, hist_1D_mc_p->GetMinimum(),xRight, hist_1D_mc_p->GetMaximum());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  //MC N
  cut_1Dcanvas->cd(3);
  TH1D* hist_1D_mc_n = (TH1D*) (MC_n_Var->get2DdxCutHisto())->ProjectionX()->Clone();
  hist_1D_mc_n->SetTitle(Form("MC N Cut Visualization 1D: %s",myCutVar.Data()));
  hist_1D_mc_n->SetStats(0);
  hist_1D_mc_n->Draw();

  // Draw vertical lines at each x value
  for (int slice = 0; slice < num_slices; slice++) {
    double xLeft = xMin_xMax_vec[slice].first;
    double xRight =  xMin_xMax_vec[slice].second;
    TLine* lineLeft = new TLine(xLeft, hist_1D_mc_n->GetMinimum(), xLeft, hist_1D_mc_n->GetMaximum());
    TLine* lineRight = new TLine(xRight, hist_1D_mc_n->GetMinimum(),xRight, hist_1D_mc_n->GetMaximum());
    lineLeft->SetLineColor(color_1);  // Set line color
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color_1);  // Set line color
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  cut_1Dcanvas->Update();

  return cut_1Dcanvas;
}

//Function that shows the number of entries in each slice
TCanvas* stability_analysis::plotNEntries(){
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();

  TCanvas* nEntries_canvas = new TCanvas(Form("N Entries: %s",myCutVar.Data()),Form("N Entries: %s",myCutVar.Data()),1600,1200);
  //convert the slice ranges to array for TGraph errors, uninitialized arrays
  double x_n[num_slices];
  double y_n[num_slices];
  double x_n_err[num_slices];
  double y_n_err[num_slices];

  //Fill the arrays
        //loop over the number of slices and plot the number of entries at the central value of the range
        for(int slice_num = 0; slice_num < num_slices; slice_num++){
	double xMin = xMin_xMax_vec[slice_num].first;
        double xMax = xMin_xMax_vec[slice_num].second;
        double xCenter = (xMin + xMax)/2;
        double xWidth = xMax - xMin;

	x_n[slice_num] = xCenter;
        y_n[slice_num] = slice_histo_data[slice_num]->GetEntries();
        x_n_err[slice_num] = xWidth/2;
        y_n_err[slice_num] = 0;
	}
  //Initialize the TGraphErrors
  TGraphErrors *nEntries_graph = new TGraphErrors(num_slices, x_n, y_n, x_n_err, y_n_err);
  utility::customizeGraph(nEntries_graph, 33, kBlue, 3,"",myCutVar.Data(),"nEntries",1.4,1.4);

  nEntries_graph->Draw("AP");

  nEntries_canvas->Update();
return nEntries_canvas;
}

//Function that plots Rsf for the different slice as a TGraphErrors at the XMin values. Return the canvas
TCanvas* stability_analysis::plotRsfTGraphError_xMin(){
//calculate some stat info for Rsf
  double Rsf_mean = utility::calculateMean(Rsf_vector);
  double Rsf_stdev = utility::calculateStDev(Rsf_vector);
  double Rsf_weight_mean = utility::calculateWeightMean(Rsf_vector,Rsf_err_vector);
  double Rsf_weight_stdev = utility::calculateWeightStDev(Rsf_vector,Rsf_err_vector);
  double Rsf_pull = utility::calculatePull(Rsf_vector,Rsf_err_vector);

  //need info from the cutvar
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();
  //convert the slice ranges to array for TGraph errors, uninitialized arrays
  double x_min[num_slices];
  double y_min[num_slices];
  double x_min_err[num_slices];
  double y_min_err[num_slices];

        //Fill the arrays
        //loop over the number of slices and plot the Rsf value at the central value of the range
        for(int slice_num = 0; slice_num < num_slices; slice_num++){

        //Get the min and max values from the data cutvar vector pair
        double xMin = xMin_xMax_vec[slice_num].first;

        x_min[slice_num] = xMin;
        y_min[slice_num] = Rsf_vector[slice_num];
        x_min_err[slice_num] = 0;
        y_min_err[slice_num] = Rsf_err_vector[slice_num];

        }
 
  //Initialize the TGraphErrors
  TGraphErrors *Rsf_xMin_graph = new TGraphErrors(num_slices, x_min, y_min, x_min_err, y_min_err);
  utility::customizeGraph(Rsf_xMin_graph, 33, kBlue, 3,"",myCutVar.Data(),"Rsf",1.4,1.4);

  //Fit straight line for the Rsf graph
  double x0_xMin = Rsf_xMin_graph ->GetX()[0];
  double xEnd_xMin = Rsf_xMin_graph->GetX()[Rsf_xMin_graph->GetN()-1];
  TF1* fit_Rsf_xMin_graph = new TF1("fit_Rsf_xMin_graph", "pol1",x0_xMin, xEnd_xMin);
  fit_Rsf_xMin_graph->SetLineColor(kRed);
  Rsf_xMin_graph->Fit(fit_Rsf_xMin_graph, "Q RN"); // N supresses the drawing of it automatically

  // Fit a pol0 to the graph 
  TF1* fit_pol0_Rsf_xMin_graph = new TF1("fit_pol0_Rsf_xMin_graph", "pol0",Rsf_xMin_graph ->GetX()[0], Rsf_xMin_graph->GetX()[Rsf_xMin_graph->GetN()-1]);
  fit_pol0_Rsf_xMin_graph->SetLineColor(kViolet);
  fit_pol0_Rsf_xMin_graph->SetLineWidth(2);
  fit_pol0_Rsf_xMin_graph->SetLineStyle(9);
  Rsf_xMin_graph->Fit(fit_pol0_Rsf_xMin_graph, "Q RN"); // N supresses the drawing of it automatically

  //Make the canvas
  TCanvas *RsfxMin_Canvas = new TCanvas(Form("Rsf xMin graph stability: %s",myCutVar.Data()),Form("Rsf xMin graph stability: %s",myCutVar.Data()),1600,1200);
  RsfxMin_Canvas->SetGrid();
  Rsf_xMin_graph->Draw("AP");
  fit_Rsf_xMin_graph->Draw("same");
  fit_pol0_Rsf_xMin_graph->Draw("same");
  RsfxMin_Canvas->Update();

  //Get fit parameters pol1
  double constant_xMin = fit_Rsf_xMin_graph->GetParameter(0);
  double constant_Error_xMin = fit_Rsf_xMin_graph->GetParError(0);
  double slope_xMin = fit_Rsf_xMin_graph->GetParameter(1);
  double slope_Error_xMin = fit_Rsf_xMin_graph->GetParError(1);
  double chi2_pol1_xMin = fit_Rsf_xMin_graph->GetChisquare();
  int ndf_pol1_xMin = fit_Rsf_xMin_graph->GetNDF();
  double chi2_ndf_pol1_xMin = chi2_pol1_xMin/ndf_pol1_xMin;

  // Get fit parameters pol0
  double constant_pol0_xMin = fit_pol0_Rsf_xMin_graph->GetParameter(0);       // The constant value
  double constantError_pol0_xMin = fit_pol0_Rsf_xMin_graph->GetParError(0);    // The error on the constan
  double chi2_pol0_xMin = fit_pol0_Rsf_xMin_graph->GetChisquare();             // The chi-squared value
  int ndf_pol0_xMin = fit_pol0_Rsf_xMin_graph->GetNDF();// The number of degrees of freedom
  double chi2_ndf_pol0_xMin = chi2_pol0_xMin / ndf_pol0_xMin;

  //Calculate the rise of the linear fit
  double y0_xMin = slope_xMin*x0_xMin + constant_xMin;
  double yEnd_xMin = slope_xMin*xEnd_xMin +constant_xMin;
  double rise_xMin = yEnd_xMin - y0_xMin;  
 
  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.025);
  latex.DrawLatex(0.14, 0.87, Form("pol0: y =  %.5f #pm %.5f, Fit #chi^{2}/ndf: %.3f",  constant_pol0_xMin,constantError_pol0_xMin,chi2_ndf_pol0_xMin));
  latex.DrawLatex(0.14, 0.84, Form("pol1: y =  %.5f x + %.5f, Fit #chi^{2}/ndf: %.3f",  slope_xMin,constant_xMin,chi2_ndf_pol1_xMin));
  latex.DrawLatex(0.14, 0.16, Form("Rise across range: %.5f", rise_xMin));
  latex.DrawLatex(0.14, 0.13, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  RsfxMin_Canvas->Update();
return RsfxMin_Canvas;
}// end plot Rsf xMin Function

//Function that plots Rsf for the different slice as a TGraphErrors at the XMax values. Return the canvas
TCanvas* stability_analysis::plotRsfTGraphError_xMax(){
//calculate some stat info for Rsf
  double Rsf_mean = utility::calculateMean(Rsf_vector);
  double Rsf_stdev = utility::calculateStDev(Rsf_vector);
  double Rsf_weight_mean = utility::calculateWeightMean(Rsf_vector,Rsf_err_vector);
  double Rsf_weight_stdev = utility::calculateWeightStDev(Rsf_vector,Rsf_err_vector);
  double Rsf_pull = utility::calculatePull(Rsf_vector,Rsf_err_vector);

  //need info from the cutvar
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();
  //convert the slice ranges to array for TGraph errors, uninitialized arrays
  double x_max[num_slices];
  double y_max[num_slices];
  double x_max_err[num_slices];
  double y_max_err[num_slices];

        //Fill the arrays
        //loop over the number of slices and plot the Rsf value at the central value of the range
        for(int slice_num = 0; slice_num < num_slices; slice_num++){

        //Get the min and max values from the data cutvar vector pair
        double xMax = xMin_xMax_vec[slice_num].second;

        x_max[slice_num] = xMax;
        y_max[slice_num] = Rsf_vector[slice_num];
        x_max_err[slice_num] = 0;
        y_max_err[slice_num] = Rsf_err_vector[slice_num];

        }

  //Initialize the TGraphErrors
  TGraphErrors *Rsf_xMax_graph = new TGraphErrors(num_slices, x_max, y_max, x_max_err, y_max_err);
  utility::customizeGraph(Rsf_xMax_graph, 33, kBlue, 3,"",myCutVar.Data(),"Rsf",1.4,1.4);

  //Fit straight line for the Rsf graph
  double x0_xMax = Rsf_xMax_graph ->GetX()[0];
  double xEnd_xMax = Rsf_xMax_graph->GetX()[Rsf_xMax_graph->GetN()-1];
  TF1* fit_Rsf_xMax_graph = new TF1("fit_Rsf_xMax_graph", "pol1",x0_xMax, xEnd_xMax);
  fit_Rsf_xMax_graph->SetLineColor(kRed);
  Rsf_xMax_graph->Fit(fit_Rsf_xMax_graph, "Q RN"); // N supresses the drawing of it automatically

  // Fit a pol0 to the graph 
  TF1* fit_pol0_Rsf_xMax_graph = new TF1("fit_pol0_Rsf_xMax_graph", "pol0",Rsf_xMax_graph ->GetX()[0], Rsf_xMax_graph->GetX()[Rsf_xMax_graph->GetN()-1]);
  fit_pol0_Rsf_xMax_graph->SetLineColor(kViolet);
  fit_pol0_Rsf_xMax_graph->SetLineWidth(2);
  fit_pol0_Rsf_xMax_graph->SetLineStyle(9);
  Rsf_xMax_graph->Fit(fit_pol0_Rsf_xMax_graph, "Q RN"); // N supresses the drawing of it automatically

  //Make the canvas
  TCanvas *RsfxMax_Canvas = new TCanvas(Form("Rsf xMax graph stability: %s",myCutVar.Data()),Form("Rsf xMax graph stability: %s",myCutVar.Data()),1600,1200);
  RsfxMax_Canvas->SetGrid();
  Rsf_xMax_graph->Draw("AP");
  fit_Rsf_xMax_graph->Draw("same");
  fit_pol0_Rsf_xMax_graph->Draw("same");
  RsfxMax_Canvas->Update();

  //Get fit parameters pol1
  double constant_xMax = fit_Rsf_xMax_graph->GetParameter(0);
  double constant_Error_xMax = fit_Rsf_xMax_graph->GetParError(0);
  double slope_xMax = fit_Rsf_xMax_graph->GetParameter(1);
  double slope_Error_xMax = fit_Rsf_xMax_graph->GetParError(1);
  double chi2_pol1_xMax = fit_Rsf_xMax_graph->GetChisquare();
  int ndf_pol1_xMax = fit_Rsf_xMax_graph->GetNDF();
  double chi2_ndf_pol1_xMax = chi2_pol1_xMax/ndf_pol1_xMax;

  // Get fit parameters pol0
  double constant_pol0_xMax = fit_pol0_Rsf_xMax_graph->GetParameter(0);       // The constant value
  double constantError_pol0_xMax = fit_pol0_Rsf_xMax_graph->GetParError(0);    // The error on the constan
  double chi2_pol0_xMax = fit_pol0_Rsf_xMax_graph->GetChisquare();             // The chi-squared value
  int ndf_pol0_xMax = fit_pol0_Rsf_xMax_graph->GetNDF();// The number of degrees of freedom
  double chi2_ndf_pol0_xMax = chi2_pol0_xMax / ndf_pol0_xMax;
  //Calculate the rise of the linear fit
  double y0_xMax = slope_xMax*x0_xMax + constant_xMax;
  double yEnd_xMax = slope_xMax*xEnd_xMax +constant_xMax;
  double rise_xMax = yEnd_xMax - y0_xMax;

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.025);
  latex.DrawLatex(0.14, 0.87, Form("pol0: y =  %.5f #pm %.5f, Fit #chi^{2}/ndf: %.3f",  constant_pol0_xMax,constantError_pol0_xMax,chi2_ndf_pol0_xMax));
  latex.DrawLatex(0.14, 0.84, Form("pol1: y =  %.5f x + %.5f, Fit #chi^{2}/ndf: %.3f",  slope_xMax,constant_xMax,chi2_ndf_pol1_xMax));
  latex.DrawLatex(0.14, 0.16, Form("Rise across range: %.5f", rise_xMax));
  latex.DrawLatex(0.14, 0.13, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  RsfxMax_Canvas->Update();
return RsfxMax_Canvas;
}// end plot Rsf xMax Function

//Function that shows the Chi2/ndf graph
TCanvas* stability_analysis::plotChi2_NDFGraph(){

//need info from the cutvar
  vector<pair<double,double>> xMin_xMax_vec = data_Var->getXMinXMaxRange();
  int num_slices = xMin_xMax_vec.size();

//// Plot Chi2/ndf
//// Make arrays that TGraphErrors can use

  double x_ch[num_slices];
  double y_ch[num_slices];
  double x_ch_err[num_slices];
  double y_ch_err[num_slices];

        for (int slice_num = 0 ; slice_num < num_slices ; slice_num++){
        //Get the min and max values from the data cutvar vector pair
        double xMin = xMin_xMax_vec[slice_num].first;
        double xMax = xMin_xMax_vec[slice_num].second;
        double xCenter = (xMin + xMax)/2;
        double xWidth = xMax - xMin;

        x_ch[slice_num] = xCenter;
        y_ch[slice_num] = ChiSq_vector[slice_num] /ndf_vector[slice_num];
        x_ch_err[slice_num] = xWidth;
        y_ch_err[slice_num] = 0;
        }

  TGraphErrors *Chi2_ndf_graph = new TGraphErrors(num_slices, x_ch, y_ch, x_ch_err, y_ch_err);
  utility::customizeGraph(Chi2_ndf_graph, 33, kBlue, 3,"","Bin Width","chi^{2}/ndf",1.4,1.4);

  //Make the canvas
  TCanvas *graphCanvas = new TCanvas(Form("Chi^{2}/ndf stability: %s",myCutVar.Data()),Form("Chi^{2}/ndf stability: %s",myCutVar.Data()),1600,1200);
  graphCanvas->SetGrid();
  Chi2_ndf_graph ->Draw("AP");
  graphCanvas->Update();
  return graphCanvas;
}
