//Author: Ezekiel Wertz
//This is a companion class to cutvar.h and cutvar.C. I ran into a challenge where I needed to be able to manipulate the information for all 3 histograms for data, mc p, and mc n.Building this into cutvar seemed difficult. So instead this class will take output, after some processing, from cutvar. Then do further cut stability study manipulation/evaluation here. This class should handle all information from data and MC at once. Will be modified as stability and systematic studies progress.

#include "../include/stability_analysis.h"

stability_analysis::stability_analysis(cutvar &dataVar,cutvar &MC_pVar,cutvar &MC_nVar, vector<TH1D*> daDataHisto, vector<TH1D*> daMCPHisto, vector<TH1D*> daMCNHisto,const string& fitFormula, TH1D* histo_p, TH1D* histo_n){

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

	//call the stability calculation Rsf function to populate relevant class variable vectors
	stability_calculateRsfQuantities(fitFormula, histo_p,histo_n);
	

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

//A function that will consider each sliced dx histogram for a given cut variable. And calculate Rsf from all 3 histograms for each slice. Then store relevant information in class vectors. Which will be accessed with getter functions.This is called by the constructor so all of the private variables are properly initialized.
void stability_analysis::stability_calculateRsfQuantities(const string& fitFormula,TH1D* histo_p, TH1D* histo_n){
	data_Histo_size = slice_histo_data.size();
	
	if(data_Var->get_dx_hist_low() != MC_p_Var->get_dx_hist_low() || data_Var->get_dx_hist_low() != MC_n_Var->get_dx_hist_low()){
	cout << "The dx hist low for the data, MC p, or MC n does not match. Look into what is going on!" << endl;
	}else if(data_Var->get_dx_hist_high() != MC_p_Var->get_dx_hist_high() || data_Var->get_dx_hist_high() != MC_n_Var->get_dx_hist_high()){
	cout << "The dx hist high for the data, MC p, or MC n does not match. Look into what is going on!" << endl;
	}

	double xmin = data_Var->get_dx_hist_low();
	double xmax = data_Var->get_dx_hist_high();

	//initial parameter guess
	double initialParameters[7]={1,1,-0.05,-0.05,1,1,-1};
	//loop over the histos in the vectors, since we require that all of the histo vectors have the same size in the constructor. It should not matter which we use as loop control
	for(int slice_num = 0; slice_num<data_Histo_size; slice_num++){
	//Store the histogram clones for that slice locally
	
	
	TH1D* histo_data = (TH1D*) slice_histo_data[slice_num]->Clone(Form("hist_data_%s_%i",myCutVar.Data(),slice_num));
	histo_p = (TH1D*) slice_histo_mc_p[slice_num]->Clone(Form("hist_p_%s_%i",myCutVar.Data(),slice_num));
	histo_n = (TH1D*) slice_histo_mc_n[slice_num]->Clone(Form("hist_n_%s_%i",myCutVar.Data(),slice_num));

	//Make a fit to the overall histogram for dx. In this case we used the MC proton and neutron shapes + a 2nd order poly background
	TF1 *myFit = new TF1(Form("overall_fit_%s_%i",myCutVar.Data(),slice_num),fitFormula.c_str(),xmin,xmax,7);

	myFit->SetParameters(initialParameters);
        // set parameter limits. SetParLimits -> (par#, min, max)
        myFit->SetParLimits(0, 0, 2000); // scale_p greater than 0
        myFit->SetParLimits(1, 0,2000); // scale_n greater than 0
        myFit->SetParLimits(2, -0.10,0.10); // shift_p less than +- 10cm
        myFit->SetParLimits(3, -0.10,0.10); // shift_n less than +- 10cm
        myFit->SetParLimits(6,-10000000,-0.0000000001); // x^2 term negative to force downward concavity 
        myFit->SetNpx(500);
	histo_data->GetXaxis()->SetRangeUser(xmin,xmax);

	//now fit the data histogram with the fit info
	histo_data->Fit(myFit,"Q R");

	//retrieve important fit results
	double scale_p  = myFit->GetParameter(0);
        double scale_p_err = myFit->GetParError(0);
        //Physics result R_sf is a fit parameter
	double Rsf  = myFit->GetParameter(1);
        double Rsf_err = myFit->GetParError(1);

	double scale_n = Rsf*scale_p;
	double scale_n_err = scale_n*sqrt(pow(Rsf_err/Rsf,2)+pow(scale_p_err/scale_p,2));

        double shift_p= myFit->GetParameter(2);
        double shift_p_err= myFit->GetParError(2);
        double shift_n = myFit->GetParameter(3);
        double shift_n_err = myFit->GetParError(3);

        double ChiSq= myFit->GetChisquare();
        double ndf = myFit->GetNDF();

	vector<double> poly_result;
	vector<double> poly_result_err;

		for (int i =0 ; i < 3; i++){
          	poly_result.push_back( myFit->GetParameter(4+i) );
          	poly_result_err.push_back(myFit->GetParError(4+i) );
        	}

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
	
	delete myFit;
	}//end for loop
}//end Rsf calculate function

