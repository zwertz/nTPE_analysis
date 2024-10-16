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
	
	}//end for loop
}//end Rsf calculate function

