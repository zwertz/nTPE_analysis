#ifndef STABILITY_ANALYSIS_H
#define STABILITY_ANALYSIS_H

//Author: Ezekiel Wertz
//This is a companion class to cutvar.h and cutvar.C. I ran into a challenge where I needed to be able to manipulate the information for all 3 histograms for data, mc p, and mc n. Building this into cutvar seemed difficult. So instead this class will take output, after some processing, from cutvar. Then do further cut stability study manipulation/evaluation here. This class should handle all information from data and MC at once. Will be modified as stability and systematic studies progress.

class stability_analysis{
private:
vector<TH1D*> slice_histo_data, slice_histo_mc_p, slice_histo_mc_n;
cutvar myCutVar;
TH1D* hist_data, hist_p, hist_n;
vector<double> scale_p_vector, scale_n_vector, scale_p_err_vector, scale_n_err_vector, shift_p_vector, shift_n_vector, shift_p_err_vector, shift_n_err_vector, ChiSq_vector, ndf_vector, Rsf_vector, Rsf_err_vector;
vector<vector<double>> poly_result_vector_vectors, poly_result_err_vector_vectors;

//A function that will consider each sliced dx histogram for a given cut variable. And calculate Rsf from all 3 histograms for each slice. Then store relevant information in class vectors. Which will be accessed with getter functions. This is called by the constructor so all of the private variables are properly initialized.
void stability_calculateRsfQuantities();

public:
//Standard constructor, it should initialize all the class variables we care about.
stability_analysis(cutvar datVar,vector<TH1D*> daDataHisto, vector<TH1D*> daMCPHisto, vector<TH1D*> daMCNHisto);

//getter functions for various class variables

vector<TH1D*> get_slice_histo_data();

vector<TH1D*> get_slice_histo_mc_p();

vector<TH1D*> get_slice_histo_mc_n();

TString getmyCutVar();

vector<double> get_scale_p_vector();

vector<double> get_scale_n_vector();

vector<double> get_scale_p_err_vector();

vector<double> get_scale_n_err_vector();

vector<double> get_shift_p_vector();

vector<double> get_shift_n_vector();

vector<double> get_shift_p_err_vector();

vector<double> get_shift_n_err_vector();

vector<double> get_ChiSq_vector();

vector<double> get_ndf_vector();

vector<double> get_Rsf_vector();

vector<double> get_Rsf_err_vector();

vector<vector<double>> get_poly_result_vector_vectors(); 

vector<vector<double>> get_poly_result_err_vector_vectors();


};//end class
#endif
