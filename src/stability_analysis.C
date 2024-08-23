//Author: Ezekiel Wertz
//This is a companion class to cutvar.h and cutvar.C. I ran into a challenge where I needed to be able to manipulate the information for all 3 histograms for data, mc p, and mc n.Building this into cutvar seemed difficult. So instead this class will take output, after some processing, from cutvar. Then do further cut stability study manipulation/evaluation here. This class should handle all information from data and MC at once. Will be modified as stability and systematic studies progress.

stability_analysis::stability_analysis(cutvar datVar,vector<TH1D*> daDataHisto, vector<TH1D*> daMCPHisto, vector<TH1D*> daMCNHisto){

	cutvar myCutvar =datVar;
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
	stability_calculateRsfQuantities();
	

}//end constructor implemenation

//Corresponding getter function implementations

vector<TH1D*> get_slice_histo_data(){return slice_histo_data;}

vector<TH1D*> get_slice_histo_mc_p(){return slice_histo_mc_p;}

vector<TH1D*> get_slice_histo_mc_n(){return slice_histo_mc_n;}

TString getmyCutVar(){return myCutVar;}

vector<double> get_scale_p_vector(){return scale_p_vector;}

vector<double> get_scale_n_vector(){return scale_n_vector;}

vector<double> get_scale_p_err_vector(){return scale_p_err_vector;}

vector<double> get_scale_n_err_vector(){return scale_n_err_vector;}

vector<double> get_shift_p_vector(){return shift_p_vector;}

vector<double> get_shift_n_vector(){return shift_n_vector;}

vector<double> get_shift_p_err_vector(){return shift_p_err_vector;}

vector<double> get_shift_n_err_vector(){return shift_n_err_vector;}

vector<double> get_ChiSq_vector(){return ChiSq_vector;}

vector<double> get_ndf_vector();

vector<double> get_Rsf_vector();

vector<double> get_Rsf_err_vector();

vector<vector<double>> get_poly_result_vector_vectors();

vector<vector<double>> get_poly_result_err_vector_vectors();

//A function that will consider each sliced dx histogram for a given cut variable. And calculate Rsf from all 3 histograms for each slice. Then store relevant information in class vectors. Which will be accessed with getter functions.This is called by the constructor so all of the private variables are properly initialized.
void stability_analysis::stability_calculateRsfQuantities(){
	data_Histo_size = slice_histo_data.size();
	double xmin = myCutVar.get_dx_hist_low();
	double xmax = myCutVar.get_dx_hist_high();

	//initial parameter guess
	double initialParameters[7]={1,1,-0.05,-0.05,1,1,-1};
	//loop over the histos in the vectors, since we require that all of the histo vectors have the same size in the constructor. It should not matter which we use as loop control
	for(int slice_num = 0; slice_num<data_Histo_size; slice_num++){
	//Store the histogram clones for that slice locally
	hist_data = (TH1D*) slice_histo_data[slice_num]->Clone(Form("hist_data_%s_%i",myCutvar.getCutVar().Data(),slice_num));
	hist_p = (TH1D*) slice_histo_mc_p[slice_num]->Clone(Form("hist_p_%s_%i",myCutvar.getCutVar().Data(),slice_num));
	hist_n = (TH1D*) slice_histo_mc_n[slice_num]->Clone(Form("hist_n_%s_%i",myCutvar.getCutVar().Data(),slice_num));

	//Make a fit to the overall histogram for dx. In this case we used the MC proton and neutron shapes + a 2nd order poly background
	TFit *myFit = new TF1(Form("overall_fit_%s_%i",myCutvar.getCutVar().Data(),slice_num),fits::mc_p_n_poly2_slice_fit,xmin,xmax,initialParameters.size());

	Fit->SetParameters(initialParameters);
        // set parameter limits. SetParLimits -> (par#, min, max)
        Fit->SetParLimits(0, 0, 2000); // scale_p greater than 0
        Fit->SetParLimits(1, 0,2000); // scale_n greater than 0
        Fit->SetParLimits(2, -0.10,0.10); // shift_p less than +- 10cm
        Fit->SetParLimits(3, -0.10,0.10); // shift_n less than +- 10cm
        Fit->SetParLimits(6,-10000000,-0.0000000001); // x^2 term negative to force downward concavity 
        Fit->SetNpx(500);
	hist_data->GetXaxis()->SetRangeUser(xmin,xmax);

	//now fit the data histogram with the fit info
	hist_data->Fit(Fit,"Q R");

	//retrieve important fit results
	double scale_p  = Fit->GetParameter(0);
        double scale_p_err = Fit->GetParError(0);
        double scale_n  = Fit->GetParameter(1);
        double scale_n_err = Fit->GetParError(1);

        double shift_p= Fit->GetParameter(2);
        double shift_p_err= Fit->GetParError(2);
        double shift_n = Fit->GetParameter(3);
        double shift_n_err = Fit->GetParError(3);

        double ChiSq= Fit->GetChisquare();
        double ndf = Fit->GetNDF();

	vector<double> poly_result;
	vector<double> poly_result_err;

		for (int i =0 ; i < 3; i++){
          	poly_result.push_back( Fit->GetParameter(4+i) );
          	poly_result_err.push_back(Fit->GetParError(4+i) );
        	}

	//Compute physics results
	double Rsf = scale_n/scale_p;
	//just adding the uncert from the fit parameters in quadrature for now. Better method may be implemented later, might need to consider the whole fit.
        double Rsf_err = Rsf * sqrt( pow( (scale_n_err / scale_n), 2) + pow( (scale_p_err / scale_p),2) ); 

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

