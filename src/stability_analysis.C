//Author: Ezekiel Wertz
//This is a companion class to cutvar.h and cutvar.C. I ran into a challenge where I needed to be able to manipulate the information for all 3 histograms for data, mc p, and mc n.Building this into cutvar seemed difficult. So instead this class will take output, after some processing, from cutvar. Then do further cut stability study manipulation/evaluation here. This class should handle all information from data and MC at once. Will be modified as stability and systematic studies progress.

stability_analysis(vector<TH1D*> daDataHisto, vector<TH1D*> daMCPHisto, vector<TH1D*> daMCNHisto){

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

}
