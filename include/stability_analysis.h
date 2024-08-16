#ifndef STABILITY_ANALYSIS_H
#define STABILITY_ANALYSIS_H

//Author: Ezekiel Wertz
//This is a companion class to cutvar.h and cutvar.C. I ran into a challenge where I needed to be able to manipulate the information for all 3 histograms for data, mc p, and mc n. Building this into cutvar seemed difficult. So instead this class will take output, after some processing, from cutvar. Then do further cut stability study manipulation/evaluation here. This class should handle all information from data and MC at once. Will be modified as stability and systematic studies progress.

class stability_analysis{
private:
vector<TH1D*> slice_histo_data, slice_histo_mc_p, slice_histo_mc_n;


public:
//Standard constructor, it should initialize all the class variables we care about.
stability_analysis(vector<TH1D*> daDataHisto, vector<TH1D*> daMCPHisto, vector<TH1D*> daMCNHisto);


};//end class
#endif
