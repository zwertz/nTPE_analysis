
//This inherently relies on the parameterization and world data made accessible from the supplementary material of the article by Ye et. al in 2018. Takes the GEn world data and combines it with the GMn parameterization info to form a plot of muGEn/GMn for the world data and my own extracted point.

#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include "GetFF.C"

const double LambdaD2 = 0.71;
double GD(double Q2){

	double GD = pow((1+ Q2/LambdaD2),-2); 
	return GD;

}

void plot_muGEnGMn_mine(){



  //vectors to store info for data
  vector<double> Q2_data;
  vector<double> dQ2_data;
  vector<double> GEn_data;
  vector<double> dGEn_data;


  vector<double> GMnGDmu_data;
  vector<double> dGMnGDmu_data;

  vector<double> muGEnGMn_data;
  vector<double> dmuGEnGMn_data;

  //vectors to store info for world fit
  vector<double> Q2_fit;
  vector<double> dQ2_fit;
  vector<double> GEnGD_fit;
  vector<double> dGEnGD_fit;

  vector<double> GMnGDmu_fit;
  vector<double> dGMnGDmu_fit;

  vector<double> muGEnGMn_fit;
  vector<double> dmuGEnGMn_fit;

  //Start io stream to read-in world GEn data. NOTE it is labeled as GEn/GD. This is wrong it is just GEn. I have verified this with the referenced data/articles.
  //We will also neglect Experiments marked with the number 7 which seem to be derived from theoretical extractions of elastic electron deuteron scattering. Need to check this topic. A version of this dat is available from the Ye Supplement.
  ifstream in_GEn("World_GEn_no7.dat");
  TString currentline;

  //Number of data points from file
  int num_points = 0;

  //loop over the entries in the file and store the data
  while( currentline.ReadLine(in_GEn) ){
  	if( !currentline.BeginsWith("#") ){
	Int_t ntokens = 0;
	TString stemp;
	std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
		if( !tokens->IsEmpty() ) {
        	ntokens = tokens->GetLast()+1;
			if(ntokens==6){
			Q2_data.push_back((((TObjString*) (*tokens)[0] )->GetString()).Atof());
			dQ2_data.push_back(0.0);
			GEn_data.push_back((((TObjString*) (*tokens)[1] )->GetString()).Atof());
			dGEn_data.push_back((((TObjString*) (*tokens)[2] )->GetString()).Atof());
			num_points++;
        		}
		}
	}
  }

  //Calculate the corresponding GMn/mu_n values and uncertainies for the same Q^2 as the data points from the parameterization
  //Then determine mun GEn/GMn from (GEn/GD)/(GMN/munGD)
  double GMnGDmu_temp, dGMnGDmu_temp;
  for(int i = 0; i<num_points; i++){
    
    GetFF(4, Q2_data[i], &GMnGDmu_temp, &dGMnGDmu_temp);
    double GD_temp = GD(Q2_data[i]);
    GMnGDmu_data.push_back(GMnGDmu_temp);
    dGMnGDmu_data.push_back(dGMnGDmu_temp);

    double munGEnGMn =(GEn_data[i]/GD_temp)/GMnGDmu_data[i];
    muGEnGMn_data.push_back(munGEnGMn);
    double munGEnGMn_error = muGEnGMn_data[i]*sqrt( dGEn_data[i]*dGEn_data[i]/GEn_data[i]/GEn_data[i] + dGMnGDmu_data[i]*dGMnGDmu_data[i]/GMnGDmu_data[i]/GMnGDmu_data[i] );
    dmuGEnGMn_data.push_back(munGEnGMn_error);

    cout << Q2_data[i] << " GEn: " << GEn_data[i] << " +/- " << dGEn_data[i] << " GMn/muGD " << GMnGDmu_data[i] << " +/- " << dGMnGDmu_data[i] << " muGEn/GMn: " << muGEnGMn_data[i] << " +/- " << dmuGEnGMn_data[i] << endl;
  }

  //Handle the information that will be used to make the world data fit
  double GEnGD_fit_temp, dGEnGD_fit_temp;
   double GMnGDmu_fit_temp, dGMnGDmu_fit_temp;
  //Lets assume we want to break 15 GeV^2 into 500 bins That means each bin is 0.03 Q^2
  int num_bins = 500; 
  for(int i = 0; i<num_bins; i++){
    dQ2_fit.push_back(0.0);
    Q2_fit.push_back(i*0.03);

    GetFF(3, Q2_fit[i], &GEnGD_fit_temp, &dGEnGD_fit_temp);
    GetFF(4, Q2_fit[i], &GMnGDmu_fit_temp, &dGMnGDmu_fit_temp);

    GEnGD_fit.push_back(GEnGD_fit_temp);
    dGEnGD_fit.push_back(dGEnGD_fit_temp);

    GMnGDmu_fit.push_back(GMnGDmu_fit_temp);
    dGMnGDmu_fit.push_back(dGMnGDmu_fit_temp);

    muGEnGMn_fit.push_back(GEnGD_fit[i]/GMnGDmu_fit[i]);
    dmuGEnGMn_fit.push_back(muGEnGMn_fit[i]*sqrt( dGEnGD_fit[i]*dGEnGD_fit[i]/GEnGD_fit[i]/GEnGD_fit[i] + dGMnGDmu_fit[i]*dGMnGDmu_fit[i]/GMnGDmu_fit[i]/GMnGDmu_fit[i] ));

  }

  //My Extracted preliminary value
  double Q2_ntpe[1] = {4.48};
  double ErrQ2_ntpe[1] = {0.0};
  double muGEnGMN_ntpe[1] = {0.652};
  double ErrmuGEnGMn_ntpe[1] = {0.16};

  //Rough TPE correction applied My Extracted preliminary value
  double Q2_TPE[1] = {4.48};
  double ErrQ2_TPE[1] = {0.0};
  double muGEnGMN_TPE[1] = {0.867};
  double ErrmuGEnGMn_TPE[1] = {0.16};


  //Convert the relevant vectors to be used as arrays
  double * Q2_data_array = Q2_data.data();
  double * dQ2_data_array = dQ2_data.data();

  double * muGEnGMn_data_array = muGEnGMn_data.data();
  double * dmuGEnGMn_data_array = dmuGEnGMn_data.data();

  double * Q2_fit_array = Q2_fit.data();
  double * dQ2_fit_array = dQ2_fit.data();

  double * muGEnGMn_fit_array = muGEnGMn_fit.data();
  double * dmuGEnGMn_fit_array = dmuGEnGMn_fit.data();


  TF1 *zero = new TF1("zero","0",0.0,15.0);
  zero->SetLineWidth(2);
  zero->SetLineColor(1);
  zero->SetLineStyle(2);

  TF1 *one = new TF1("one","1",0.0,15.0);
  one->SetLineWidth(2);
  one->SetLineColor(1);
  one->SetLineStyle(2);

  //Graph my value info
  TGraphErrors* graph_muGEnGMN_ntpe = new TGraphErrors(1, Q2_ntpe, muGEnGMN_ntpe, ErrQ2_ntpe, ErrmuGEnGMn_ntpe);
  graph_muGEnGMN_ntpe->SetMarkerColor(2);
  graph_muGEnGMN_ntpe->SetMarkerStyle(21);
  graph_muGEnGMN_ntpe->SetLineColor(2);
  graph_muGEnGMN_ntpe->SetMarkerSize(2);
  graph_muGEnGMN_ntpe->SetLineWidth(1);

  //Graph my rough TPE corrected value
  TGraphErrors* graph_muGEnGMN_TPE = new TGraphErrors(1, Q2_TPE, muGEnGMN_TPE, ErrQ2_TPE, ErrmuGEnGMn_TPE);
  graph_muGEnGMN_TPE->SetMarkerColor(4);
  graph_muGEnGMN_TPE->SetMarkerStyle(22);
  graph_muGEnGMN_TPE->SetLineColor(4);
  graph_muGEnGMN_TPE->SetMarkerSize(2.65);
  graph_muGEnGMN_TPE->SetLineWidth(1);

  //Graph rest of world data points
  TGraphErrors* graph_world_data = new TGraphErrors(num_points, Q2_data_array, muGEnGMn_data_array, dQ2_data_array, dmuGEnGMn_data_array);
  graph_world_data->SetMarkerColor(1);
  graph_world_data->SetMarkerStyle(20);
  graph_world_data->SetMarkerSize(2);
  graph_world_data->SetLineWidth(1);

  //Graph the world data fit
  TGraph* graph_world_fit = new TGraph(num_bins, Q2_fit_array, muGEnGMn_fit_array);
  graph_world_fit->SetLineColor(1);
  graph_world_fit->SetLineWidth(2);

  //Error band on the world data fit
  TGraphErrors* graph_world_fit_Err = new TGraphErrors(num_bins, Q2_fit_array, muGEnGMn_fit_array, dQ2_fit_array, dmuGEnGMn_fit_array);
  graph_world_fit_Err->SetFillColor(kGray);
  graph_world_fit_Err->SetFillColorAlpha(kGray,0.5);
  graph_world_fit_Err->SetLineWidth(2);

  TLegend* Leg1 = new TLegend(0.15, 0.15, 0.55, 0.3);
  Leg1->AddEntry(graph_world_data, "World Data (polarization)", "PE");
  Leg1->AddEntry(graph_muGEnGMN_ntpe, "This Work (Rosenbluth)", "PE");
  Leg1->SetMargin(0.125);
  Leg1->SetFillStyle(0);
  Leg1->SetBorderSize(0);

  TLegend* Leg2 = new TLegend(0.55, 0.15, 0.85, 0.3);
  Leg2->AddEntry(graph_world_fit_Err, "Global Fit (Ye 2018)", "LF");
  Leg2->SetMargin(0.125);
  Leg2->SetFillStyle(0);
  Leg2->SetBorderSize(0);

  TLegend* Leg3 = new TLegend(0.15, 0.15, 0.55, 0.3);
  Leg3->AddEntry(graph_world_data, "World Data (polarization)", "PE");
  Leg3->AddEntry(graph_muGEnGMN_TPE, "This Work, #approx TPE corrected", "PE");
  Leg3->SetMargin(0.125);
  Leg3->SetFillStyle(0);
  Leg3->SetBorderSize(0);

  TLegend* Leg4 = new TLegend(0.55, 0.15, 0.85, 0.3);
  Leg4->AddEntry(graph_world_fit_Err, "Global Fit (Ye 2018)", "LF");
  Leg4->SetMargin(0.125);
  Leg4->SetFillStyle(0);
  Leg4->SetBorderSize(0);


  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //My value with the world data
  TCanvas *c1 = new TCanvas("","",1600,1200);
  c1->cd();
  graph_world_fit->GetXaxis()->SetRangeUser(0.0,15.0);
  graph_world_fit->GetXaxis()->SetNdivisions(504);
  graph_world_fit->GetYaxis()->SetRangeUser(-0.5,1.5);
  graph_world_fit->GetYaxis()->SetNdivisions(505);
  graph_world_fit->Draw("AC");
  graph_world_fit_Err->Draw("C3, same");
  zero->Draw("C, same");
  one->Draw("C, same");
  graph_world_fit->Draw("C, same");
  graph_world_data->Draw("P, same");
  graph_muGEnGMN_ntpe->Draw("P, same");
  graph_world_fit->SetTitle(";Q^{2} (GeV/c)^{2};#mu_{n} G_{E}^{n}/G_{M}^{n}");
  graph_world_fit->GetXaxis()->CenterTitle();
  graph_world_fit->GetYaxis()->CenterTitle();
  graph_world_fit->GetXaxis()->SetTitleOffset(0.65);
  graph_world_fit->GetYaxis()->SetTitleOffset(0.7);
  graph_world_fit->GetXaxis()->SetTitleSize(0.06);
  graph_world_fit->GetYaxis()->SetTitleSize(0.06);
  graph_world_fit->GetXaxis()->SetLabelSize(0.04);
  graph_world_fit->GetYaxis()->SetLabelSize(0.04);
  Leg1->Draw("same");
  Leg2->Draw("same");

  //My rough TPE corrected value with the world data
  TCanvas *c2 = new TCanvas("","",1600,1200);
  c2->cd();
  graph_world_fit->GetXaxis()->SetRangeUser(0.0,15.0);
  graph_world_fit->GetXaxis()->SetNdivisions(504);
  graph_world_fit->GetYaxis()->SetRangeUser(-0.5,1.5);
  graph_world_fit->GetYaxis()->SetNdivisions(505);
  graph_world_fit->Draw("AC");
  graph_world_fit_Err->Draw("C3, same");
  zero->Draw("C, same");
  one->Draw("C, same");
  graph_world_fit->Draw("C, same");
  graph_world_data->Draw("P, same");
  graph_muGEnGMN_TPE->Draw("P, same");
  graph_world_fit->SetTitle(";Q^{2} (GeV/c)^{2};#mu_{n} G_{E}^{n}/G_{M}^{n}");
  graph_world_fit->GetXaxis()->CenterTitle();
  graph_world_fit->GetYaxis()->CenterTitle();
  graph_world_fit->GetXaxis()->SetTitleOffset(0.65);
  graph_world_fit->GetYaxis()->SetTitleOffset(0.7);
  graph_world_fit->GetXaxis()->SetTitleSize(0.06);
  graph_world_fit->GetYaxis()->SetTitleSize(0.06);
  graph_world_fit->GetXaxis()->SetLabelSize(0.04);
  graph_world_fit->GetYaxis()->SetLabelSize(0.04);
  Leg3->Draw("same");
  Leg4->Draw("same");

  //Add a save to a pdf file
  //Outfile
  TString plotname = "Zeke_muGEnGMn_plots.pdf";
  TString start = Form("%s%s",plotname.Data(),"(");
  //middle is the same as the name
  TString end = Form("%s%s",plotname.Data(),")");
  c1->Print(start.Data(),"pdf");
  c2->Print(end.Data(),"pdf");


}//end main

