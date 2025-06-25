//This inherently relies on the parameterization and world data made accessible from the supplementary material of the article by Ye et. al in 2018.Takes the GMn world data and plots it and my extracted values.
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include "GetFF.C"

void plot_GMnmuGD_mine(){

//vectors to store info for data
  vector<double> Q2_data;
  vector<double> dQ2_data;
  
  vector<double> GMnGDmu_data;
  vector<double> dGMnGDmu_data;

  //vectors to store info for world fit
  vector<double> Q2_fit;
  vector<double> dQ2_fit;

  vector<double> GMnGDmu_fit;
  vector<double> dGMnGDmu_fit;

  //Start io stream to read-in world GMn data. NOTE it is labeled as GMn/mun GD
  ifstream in_GMn("World_GMn_50QCD.dat");
  TString currentline;

  //Number of data points from file
  int num_points = 0;

  //loop over the entries in the file and store the data
  while( currentline.ReadLine(in_GMn) ){
        if( !currentline.BeginsWith("#") ){
        Int_t ntokens = 0;
        TString stemp;
        std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
                if( !tokens->IsEmpty() ) {
                ntokens = tokens->GetLast()+1;
		
		//if(ntokens==6){
                        if(ntokens==3){
			Q2_data.push_back((((TObjString*) (*tokens)[0] )->GetString()).Atof());
                        dQ2_data.push_back(0.0);
                        GMnGDmu_data.push_back((((TObjString*) (*tokens)[1] )->GetString()).Atof());
                        dGMnGDmu_data.push_back((((TObjString*) (*tokens)[2] )->GetString()).Atof());
                        num_points++;
                        }
                }
        }
  }

//Handle the information that will be used to make the world data fit
   double GMnGDmu_fit_temp, dGMnGDmu_fit_temp;
  //Lets assume we want to break 15 GeV^2 into 500 bins That means each bin is 0.03 Q^2
  int num_bins = 500;
  for(int i = 0; i<num_bins; i++){
    dQ2_fit.push_back(0.0);
    Q2_fit.push_back(i*0.03);

    GetFF(4, Q2_fit[i], &GMnGDmu_fit_temp, &dGMnGDmu_fit_temp);


    GMnGDmu_fit.push_back(GMnGDmu_fit_temp);
    dGMnGDmu_fit.push_back(dGMnGDmu_fit_temp);

  }
  //My Extracted preliminary value
  double Q2_SBS850[1] = {4.48};
  double ErrQ2_SBS850[1] = {0.0};
  double GMNmuGD_SBS850[1] = {0.9582};
  double ErrGMNmuGD_SBS850[1] = {0.0142};

  double Q2_SBS870[1] = {4.48};
  double ErrQ2_SBS870[1] = {0.0};
  double GMNmuGD_SBS870[1] = {0.9556};
  double ErrGMNmuGD_SBS870[1] = {0.0134};

  double Q2_SBS8100[1] = {4.48};
  double ErrQ2_SBS8100[1] = {0.0};
  double GMNmuGD_SBS8100[1] = {0.9518};
  double ErrGMNmuGD_SBS8100[1] = {0.0138};

  //True value
  //double Q2_SBS8WM[1] = {4.48};
  //Ofset for clarity
  double Q2_SBS8WM[1] = {4.489};
  double ErrQ2_SBS8WM[1] = {0.0};
  double GMNmuGD_SBS8WM[1] = {0.9546};
  double ErrGMNmuGD_SBS8WM[1] = {0.0132};


  //True value
  //double Q2_SBS970[1] = {4.476};
  //Ofset for clarity
  double Q2_SBS970[1] = {4.476};
  double ErrQ2_SBS970[1] = {0.0};
  double GMNmuGD_SBS970[1] = {0.9563};
  double ErrGMNmuGD_SBS970[1] = {0.0110 };

  //Convert the relevant vectors to be used as arrays
  double * Q2_data_array = Q2_data.data();
  double * dQ2_data_array = dQ2_data.data();

  double * GMnGDmu_data_array = GMnGDmu_data.data();
  double * dGMnGDmu_data_array = dGMnGDmu_data.data();

  double * Q2_fit_array = Q2_fit.data();
  double * dQ2_fit_array = dQ2_fit.data();

  double * GMnGDmu_fit_array = GMnGDmu_fit.data();
  double * dGMnGDmu_fit_array = dGMnGDmu_fit.data();

   TF1 *one = new TF1("one","1",0.0,15.0);
  one->SetLineWidth(2);
  one->SetLineColor(1);
  one->SetLineStyle(2);

  TF1 *one_sec = new TF1("one_sec","1",3.5,5.0);
  one_sec->SetLineWidth(2);
  one_sec->SetLineColor(1);
  one_sec->SetLineStyle(2);

  //Graph my value info one by one
  TGraphErrors* graph_SBS850 = new TGraphErrors(1, Q2_SBS850, GMNmuGD_SBS850, ErrQ2_SBS850, ErrGMNmuGD_SBS850);
  graph_SBS850->SetMarkerColor(4);
  graph_SBS850->SetMarkerStyle(24);
  graph_SBS850->SetLineColor(4);
  graph_SBS850->SetMarkerSize(1.5);
  graph_SBS850->SetLineWidth(1);

  TGraphErrors* graph_SBS870 = new TGraphErrors(1, Q2_SBS870, GMNmuGD_SBS870, ErrQ2_SBS870, ErrGMNmuGD_SBS870);
  graph_SBS870->SetMarkerColor(4);
  graph_SBS870->SetMarkerStyle(25);
  graph_SBS870->SetLineColor(4);
  graph_SBS870->SetMarkerSize(1.5);
  graph_SBS870->SetLineWidth(1);

  TGraphErrors* graph_SBS8100 = new TGraphErrors(1, Q2_SBS8100, GMNmuGD_SBS8100, ErrQ2_SBS8100, ErrGMNmuGD_SBS8100);
  graph_SBS8100->SetMarkerColor(4);
  graph_SBS8100->SetMarkerStyle(26);
  graph_SBS8100->SetLineColor(4);
  graph_SBS8100->SetMarkerSize(1.5);
  graph_SBS8100->SetLineWidth(1);

  TGraphErrors* graph_SBS8WM = new TGraphErrors(1, Q2_SBS8WM, GMNmuGD_SBS8WM, ErrQ2_SBS8WM, ErrGMNmuGD_SBS8WM);
  graph_SBS8WM->SetMarkerColor(6);
  graph_SBS8WM->SetMarkerStyle(20);
  graph_SBS8WM->SetLineColor(6);
  graph_SBS8WM->SetMarkerSize(1.5);
  graph_SBS8WM->SetLineWidth(1);

  TGraphErrors* graph_SBS970 = new TGraphErrors(1, Q2_SBS970, GMNmuGD_SBS970, ErrQ2_SBS970, ErrGMNmuGD_SBS970);
  graph_SBS970->SetMarkerColor(6);
  graph_SBS970->SetMarkerStyle(23);
  graph_SBS970->SetLineColor(6);
  graph_SBS970->SetMarkerSize(1.5);
  graph_SBS970->SetLineWidth(1);

  //Graph rest of world data points
  TGraphErrors* graph_world_data = new TGraphErrors(num_points, Q2_data_array, GMnGDmu_data_array, dQ2_data_array, dGMnGDmu_data_array);
  graph_world_data->SetMarkerColor(1);
  graph_world_data->SetMarkerStyle(20);
  graph_world_data->SetMarkerSize(1.5);
  graph_world_data->SetLineWidth(1);

  //Graph the world data fit
  TGraph* graph_world_fit = new TGraph(num_bins, Q2_fit_array, GMnGDmu_fit_array);
  graph_world_fit->SetLineColor(1);
  graph_world_fit->SetLineWidth(2);

  TGraph* graph_world_fit_zoom = new TGraph(num_bins, Q2_fit_array, GMnGDmu_fit_array);
  graph_world_fit_zoom->SetLineColor(1);
  graph_world_fit_zoom->SetLineWidth(2);


  //Error band on the world data fit
  TGraphErrors* graph_world_fit_Err = new TGraphErrors(num_bins, Q2_fit_array, GMnGDmu_fit_array, dQ2_fit_array, dGMnGDmu_fit_array);
  graph_world_fit_Err->SetFillColor(kGray);
  graph_world_fit_Err->SetFillColorAlpha(kGray,0.5);
  graph_world_fit_Err->SetLineWidth(2);

  TLegend* Leg1 = new TLegend(0.15, 0.15, 0.6, 0.3);
  Leg1->AddEntry(graph_world_data, "World Data", "PE");
  //Leg1->AddEntry(graph_SBS850, "This Work: SBS-8 50\% field", "PE");
  //Leg1->AddEntry(graph_SBS870, "This Work: SBS-8 70\% field", "PE");
  //Leg1->AddEntry(graph_SBS8100, "This Work: SBS-8 100\% field", "PE");
  Leg1->AddEntry(graph_SBS8WM, "This Work: SBS-8 Weighted Mean", "PE");
  Leg1->AddEntry(graph_SBS970, "This Work: SBS-9 70\% field", "PE");
  Leg1->SetMargin(0.125);
  Leg1->SetFillStyle(0);
  Leg1->SetBorderSize(0);

  TLegend* Leg3 = new TLegend(0.15, 0.12, 0.6, 0.27);
  Leg3->AddEntry(graph_world_data, "World Data", "PE");
  //Leg3->AddEntry(graph_SBS850, "This Work: SBS-8 50\% field", "PE");
  //Leg3->AddEntry(graph_SBS870, "This Work: SBS-8 70\% field", "PE");
  //Leg3->AddEntry(graph_SBS8100, "This Work: SBS-8 100\% field", "PE");
  Leg3->AddEntry(graph_SBS8WM, "This Work: SBS-8 Weighted Mean", "PE");
  Leg3->AddEntry(graph_SBS970, "This Work: SBS-9 70\% field", "PE");
  Leg3->SetMargin(0.125);
  Leg3->SetFillStyle(0);
  Leg3->SetBorderSize(0);

  TLegend* Leg2 = new TLegend(0.6, 0.15, 0.85, 0.3);
  Leg2->AddEntry(graph_world_fit_Err, "Global Fit (Ye 2018)", "LF");
  Leg2->SetMargin(0.125);
  Leg2->SetFillStyle(0);
  Leg2->SetBorderSize(0);


  TLegend* Leg4 = new TLegend(0.6, 0.12, 0.85, 0.27);
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
  graph_world_fit->GetYaxis()->SetRangeUser(0.0,1.25);
  graph_world_fit->GetYaxis()->SetNdivisions(506);
  graph_world_fit->Draw("AC");
  graph_world_fit_Err->Draw("C3, same");
  one->Draw("C, same");
  graph_world_fit->Draw("C, same");
  graph_world_data->Draw("P, same");
  //graph_SBS850->Draw("P, same");
  //graph_SBS870->Draw("P, same");
  //graph_SBS8100->Draw("P, same");
  graph_SBS8WM->Draw("P, same");
  graph_SBS970->Draw("P, same");
  graph_world_fit->SetTitle(";Q^{2} (GeV/c)^{2}; G_{M}^{n}/(#mu_{n} G_{D})");
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

  TCanvas *c2 = new TCanvas("","",1600,1200);
  c2->cd();
  graph_world_fit_zoom->GetXaxis()->SetRangeUser(3.5,5.0);
  graph_world_fit_zoom->GetXaxis()->SetNdivisions(504);
  graph_world_fit_zoom->GetYaxis()->SetRangeUser(0.8,1.2);
  graph_world_fit_zoom->GetYaxis()->SetNdivisions(506);
  graph_world_fit_zoom->Draw("AC");
  graph_world_fit_Err->Draw("C3, same");
  one_sec->Draw("C, same");
  graph_world_fit_zoom->Draw("C, same");
  graph_world_data->Draw("P, same");
  //graph_SBS850->Draw("P, same");
  //graph_SBS870->Draw("P, same");
  //graph_SBS8100->Draw("P, same");
  graph_SBS8WM->Draw("P, same");
  graph_SBS970->Draw("P, same");
  graph_world_fit_zoom->SetTitle(";Q^{2} (GeV/c)^{2}; G_{M}^{n}/(#mu_{n} G_{D})");
  graph_world_fit_zoom->GetXaxis()->CenterTitle();
  graph_world_fit_zoom->GetYaxis()->CenterTitle();
  graph_world_fit_zoom->GetXaxis()->SetTitleOffset(0.65);
  graph_world_fit_zoom->GetYaxis()->SetTitleOffset(0.7);
  graph_world_fit_zoom->GetXaxis()->SetTitleSize(0.06);
  graph_world_fit_zoom->GetYaxis()->SetTitleSize(0.06);
  graph_world_fit_zoom->GetXaxis()->SetLabelSize(0.04);
  graph_world_fit_zoom->GetYaxis()->SetLabelSize(0.04);
  Leg3->Draw("same");
  Leg4->Draw("same");

  //Add a save to a pdf file
  //Outfile
  TString plotname = "Zeke_GMnmuGD_plots.pdf";
  TString start = Form("%s%s",plotname.Data(),"(");
  //middle is the same as the name
  TString end = Form("%s%s",plotname.Data(),")");
  //c1->Print(plotname.Data(),"pdf");
  c1->Print(start.Data(),"pdf");
  c2->Print(end.Data(),"pdf");

}//end main
