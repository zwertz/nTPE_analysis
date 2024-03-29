#include "../include/plots.h"

//Author Ezekiel Wertz
//Implementations for plots related functions


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

//make the canvas to put fiducial/acceptance region check on
TCanvas* plotAcceptance_Fid_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TLine *LineFidPro,TH2D *hxy_globcut,TH2D *hxy_expect_cut,TH2D *hxy_expect_failedfid){

TCanvas* myCan = new TCanvas(name,name,1600,1200); 
myCan->Divide(3,1);

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

//second plot
myCan->cd(2);
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
//third plot
myCan->cd(3);
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

return myCan;
}

//make the canvas to put fiducial neutron and proton hypothesis check on
TCanvas* plotFid_Hypothesis_Check(const char *name,vector<TLine*> Lines_pos,vector<TLine*> Lines_aa,vector<TLine*> Lines_Fid,TLine *LineFidPro,TH2D *hxy_expect_fidcutn,TH2D *hxy_expect_fidcutp){

TCanvas* myCan = new TCanvas(name,name,1600,1200);
myCan->Divide(2,1);

//first plot
myCan->cd(1);
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

//second plot
myCan->cd(2);
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

return myCan;
}


}//end namespace
