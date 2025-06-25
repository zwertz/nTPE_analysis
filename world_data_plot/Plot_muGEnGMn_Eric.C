/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//// Parameterized Form Factor Central Value and Error
/////////////////////////////////////////////////////////
//// ID = 1 for GEp, 2 for GMp, 3 for GEn, 4 for GMn, 
//// Q2 in GeV^2
////
// The parameterization formula returns the uncertainty devided by G(0)*GD, where 
//  GD(Q2) = 1./(1+Q2/0.71)^2
// and GEp(0) = 1, GMp(0) = 2.79284356, GEn(0) = 1, GMn(0) = -1.91304272,
//
// The parameterization formula for the Form Factor value is:
//  $$ GN(z) = sum_{i=0}^{N=12}(a_i * z^i) 
// Note that the return value has been divided by (G(Q2=0)*G_Dip)
//
// The parameterization formula for the Form Factor error is:
// $$ log_{10}\frac{\delta G}{G_D} = (L+c_0)\Theta_a(L_1-L) 
//                                 +\sum_{i=1}^{N}(c_i+d_i L)[\Theta_a(L_i-L)-\Theta_a(L_{i+1}-L)]
//                                 +log_{10}(E_{\inf})\Theta_a(L-L_{N+1})$$
// where $L=log_{10}(Q^2)$, $\Theta_{a}(x)=[1+10^{-ax}]^{-1}$. $a=1$.

int GetFF(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err){/*{{{*/
    // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
    if (kID<1 || kID>4){
        cerr<<"*** ERROR***, kID is not any of [1->GEp, 2->GMp, 3->GEn, 4->GMn]"<<endl;
        GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
        return -1;
    }
    ////////////////////////////////////////////////
    //// z-Expansion Parameters for Form Factor Values
    /////////////////////////////////////////////////*{{{*/
    const double GN_Coef_Fit[4][13] ={
        {0.239163298067, -1.10985857441, 1.44438081306, 0.479569465603, -2.28689474187,  1.12663298498, 1.25061984354,-3.63102047159, 4.08221702379,  0.504097346499,  -5.08512046051,  3.96774254395,-0.981529071103}, /*GEp*/
        {0.264142994136, -1.09530612212, 1.21855378178, 0.661136493537, -1.40567892503, -1.35641843888, 1.44702915534, 4.2356697359, -5.33404565341, -2.91630052096,    8.70740306757, -5.70699994375, 1.28081437589}, /*GMp*/
        {0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, 0.300445258602,-0.661888687179,-0.175639769687, 0.624691724461,-0.077684299367,-0.236003975259, 0.090401973470, 0.0, 0.0}, /*GEn*/
        {0.257758326959,-1.079540642058, 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208, 2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0} /*GMn*/
    };/*}}}*/

    ////////////////////////////////////////////////
    //// Parameters for Form Factor Errors
    ////////////////////////////////////////////////   /*{{{*/
    const double parL[4][2] ={
        {-0.97775297,  0.99685273}, //GEp
        {-0.68452707,  0.99709151}, //GMp
        {-2.02311829, 1.00066282}, //GEn
        {-0.20765505, 0.99767103}, //GMn
    };
    const double parM[4][15] = {
        {  -1.97750308e+00,  -4.46566998e-01,   2.94508717e-01,   1.54467525e+00,
            9.05268347e-01,  -6.00008111e-01,  -1.10732394e+00,  -9.85982716e-02,
            4.63035988e-01,   1.37729116e-01,  -7.82991627e-02,  -3.63056932e-02,
            2.64219326e-03,   3.13261383e-03,   3.89593858e-04}, //GEp

        {  -1.76549673e+00,   1.67218457e-01,  -1.20542733e+00,  -4.72244127e-01,
            1.41548871e+00,   6.61320779e-01,  -8.16422909e-01,  -3.73804477e-01,
            2.62223992e-01,   1.28886639e-01,  -3.90901510e-02,  -2.44995181e-02,
            8.34270064e-04,   1.88226433e-03,   2.43073327e-04}, //GMp

        {  -2.07343771e+00,   1.13218347e+00,   1.03946682e+00,  -2.79708561e-01,
           -3.39166129e-01,   1.98498974e-01,  -1.45403679e-01,  -1.21705930e-01,
            1.14234312e-01,   5.69989513e-02,  -2.33664051e-02,  -1.35740738e-02,
            7.84044667e-04,   1.19890550e-03,   1.55012141e-04}, //GEn

        {  -2.07087611e+00,   4.32385770e-02,  -3.28705077e-01,   5.08142662e-01,
            1.89103676e+00,   1.36784324e-01,  -1.47078994e+00,  -3.54336795e-01,
            4.98368396e-01,   1.77178596e-01,  -7.34859451e-02,  -3.72184066e-02,
            1.97024963e-03,   2.88676628e-03,   3.57964735e-04} //GMn
    };
    const double parH[4][3] = {
        {0.78584754,  1.89052183, -0.4104746}, //GEp
        {0.80374002,  1.98005828, -0.69700928}, //GMp
        {0.4553596,  1.95063341,  0.32421279}, //GEn
        {0.50859057, 1.96863291,  0.2321395} //GMn 
    };
    /*}}}*/

    //// Apply the z-expansion formula for form factor/*{{{*/
    const double tcut = 0.0779191396 ;
    const double t0 = -0.7 ;
    double z = (sqrt(tcut+kQ2)-sqrt(tcut-t0))/(sqrt(tcut+kQ2)+sqrt(tcut-t0)) ;
    double GNQ2 = 0.0;
    for (int i=0;i<13;i++) GNQ2 += GN_Coef_Fit[kID-1][i] * pow(z, i);
    double GDip= pow(1./(1. + kQ2/0.71), 2);
    GNGD_Fit[0] = GNQ2 / GDip; //Note that Coef_Fit have been divided by mu_p or mu_n
    /*}}}*/

    //// Apply the parameterization formula for error/*{{{*/
    double lnQ2 = log10(kQ2);
    double lnGNGD_Err=0.0;
    if (kQ2<1e-3)
        lnGNGD_Err = parL[kID-1][0] + parL[kID-1][1]*lnQ2;
    else if (kQ2>1e2)
        lnGNGD_Err = parH[kID-1][0]*sqrt(lnQ2-parH[kID-1][1]) + parH[kID-1][2];
    else{
        for (int i=0; i<15;i++) lnGNGD_Err += parM[kID-1][i] * pow(lnQ2,i);
    }
    GNGD_Err[0] = pow(10.,(lnGNGD_Err));    //LOG10(dG/G(0)/GD);
    /*}}}*/

    return 0;
}/*}}}*/

const double LambdaD2 = 0.71;
double GD(double Q2){
  return( 1.0/(1+Q2/LambdaD2)/(1+Q2/LambdaD2) );
  
}

void PlotmuGEnGMn()
{
  double Q2data[100];
  double dQ2data[100];
  double GEnGDdata[100];
  double dGEnGDdata[100];

  double GMnGDmudata[100];
  double dGMnGDmudata[100];
  
  double muGEnGMndata[100];
  double mudGEnGMndata[100];

  double Q2fit[200];
  double dQ2fit[200];
  double GEnGDfit[200];
  double dGEnGDfit[200];

  double GMnGDmufit[200];
  double dGMnGDmufit[200];
  
  double muGEnGMnfit[200];
  double mudGEnGMnfit[200];
  
  ifstream in_gen("/w/halla-scshelf2102/sbs/efuchey/software_2023/work/PlotGEn/supplemental_materials/data/neutron/World_GEn.dat");
  TString currentline;

  int npts = 0;
  while( currentline.ReadLine(in_gen) ){
    if( !currentline.BeginsWith("#") ){
      Int_t ntokens = 0;
      TString stemp;
      std::unique_ptr<TObjArray> tokens( currentline.Tokenize(", \t") );
      if( !tokens->IsEmpty() ) {
	ntokens = tokens->GetLast()+1;
	if(ntokens==6){
	  stemp = ( (TObjString*) (*tokens)[0] )->GetString();
	  Q2data[npts] = stemp.Atof();
	  
	  stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  GEnGDdata[npts] = stemp.Atof();
	  
	  stemp = ( (TObjString*) (*tokens)[2] )->GetString();
	  dGEnGDdata[npts] = stemp.Atof();
	  npts++;
	}
      }
    }
  }
  const Int_t Npoints = npts;
  double GMnGDmu_temp, dGMnGDmu_temp;
  double GD_temp;
  for(int i = 0; i<Npoints; i++){
    dQ2data[i] = 0.0;
    GetFF(4, Q2data[i], &GMnGDmu_temp, &dGMnGDmu_temp);
    GD_temp = GD(Q2data[i]);
    GMnGDmudata[i] = GMnGDmu_temp;
    dGMnGDmudata[i] = dGMnGDmu_temp;
    
    cout << Q2data[i] << " " << GEnGDdata[i] << " " << dGEnGDdata[i] << " " << GMnGDmudata[i] << " " << dGMnGDmudata[i] << endl;

    muGEnGMndata[i] = GEnGDdata[i]/GD_temp/GMnGDmudata[i];
    mudGEnGMndata[i] = muGEnGMndata[i]*sqrt( dGEnGDdata[i]*dGEnGDdata[i]/GEnGDdata[i]/GEnGDdata[i] + dGMnGDmudata[i]*dGMnGDmudata[i]/GMnGDmudata[i]/GMnGDmudata[i] );
  }

  double GEnGD_temp, dGEnGD_temp;
  for(int i = 0; i<200; i++){
    dQ2fit[i] = 0.0;
    Q2fit[i] = i*0.05;
    GetFF(3, Q2fit[i], &GEnGD_temp, &dGEnGD_temp);
    GetFF(4, Q2fit[i], &GMnGDmu_temp, &dGMnGDmu_temp);
    
    GEnGDfit[i] = GEnGD_temp;
    dGEnGDfit[i] = dGEnGD_temp;
    
    GMnGDmufit[i] = GMnGDmu_temp;
    dGMnGDmufit[i] = dGMnGDmu_temp;
    
    muGEnGMnfit[i] = GEnGDfit[i]/GMnGDmufit[i];
    mudGEnGMnfit[i] = muGEnGMnfit[i]*sqrt( dGEnGDfit[i]*dGEnGDfit[i]/GEnGDfit[i]/GEnGDfit[i] + dGMnGDmufit[i]*dGMnGDmufit[i]/GMnGDmufit[i]/GMnGDmufit[i] );
    
  }

  double Q2_ntpe[1] = {4.5};
  double ErrQ2_ntpe[1] = {0.0};
  double GENslope_ntpe[1] = {0.5};
  double ErrGENslope_ntpe[1] = {0.1};
  
  TGraphErrors* grGEn_nTPE = new TGraphErrors(1, Q2_ntpe, GENslope_ntpe, ErrQ2_ntpe, ErrGENslope_ntpe);
  grGEn_nTPE->SetMarkerColor(8);
  grGEn_nTPE->SetMarkerStyle(21);
  
  TGraphErrors* grGEndata = new TGraphErrors(Npoints, Q2data, muGEnGMndata, dQ2data, mudGEnGMndata);
  grGEndata->SetMarkerColor(1);
  grGEndata->SetMarkerStyle(20);
  
  TGraph* grGEnfit = new TGraph(200, Q2fit, muGEnGMnfit);
  grGEnfit->SetLineColor(1);
  grGEnfit->SetLineWidth(2);
  
  TGraphErrors* grGEnfitErr = new TGraphErrors(200, Q2fit, muGEnGMnfit, dQ2fit, mudGEnGMnfit);
  grGEnfitErr->SetFillColor(kGray);
  grGEnfitErr->SetLineWidth(2);

  TLegend* Leg = new TLegend(0.45, 0.15, 0.85, 0.4);
  Leg->AddEntry(grGEndata, "GEn world data", "PE");
  Leg->AddEntry(grGEnfitErr, "Ye at al, [PLB 777 (2018) 8-15]", "LF");
  Leg->AddEntry(grGEn_nTPE, "GEn from Rosenbluth slope", "PE");
  Leg->SetMargin(0.125);
  Leg->SetFillStyle(0);
  Leg->SetBorderSize(0);
  
  TCanvas* C1 = new TCanvas();
  C1->SetGridx();
  C1->SetGridy();
  C1->cd();
  
  grGEnfit->Draw("AC");
  grGEnfitErr->Draw("C3, same");
  grGEnfit->Draw("C, same");
  grGEndata->Draw("P, same");
  grGEnfit->SetMinimum(0.0);
  grGEnfit->SetMaximum(1.0);
  grGEnfit->SetTitle(";Q^{2} (GeV/c)^{2};#mu G_{E}^{n}/G_{M}^{n}");
  grGEnfit->GetXaxis()->CenterTitle();
  grGEnfit->GetYaxis()->CenterTitle();
  grGEnfit->GetXaxis()->SetTitleOffset(0.7);
  grGEnfit->GetYaxis()->SetTitleOffset(0.7);
  grGEnfit->GetXaxis()->SetTitleSize(0.06);
  grGEnfit->GetYaxis()->SetTitleSize(0.06);
  grGEnfit->GetXaxis()->SetLabelSize(0.04);
  grGEnfit->GetYaxis()->SetLabelSize(0.04);
  Leg->Draw("same");
}
