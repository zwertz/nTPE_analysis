//Author Ezekiel Wertz
//Implementation of fit functions that are common throughout analysis
//Largely outdated and depreciated due to the fit_histogram class. Ideally I should come back to this class and see if anything still uses it.

namespace fits{

//Get the variance on a fit from the quadrature sum of the parameter errors
double getFitError(TF1* fit){
double SumSquare = 0.0;

	//loop over the array of parameter errors. Square them and then add that to the sum of squares
	for(int i=0; i < fit->GetNpar();i++ ){
	double ParError = fit->GetParError(i);
	//cout << "Par Error: "<< i << " " << ParError << endl;
	SumSquare += ParError*ParError;
	}
	//As of here SumSquare should be the Param Errors Squared then Summed
	//Just take the Square root of SumSquare to get Quad Sum
	double QuadSum = sqrt(SumSquare);
return QuadSum;
}

//Assume the input fit ptr is for the total fit
//We want the variance on the fit for just the background. Get the sub matrix of the covariance matrix for the background fit. Return the fit error function.
double FitErrorFunc(double x,TFitResultPtr fit_ptr){
TMatrixD fit_covariance = fit_ptr->GetCovarianceMatrix();
TMatrixD bg_mat = fit_covariance.GetSub(4,8,4,8);
TVectorD diagonal = TMatrixDDiag(bg_mat);

//Check that Vij = Vji for all cases
for(int i = 0; i < bg_mat.GetNrows(); i++){
	for(int j = 0; j < bg_mat.GetNrows(); j++){
	double Vij = bg_mat(i,j);
	double Vji = bg_mat(j,i);

		if( Vij != Vji){
		cout << "Warning Vij: " << Vij << " does not equal Vji: " << Vji <<" figure out what is going on!" << endl;
		}
	}
}

//Based on error propagation for a fit function per bin. Based on David's notes and stats book. General form is available there. This is for the specific case of a fourth order polynomial
double Vaa = bg_mat(0,0); //sigma_a^2
double Vbb = bg_mat(1,1); //sigma_b^2
double Vcc = bg_mat(2,2); //sigma_c^2
double Vdd = bg_mat(3,3); //sigma_d^2
double Vee = bg_mat(4,4); //sigma_e^2

double Vab = bg_mat(0,1);
double Vac = bg_mat(0,2);
double Vad = bg_mat(0,3);
double Vae = bg_mat(0,4);

double Vbc = bg_mat(1,2);
double Vbd = bg_mat(1,3);
double Vbe = bg_mat(1,4);

double Vcd = bg_mat(2,3);
double Vce = bg_mat(2,4);
double Vde = bg_mat(3,4);

//Correct formula but causes dereferencing
long double fit_error =TMath::Sqrt(Vaa + pow(x,2)*Vbb + pow(x,4)*Vcc + pow(x,6)*Vdd + pow(x,8)*Vee + 2*x*Vab + 2*pow(x,2)*Vac + 2*pow(x,3)*Vad + 2*pow(x,4)*Vae + 2*pow(x,3)*Vbc + 2*pow(x,4)*Vbd + 2*pow(x,5)*Vbe + 2*pow(x,5)*Vcd + 2*pow(x,6)*Vce + 2*pow(x,7)*Vde  );

return fit_error;
}

}//end namepspace
