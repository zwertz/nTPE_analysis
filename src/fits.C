//Author Ezekiel Wertz
//Implementation of fit functions that are common throughout analysis

namespace fits{

//4th-order polynomial
double poly4_fit(double *x, double *param){

double yint = param[0];
double p1 = param[1];
double p2 = param[2];
double p3 = param[3];
double p4 = param[4];

double func = p4*(pow(x[0],4))+p3*(pow(x[0],3))+ p2*(pow(x[0],2))+p1*(x[0])+yint;
return func;
}

//function to fit the histogram then fine fit and return all desired parameters.
vector<pair<double,double>> fitAndFineFit(TH1D* histogram, const string& fitName, const string& fitFormula, int paramCount, double hcalfit_low, double hcalfit_high, pair<double,double>& fitqual, const string& fitOptions = "RBMQ0"){

//make the vector we will eventually return
vector<pair<double,double>> paramsAndErrs(paramCount);

//Make a fit based on the information provided
TF1* datFit = new TF1(fitName.c_str(),fitFormula.c_str(),hcalfit_low,hcalfit_high,paramCount);

	//Reset the parameters for the fit
	for(int i=0; i<paramCount; ++i){
	datFit->SetParameter(i,0);
	datFit->SetParError(i,0);
	}

//Fit the provided histogram with the fit we have
histogram->Fit(datFit,fitOptions.c_str());

//vector to hold fine fit initial parameters
vector<double> fineFitParams(paramCount);
	//loop over the information from the fit and store it in our parameter and error vector
	for(int j=0; j<paramCount; ++j){
	paramsAndErrs[j].first = datFit->GetParameter(j);
	paramsAndErrs[j].second = datFit->GetParError(j);
	fineFitParams[j] = paramsAndErrs[j].first;
	}

//Fine fit to get a better set of parameters
TF1* datFineFit = new TF1((fitName+"_fine").c_str(),fitFormula.c_str(),hcalfit_low,hcalfit_high,paramCount);

//set the fine fit intiall parameters to be the ones we just found
datFineFit->SetParameters(fineFitParams.data());

//Fit the histogram again but with these better intial parameters
histogram->Fit(datFineFit,fitOptions.c_str());

	//update the parameters and errors with the ones we just found from fine fit
	for(int k=0; k<paramCount; ++k){
	paramsAndErrs[k].first = datFineFit->GetParameter(k);
        paramsAndErrs[k].second = datFineFit->GetParError(k);
	//Diagnostic
	cout << datFineFit->GetParameter(k) << endl;
	}

//don't like these lines because it kind of circumvents the idea of a return
//get the Chisquare and NDF from the fit which is useful for analysis
fitqual.first = datFineFit->GetChisquare();
fitqual.second = datFineFit->GetNDF();



//cleanup
delete datFit;
delete datFineFit;
return paramsAndErrs;	
}

//Get the variance on a fit from the quadrature sum of the parameter errors
double getFitError(TF1* fit){
double SumSquare = 0.0;

	//loop over the array of parameter errors. Square them and then add that to the sum of squares
	for(int i=0; i < fit->GetNpar();i++ ){
	double ParError = fit->GetParError(i);
	cout << "Par Error: "<< i << " " << ParError << endl;
	SumSquare += ParError*ParError;
	}
	//As of here SumSquare should be the Param Errors Squared then Summed
	//Just take the Square root of SumSquare to get Quad Sum
	double QuadSum = sqrt(SumSquare);
return QuadSum;
}

//Assume the input fit ptr is for the total fit
//Get the variance on a fit from the quadrature sum of the diagonal entries. Which the diagonal entries should be the variance on each parameter.
double getFitError(TFitResultPtr fit_ptr){
TMatrixD fit_covariance = fit_ptr->GetCovarianceMatrix();
TMatrixD bg_mat = fit_covariance.GetSub(4,8,4,8);
TVectorD diagonal = TMatrixDDiag(bg_mat);
double SumSquare = 0.0;
	
	//loop over the elements of the vector. Sum the variances
	for(int j=0; j < diagonal.GetNrows();j++){
	double value = diagonal[j];
	cout << "Error from Covariance Matrix: " << j << " " << value << endl;
	SumSquare += value;
	}

	double QuadSum =sqrt(SumSquare);
return QuadSum;
}

}//end namepspace
