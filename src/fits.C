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




}//end namepspace
