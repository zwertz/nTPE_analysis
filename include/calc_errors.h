#ifndef CALC_ERRORS_H
#define CALC_ERRORS_H

//Imported from John Boyds files. Really this should be both a header and a source file. For now will maintain current standing. Has some useful functions will reorganize it later, if necessary. So far this header file has not been incorporated into this analysis framework. And likely will not be.

//Error propagation for addition and subtraction:
double CalculateErrorAdditionSubtraction( double val1, double val1_err, double val2, double val2_err){
    double error = std::sqrt(std::pow(val1_err, 2) + std::pow(val2_err, 2));
    return error;
}

//Error propagation for multiplication and division:
double CalculateErrorMultiplicationDivision( double val1, double val1_err, double val2, double val2_err, double eval ){
    double relativeError = std::sqrt(std::pow(val1_err / val1, 2) + std::pow(val2_err / val2, 2));
    double finalError = eval*relativeError;

    return finalError;
}

double calculateVectorRSS( const std::vector<double>& values ){
    int n = values.size();
    double sum = 0.0;

    for( int i = 0; i < n; i++ ){
        sum += values[i]*values[i];
    }

    return sqrt(sum);
}

double calculateVectorRSSinRange( const std::vector<double>& values, int start_index, int end_index ){
    int n = (end_index - start_index) + 1;
    double sum = 0.0;

    for( int i = start_index; i <= end_index; i++ ){
        sum += values[i]*values[i];
    }

    return sqrt(sum);
}

// Function to calculate the error in a fit function result at a given x value
double calculateFitFunctionError(TF1* fitFunction, const std::vector<double>& coefficients, const std::vector<double>& coefficientErrors, double xValue) {
    // Set the coefficients with their respective errors
    for (int i = 0; i < coefficients.size(); ++i) {
        fitFunction->SetParameter(i, coefficients[i]);
        fitFunction->SetParError(i, coefficientErrors[i]);
    }

    // Calculate the derivative of the fit function with respect to each coefficient using a finite difference method
    const double epsilon = 1e-5;  // Small finite difference step
    std::vector<double> derivatives;

    for (int i = 0; i < coefficients.size(); ++i) {
        // Calculate the forward difference
        double originalValue = fitFunction->Eval(xValue);
        fitFunction->SetParameter(i, coefficients[i] + epsilon);
        double newValue = fitFunction->Eval(xValue);
        fitFunction->SetParameter(i, coefficients[i]);

        double derivative = (newValue - originalValue) / epsilon;
        derivatives.push_back(derivative);
    }

    // Calculate the error in the result using error propagation
    double error = 0.0;
    for (int i = 0; i < coefficients.size(); ++i) {
        error += TMath::Power(derivatives[i] * coefficientErrors[i], 2);
    }

    return TMath::Sqrt(error);
}

double calculateVectorStandardDeviation(const std::vector<double>& data) {
    double mean = 0.0;
    for (double value : data) {
        mean += value;
    }
    mean /= data.size();

    double variance = 0.0;
    for (double value : data) {
        variance += (value - mean) * (value - mean);
    }
    variance /= data.size();

    return std::sqrt(variance);
}

double calculateVectorStandardDeviationInRange(const std::vector<double>& data, int start_index, int end_index) {
    double mean = 0.0;
    for( int i = start_index; i <= end_index; i++){
        double value = data[i];
        mean += value;
    }
    double data_size = (end_index - start_index) + 1.0;
    mean /= data_size;

    double variance = 0.0;
    for( int i = start_index; i <= end_index; i++){
        double value = data[i];
        variance += (value - mean) * (value - mean);
    }
    variance /= (data_size);

    return std::sqrt(variance);
}

double calculateVectorMean(const std::vector<double>& data) {
    double sum = 0.0;
    for (double value : data) {
        sum += value;
    }
    return sum / data.size();
}
double calculateVectorMeanInRange(const std::vector<double>& data, int start_index, int end_index) {
    double sum = 0.0;
    double data_size = (end_index - start_index) + 1.0;
    for( int i = start_index; i <= end_index; i++){
        double value = data[i];
        sum += value;
    }
    return sum / (data_size);
}

double calculatePolynomialValue(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (int i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * TMath::Power(x, i);
    }
    return result;
}

// Function to calculate the error using error propagation
double calculatePolynomialError(const std::vector<double>& coefficients, const std::vector<double>& coefficientErrors, double x) {
    double error = 0.0;

    // Calculate the polynomial value at x
    double polynomialValue = calculatePolynomialValue(coefficients, x);

    for (int i = 1; i < coefficients.size(); ++i) {
        // Calculate the partial derivative of the polynomial with respect to coefficient a_i
        double derivative = TMath::Power(x, i);

        // Calculate the error contribution for this coefficient
        double errorContribution = TMath::Power(derivative * coefficientErrors[i], 2);

        error += errorContribution;
    }

    // Take the square root of the sum to get the total error
    return TMath::Sqrt(error);
}

//Lets make a general error calculation function for any polynomial evaluated at x. 
double CalculateFunctionErrorAtX( vector<double> coeff, vector<double> coeff_err, double x_eval) {
    // Create a TFormula object to define your fit function (Nth order polynomial)
    //We'll get the order of the polynomial from the size of coeff and coeff_err
    TString function_str = "[0]";
    int polN = coeff.size();

    for( int i = 1; i < polN; i++ ){
        function_str = Form("%s + [%i]", function_str.Data(), i );
        for( int j = 0; j < i; j++ ){
            function_str = Form("%s*x",function_str.Data());
        }
    }
    cout << "-----------------------------------------------------" << endl;
    cout << "Evaluating the pol" << coeff.size() << " function: " << function_str.Data() << endl;
    cout << "-----------------------------------------------------" << endl;

    //Define a formula from this
    TFormula fit_function("fit_function", function_str.Data() );
    // TFormula f("fit_function", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x");
   // Set the coefficients
    for( int i = 0; i < polN; i++ ){
        fit_function.SetParameter( i, coeff[i] );
    }

    // Calculate the value of the function at the point x
    double f_x = fit_function.Eval(x_eval);

    // Calculate the partial derivatives with respect to the coefficients
    // vector<double> df_dN = {0.0, 1.0}; //First derivative of the polN
    vector<double> df_dN = {}; //First derivative of the polN (dF(x)/d(A)) is 1 but we can calculate that on next line.
    for( int power = 0; power < polN; power++ ){
        df_dN.push_back( pow(x_eval, power) );
    }

    // Calculate the error using the error propagation formula for a polynomial expansion
    double sigma_f_x = 0;

    for( int i = 0; i < polN; i++){
        sigma_f_x += pow( df_dN[i]*coeff_err[i], 2 );
    }

    // Now, sigma_f_x contains the error at the point x
    //We need the sqrt though...
    cout << "Error at " << x_eval << " = " << sqrt(sigma_f_x) << endl;

    return sqrt(sigma_f_x);
}

// double calc_GMn_err( double CErr_tau_n, double CErr_eps_n, double CErr_scale_n, double CErr_scale_n_err, double CErr_scale_p, double CErr_scale_p_err, double CErr_sigma_n_MC, double CErr_sigma_n_MC_err, double CErr_GEn, double CErr_GEn_err, double return_err = 1 ){
//     double R_scale, R_scale_err, R_sigma_term, R_sigma_term_err;
//     double eps_GEn_term, eps_GEn_term_err;
//     double GMn_squared, GMn_err_squared, CErr_GMn, CErr_GMn_err;

//     double denominator;
//     double delta_R_term, delta_R_numer;
//     double delta_sig_term, delta_sig_numer;
//     double delta_G_term, delta_G_numer;

//     R_scale = CErr_scale_n/CErr_scale_p;
//     R_scale_err = CalculateErrorMultiplicationDivision( CErr_scale_n, CErr_scale_n_err, CErr_scale_p, CErr_scale_p_err, R_scale );

//     denominator = sqrt( (1/CErr_tau_n)*( (R_scale*CErr_sigma_n_MC) - (CErr_eps_n*CErr_GEn*CErr_GEn)));

//     delta_R_numer = CErr_sigma_n_MC*R_scale_err;
//     delta_R_term = delta_R_numer/denominator;

//     delta_sig_numer = R_scale*CErr_sigma_n_MC_err;
//     delta_sig_term = delta_sig_numer/denominator;

//     delta_G_numer = 2*CErr_eps_n*CErr_GEn*CErr_GEn_err;
//     delta_G_term = delta_G_numer/denominator;

//     GMn_err_squared = pow(delta_R_term, 2) + pow(delta_sig_term, 2) + pow(delta_G_term, 2);
//     CErr_GMn_err = sqrt( GMn_err_squared );

//     GMn_squared = (1/CErr_tau_n)*( (R_scale*CErr_sigma_n_MC) - (CErr_eps_n*CErr_GEn*CErr_GEn));
//     CErr_GMn = sqrt( GMn_squared );

//     if( return_err ){
//         return CErr_GMn_err;        
//     }
//     else{
//         return CErr_GMn;
//     }
// }

double calc_GMn_err( double CErr_tau_n, double CErr_eps_n, double CErr_scale_n, double CErr_scale_n_err, double CErr_scale_p, double CErr_scale_p_err, double CErr_sigma_n_MC, double CErr_sigma_n_MC_err, double CErr_GEn, double CErr_GEn_err, double return_err = 1 ){
    double F_scale, F_scale_err, R_sigma_term, R_sigma_term_err;
    double eps_GEn_term, eps_GEn_term_err;
    double GMn_squared, GMn_err_squared, CErr_GMn, CErr_GMn_err;

    double denominator;
    double delta_R_term, delta_R_numer;
    double delta_sig_term, delta_sig_numer;
    double delta_G_term, delta_G_numer;

    F_scale = CErr_scale_n/CErr_scale_p;
    F_scale_err = CalculateErrorMultiplicationDivision(CErr_scale_n, CErr_scale_n_err, CErr_scale_p, CErr_scale_p_err, F_scale);

    //Using formula for GMn-squared
    R_sigma_term_err = (F_scale*CErr_sigma_n_MC)*sqrt( pow( (F_scale_err/F_scale), 2) + pow( (CErr_sigma_n_MC_err/CErr_sigma_n_MC), 2));

    eps_GEn_term_err = 2.0*CErr_eps_n*CErr_GEn*CErr_GEn_err;

    CErr_GMn_err = R_sigma_term_err + eps_GEn_term_err;

    //Using formula for GMn-squared
    CErr_GMn = (1.0/CErr_tau_n)*( (F_scale*CErr_sigma_n_MC) - (CErr_eps_n*CErr_GEn*CErr_GEn) );

    if( return_err ){
        return sqrt(CErr_GMn_err);        
    }
    else{
        return sqrt(CErr_GMn);
    }
}

double calc_arr07_GMn_err( double CErr_R_np, double CErr_R_np_err, double CErr_F_np, double CErr_F_np_err, double CErr_sigma_R_p, double CErr_sigma_R_p_err, double CErr_tau ){
    // Statistical error only, so....
    // Considering the error in sigma_R_p from Arr07 to only come from systematics... they don't quote statistics.
    // So, we use sigma_R_p in the calculatino of GMn but not in the error propagation 
    CErr_sigma_R_p_err = 0.0;

    cout << "*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&" << endl;
    cout << "R_np: " << CErr_R_np << ", err: " << CErr_R_np_err << " (set to 0 for stat error)" << endl;
    cout << "F_np: " << CErr_F_np << ", err: " << CErr_F_np_err << endl;
    cout << "Sigma_R_p: " << CErr_sigma_R_p << ", err: " << CErr_sigma_R_p_err << endl;
    cout << "tau: " << CErr_tau << endl;

    double GMn_eval = sqrt( ((CErr_R_np*CErr_F_np*CErr_sigma_R_p)/(2*CErr_tau))*(1 + CErr_F_np) );

    double first_term = ((CErr_R_np*CErr_F_np*CErr_F_np*CErr_sigma_R_p)/(2*CErr_tau));

    double second_term = ((CErr_R_np*CErr_F_np*CErr_sigma_R_p)/(2*CErr_tau));

    double first_term_err = first_term*( sqrt( pow((CErr_R_np_err/CErr_R_np), 2) + 2*pow((CErr_F_np_err/CErr_F_np), 2) + pow((CErr_sigma_R_p_err/CErr_sigma_R_p), 2) ));

    double second_term_err = second_term*( sqrt( pow((CErr_R_np_err/CErr_R_np), 2) + pow((CErr_F_np_err/CErr_F_np), 2) + pow((CErr_sigma_R_p_err/CErr_sigma_R_p), 2) ));

    double addition_err = sqrt( pow(first_term_err, 2) + pow(second_term_err, 2) );

    addition_err = sqrt(addition_err);

    cout << "*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&" << endl;
    cout << "GMn_eval: " << GMn_eval << endl;
    cout << "First term: " << first_term << ", err: " << first_term_err  << endl;
    cout << "Second term: " << second_term << ", err: " << second_term_err << endl;
    cout << "addition_err: " << addition_err << endl;
    cout << "*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&" << endl;

    return addition_err;
}



double calc_sigma_R_n_bosted_Err(double CErr_eps_n, double CErr_GEn_bosted, double CErr_GEn_bosted_err, double CErr_tau_n, double CErr_GMn_bosted, double CErr_GMn_bosted_err, int return_err = 1 ){
    double GEn_term, GEn_term_err, GMn_term, GMn_term_err;
    double Sigma_R_n_bosted, Sigma_R_n_bosted_err;

    GEn_term = CErr_eps_n*CErr_GEn_bosted*CErr_GEn_bosted;
    GEn_term_err = CErr_eps_n*2*CErr_GEn_bosted*CErr_GEn_bosted_err;

    GMn_term = CErr_tau_n*CErr_GMn_bosted*CErr_GMn_bosted;
    GMn_term_err = CErr_tau_n*2*CErr_GMn_bosted*CErr_GMn_bosted_err;

    Sigma_R_n_bosted = GEn_term + GMn_term;
    Sigma_R_n_bosted_err = sqrt( pow(GEn_term_err, 2) + pow(GMn_term_err, 2) );

    if( return_err ){
        return Sigma_R_n_bosted_err;       
    }
    else{
        return Sigma_R_n_bosted;
    }
}

double calc_sigma_R_N_Err(double eps_N, double GEN, double GEN_err, double tau_N, double GMN, double GMN_err, int return_err = 1 ){
    double GEN_term, GEN_term_err, GMN_term, GMN_term_err;
    double Sigma_R_N, Sigma_R_N_err;

    GEN_term = eps_N*GEN*GEN;
    GEN_term_err = eps_N*2*GEN*GEN_err;

    GMN_term = tau_N*GMN*GMN;
    GMN_term_err = tau_N*2*GMN*GMN_err;

    Sigma_R_N = GEN_term + GMN_term;
    Sigma_R_N_err = sqrt( pow(GEN_term_err, 2) + pow(GMN_term_err, 2) );

    if( return_err ){
        return Sigma_R_N_err;       
    }
    else{
        return Sigma_R_N;
    }
}


double calc_arr07_sigma_R_n_err( double Fnp_ratio, double Fnp_ratio_err, double Rnp_ratio, double Rnp_ratio_err, double sigma_R_p, double sigma_R_p_err ){
    double eval = Fnp_ratio_err*Rnp_ratio*sigma_R_p;

    double Fnp_err = Fnp_ratio_err/Fnp_ratio;
    double Rnp_err = Rnp_ratio_err/Rnp_ratio;
    double sigma_R_err = sigma_R_p_err/sigma_R_p;

    double error = sqrt( pow(Fnp_err, 2) + pow(Rnp_err, 2) + pow(sigma_R_err, 2) ); 

    double final_err = error*eval;

    return final_err;
}

#endif
