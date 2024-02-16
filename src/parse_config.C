//parse_config.C
//Author: Ezekiel Wertz
//The companion implementation for the header file
//implementation for parsing my config files. Should be able to handle all analysis types.

#include "../include/parse_config.h"
#include <iostream>
#include "TCut.h"
#include <vector>

  //constructor implementation
  parse_config::parse_config(const char *setup_file_name){
   ifstream setupfile(setup_file_name);
   	//check if there is a problem opening the file
        if(setupfile.fail()){
	 cout << "ERROR: There was a problem with the setup file " << setup_file_name << ". Figure it out nerd!" << endl;
        return;
	}
    TString myline;
    //Look for the input files. For MC we will have proton or neutron identifiers. If its data we will just have run numbers
    while(myline.ReadLine(setupfile) && !myline.BeginsWith("endrun") && !myline.BeginsWith("endfile")){  
	if(!myline.BeginsWith("#")){
	TObjArray *demObjs = myline.Tokenize(" ");
	int numObjs = demObjs->GetEntries();
        
		if(numObjs>1){
		TString temp = ((TObjString*) (*demObjs)[0])->GetString();
		bool myBool = utility::check_number(temp.Data());
		//cout << myBool << endl;
		//need to handle if we are dealing with MC or data
			//We found a number
			if(myBool){
				for(int i=0;i < numObjs;i++){
                        	//store the run numbers in the vector as ints
                        	int runnum_temp = (((TObjString*) (*demObjs)[i])->GetString()).Atoi();
                        	runnums.push_back(runnum_temp);
                        	//cout << runnum_temp << " ";        
                        	}//end for loop
			//We found a string
			}else{
			TString key = ((TObjString*) (*demObjs)[0])->GetString();
                	TString val = ((TObjString*) (*demObjs)[1])->GetString();
                        	if(key == "proton"){
                        	proton_root_file = val;
                        	//cout << "Proton File: " << proton_root_file << endl;
                        	}else if(key == "neutron"){
                        	neutron_root_file=val;
                        	//cout << "Neutron File: " << neutron_root_file << endl;
				}else if(key == "MC_file"){
                        	MC_file = val;
                        	//cout << "MC File: " << MC_file << endl;   
				}else{
                        	//We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
                        	cout << "Error:Found a key that this script can't handle. Fix that! "<< key << endl;
                        	return;
                        	}
			}//end conditional
    		}else{
    		//We either got an empty line or 1 element.
		cout << "Error:Line does not have the right number of elements. Look at the config file!" << endl;
    		return;
    		}//end conditional
	}//end conditional
    }//end of while loop
    //while loop for getting the rest of the parameters
    while(myline.ReadLine(setupfile)){
    	if(!myline.BeginsWith("#")){
        TObjArray *daObjs = myline.Tokenize(" ");
        int nObjs = daObjs->GetEntries();
        	if(nObjs >1){
                TString key = ((TObjString*) (*daObjs)[0])->GetString();
                TString val = ((TObjString*) (*daObjs)[1])->GetString();
			if(key == "exp"){
                        Exp = val;
                        //cout << "Experiment: " << Exp << endl;
                        }else if(key == "globalcut"){
			globalcut += val;
			cout << "Applying the following global cut to all data: " <<globalcut <<endl;
			}else if(key == "kin"){
                        kin = val;
                        //cout << "Kinematic " << kin << endl;
                        }else if(key == "data_map_name"){
                        data_file_name = val;
                        //cout << "Data File " << data_fiel_name << endl;
                        }else if(key == "kinematic_name"){
                        kinematic_file_name = val;
                        //cout << "Kinematic File " << kinematic_file_name << endl;
                        }else if(key == "pass"){
                        pass = val;
			//cout << "Pass " << pass << endl;
			}else if(key == "SBS_field"){
                        SBS_field = val.Atoi();
                        //cout << "SBS Field " << SBS_field << endl;
                        }else if(key == "W2_mean"){
                        W2_mean = val.Atof();
                        //cout << "W2 mean " << W2_mean << endl;
                        }else if(key == "W2_sigma"){
                        W2_sigma = val.Atof();
                        //cout << "W2 sig " << W2_sig << endl;
                        }else if(key == "W2_sigfac"){
                        W2_sigfac = val.Atof();
                        //cout << "W2 sig fac " << W2_sigfac << endl;
                        }else if(key == "targ"){
                        targ = val;
                        //cout << "Target " << targ << endl;
                        } else if(key == "MAXNTRACKS"){
                        MAXNTRACKS = val.Atoi();
                        //cout << "Max Number of Tracks per event " << MAXNTRACKS << endl;
                        }else if(key == "dxO_n"){
                        dxO_n = val.Atof();
                        //cout << "x-position of neutron spot " << dxO_n << endl;
                        } else if(key == "dyO_n"){
                        dyO_n = val.Atof();
                        //cout << "y-position of neutron spot " << dyO_n << endl;
                        }else if(key == "dxsig_n"){
                        dxsig_n = val.Atof();
                        //cout << "x sigma of neutron spot " << dxsig_n << endl;
                        }else if(key == "dysig_n"){
                        dysig_n = val.Atof();
                        //cout << "y sigma of neutron spot " << dysig_n << endl;
                        }else if(key == "dxO_p"){
                        dxO_p = val.Atof();
                        //cout << "x-position of proton spot " << dxO_p << endl;
                        }else if(key == "dyO_p"){
                        dyO_p = val.Atof();
                        //cout << "y-position of proton spot " << dyO_p << endl;
                        }else if(key == "dxsig_p"){
                        dxsig_p = val.Atof();
                        //cout << "x sigma of proton spot " << dxsig_p << endl;
                        }else if(key == "dysig_p"){
                        dysig_p = val.Atof();
                        //cout << "y sigma of proton spot " << dysig_p << endl;
                        }else if(key == "dx_pn"){
                        dx_pn = val.Atof();
                        //cout << "max x difference between peaks " << dx_pn << endl;
                        }else if(key == "useAlshield"){
                        useAlshield = val.Atoi();
                        //cout << "Use Al shield " << useAlshield << endl;
                        }else if(key == "dx_low"){
                        dx_low = val.Atof();
                        //cout << "dx plot lower bound " << dx_low << endl;
                        }else if(key == "dx_high"){
                        dx_high = val.Atof();
                        //cout << "dx plot higher bound " << dx_high << endl;
                        }else if(key == "dy_low"){
                        dy_low = val.Atof();
                        //cout << "dy plot lower bound " << dy_low << endl;
                        }else if(key == "dy_high"){
                        dy_high = val.Atof();
                        //cout << "dy plot higher bound " << dy_high << endl;
                        }else if(key == "dxsig_n_fac"){
                        dxsig_n_fac = val.Atof();
                        //cout << "dx sigma factor for neutron " << dxsig_n_fac << endl;
                        }else if(key == "dxsig_p_fac"){
                        dxsig_p_fac = val.Atof();
                        //cout << "dx sigma factor for proton " << dxsig_p_fac << endl;
                        }else if(key == "dysig_n_fac"){
                        dysig_n_fac = val.Atof();
                        //cout << "dy sigma factor for neutron " << dysig_n_fac << endl;
                        }else if(key == "dysig_p_fac"){
                        dysig_p_fac = val.Atof();
                        //cout << "dy sigma factor for proton " << dysig_p_fac << endl;
                        }else if(key == "proton_thresh_fac"){
                        proton_thresh_fac = val.Atoi();
                        //cout << "Proton thresh factor: " << proton_thresh_fac << endl;
                        }else if(key == "neutron_thresh_fac"){
                        neutron_thresh_fac = val.Atoi();
                        //cout << "Neutron thresh factor: " << neutron_thresh_fac << endl;
                        }else if(key == "pmin"){
                        pmin = val.Atof();
                        //cout << "pmin: " << pmin << endl;
                        }else if(key == "pmax"){
                        pmax = val.Atof();
                        //cout << "pmax: " << pmax << endl;
                        }else if(key == "Emin"){
                        Emin = val.Atof();
                        //cout << "Emin: " << Emin << endl;
                        }else if(key == "Emax"){
                        Emax = val.Atof();
                        //cout << "Emax: " << Emax << endl;
                        }else if(key == "num_bin"){
                        num_bin = val.Atoi();
                        //cout << "Num bins: " << num_bin << endl;
                        }else if(key == "coin_mean"){
                        coin_mean = val.Atof();
                        //cout << "coin mean " << coin_mean << endl;
                        }else if(key == "coin_sigma"){
                        coin_sigma = val.Atof();
                        //cout << "coin sigma " << coin_sigma << endl;
                        }else if(key == "coin_sig_fac"){
                        coin_sig_fac = val.Atof();
                        //cout << "coin sig fac " << coin_sig_fac << endl;
                        }else if(key == "binfac"){
                        binfac = val.Atof();
                        //cout << "binfac " << binfac << endl;
                        }else if(key == "W2fitmax"){
                        W2fitmax = val.Atof();
                        //cout << "W2fitmax " << W2fitmax << endl;
                        }else if(key == "W2fitmaxwide"){
                        W2fitmaxwide = val.Atof();
                        //cout << "W2fitmaxwide " << W2fitmaxwide << endl;
                        }else if(key == "hcalemin"){
                        hcalemin = val.Atof();
                        //cout << "hcalemin " << hcalemin << endl;
                        }else if(key == "thetapq_low"){
                        thetapq_low = val.Atof();
                        //cout << "thetapq_low " << thetapq_low << endl;
                        }else if(key == "thetapq_high"){
                        thetapq_high = val.Atof();
                        //cout << "thetapq_high " << thetapq_high << endl;
                        }else if(key == "sf"){
                        sf = val.Atof();
                        //cout << "Scale Field " << sf << endl;
                        }                        
                	else{
			//We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
			cout << "Error:Found a key that this script can't handle. Fix that! "<< key << endl;
                        return;
                	}
        	//remove the objects to ensure a new comes through and no runaway
        	delete daObjs;
    		}//end conditional
    		else{
   		//We either got an empty line or 1 element.
    		cout << "Line does not have the right number of elements. Look at the config file!" << endl;
    		return;
    		}//end conditinal
    	}//end conditional
    }//end while
    if(runnums.empty() && ((proton_root_file.Length() == 0) &&  (neutron_root_file.Length() ==0))){
    // if there are data or MC files we should return true and therefore throw an error
    cout << "Error: No data files in the config file, I can't do anything if I don't know where the data lives!" << endl;
        return;
    }
  }//end constructor


  //Implement get functions for each variable

  TString parse_config::getExp(){ return Exp; }

  TString parse_config::getKin(){ return kin; }

  TString parse_config::getKinFileName(){ return kinematic_file_name; }

  TString parse_config::getDataFileName(){ return data_file_name; }

  TString parse_config::getTarg(){ return targ; }

  TString parse_config::getProtonFileName(){ return proton_root_file; }

  TString parse_config::getNeutronFileName(){ return neutron_root_file; }

  TString parse_config::getMCFileName(){ return MC_file; }

  TString parse_config::getPass(){ return pass; }

  int parse_config::getSBSField(){ return SBS_field; }

  int parse_config::getAlshield(){ return useAlshield; }

  int parse_config::getMAXNTRACKS(){ return MAXNTRACKS; }

  double parse_config::get_dxOn(){ return dxO_n; }

  double parse_config::get_dyOn(){ return dyO_n; }

  double parse_config::get_dxsign(){ return dxsig_n; }

  double parse_config::get_dysign(){ return dysig_n; }

  double parse_config::get_dxOp(){ return dxO_p; }

  double parse_config::get_dyOp(){ return dyO_p; }

  double parse_config::get_dxsigp(){ return dxsig_p; }

  double parse_config::get_dysigp(){ return dysig_p; }

  double parse_config::get_dxpn(){ return dx_pn; }

  double parse_config::getW2Mean(){ return W2_mean; }

  double parse_config::getW2Sigma(){ return W2_sigma; }

  double parse_config::getW2SigFac(){ return W2_sigfac; }

  double parse_config::get_dxLow(){ return dx_low; }

  double parse_config::get_dxHigh(){ return dx_high; }

  double parse_config::get_dyLow(){ return dy_low; }

  double parse_config::get_dyHigh(){ return dy_high; }

  double parse_config::get_dxSignFac(){ return dxsig_n_fac; }

  double parse_config::get_dxSigpFac(){ return dxsig_p_fac; }

  double parse_config::get_dySignFac(){ return dysig_n_fac; }

  double parse_config::get_dySigpFac(){ return dysig_p_fac; }

  double parse_config::getCoinMean(){ return coin_mean; }
 
  double parse_config::getCoinSig(){ return coin_sigma; }

  double parse_config::getCoinSigFac(){ return coin_sig_fac; }

  double parse_config::getHCaleMin(){ return hcalemin; }

  double parse_config::getProtonThreshFac(){ return proton_thresh_fac; }

  double parse_config::getNeutronThreshFac(){ return neutron_thresh_fac; }

  double parse_config::getNumBin(){ return num_bin; }

  double parse_config::getPmax(){ return pmax; }

  double parse_config::getPmin(){ return pmin; }

  double parse_config::getEmin(){ return Emin; }

  double parse_config::getEmax(){ return Emax; }

  double parse_config::getThetapqLow(){ return thetapq_low; }

  double parse_config::getThetapqHigh(){ return thetapq_high; }

  double parse_config::getW2FitMax(){ return W2fitmax; }

  double parse_config::getW2FitMaxWide(){ return W2fitmaxwide; }

  double parse_config::getBinFac(){ return binfac; }

  double parse_config::get_sf(){ return sf; }

  TCut parse_config::getGlobalCut(){ return globalcut; }

  vector<int> parse_config::getRunNums(){ return runnums; }

  void parse_config::printRunNums(){
  	for(int j=0; j< runnums.size(); j++){
		if(j==runnums.size()-1){
        	cout << runnums[j] << endl;
        	}else{ 
		cout << runnums[j] << ", ";
		}
  	}
  }

  //only work on LD2 data
  void parse_config::printDataYields(){
  cout << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl
       << Form("Global Cuts: %s",getGlobalCut().GetTitle())										  << endl
       << "RunNums: "                                                    							   	         ;
       printRunNums();
  cout << Form("Experiment: %s",(getExp()).Data())              									  << endl
       << Form("Kinematic: %s",(getKin()).Data())        										  << endl
       << Form("Data File Map: %s", getDataFileName().Data())										  << endl
       << Form("Kinematic File Name: %s",getKinFileName().Data())									  << endl
       << Form("SBS Field: %i",getSBSField())									                          << endl
       << Form("W2 Mean: %f",getW2Mean())												  << endl
       << Form("W2 Sigma: %f",getW2Sigma())												  << endl
       << Form("W2 Sigma Fac: %f",getW2SigFac())											  << endl
       << Form("Target: %s",getTarg().Data())												  << endl
       << Form("Maxntracks: %i",getMAXNTRACKS())											  << endl
       << Form("dxO_n: %f",get_dxOn())													  << endl
       << Form("dyO_n: %f",get_dyOn())                                                                                                    << endl
       << Form("dxsign_n: %f",get_dxsign())												  << endl
       << Form("dysign_n: %f",get_dysign())                                                                                               << endl
       << Form("dxO_p: %f",get_dxOp())                                                                                                    << endl
       << Form("dyO_p: %f",get_dyOp())                                                                                                    << endl
       << Form("dxsign_p: %f",get_dxsigp())                                                                                               << endl
       << Form("dysign_p: %f",get_dysigp())                                                                                               << endl
       << Form("dx_pn: %f",get_dxpn())													  << endl
       << Form("dx low: %f",get_dxLow())												  << endl
       << Form("dx high: %f",get_dxHigh())												  << endl
       << Form("dy low: %f",get_dyLow())                                                                                                  << endl
       << Form("dy high: %f",get_dyHigh())                                                                                                << endl
       << Form("use Al shield: %i",getAlshield())											  << endl
       << Form("dx sig n fac: %f", get_dxSignFac())											  << endl
       << Form("dx sig p fac: %f", get_dxSigpFac())                                                                                       << endl
       << Form("dy sig n fac: %f", get_dySignFac())                                                                                       << endl
       << Form("dy sig p fac: %f", get_dySigpFac())                                                                                       << endl
       << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl;
  }

  //only work on HCal Efficiency LH2 data
  //need to test
  void parse_config::printDataHCalEff(){
  cout << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl
       << Form("Global Cuts: %s",getGlobalCut().GetTitle())                                                                               << endl
       << "RunNums: "                                                                                                                            ;
       printRunNums();
  cout << Form("Experiment: %s",(getExp()).Data())                                                                                        << endl
       << Form("Kinematic: %s",(getKin()).Data())                                                                                         << endl
       << Form("Data File Map: %s", getDataFileName().Data())                                                                             << endl
       << Form("Kinematic File Name: %s",getKinFileName().Data())                                                                         << endl
       << Form("SBS Field: %i",getSBSField())                                                                                             << endl
       << Form("W2 Mean: %f",getW2Mean())                                                                                                 << endl
       << Form("W2 Sigma: %f",getW2Sigma())                                                                                               << endl
       << Form("W2 Sigma Fac: %f",getW2SigFac())                                                                                          << endl
       << Form("Target: %s",getTarg().Data())                                                                                             << endl
       << Form("Coin Mean: %f",getCoinMean())												  << endl
       << Form("Coin Sigma: %f",getCoinSig())												  << endl
       << Form("Coin Sig Fac: %f",getCoinSigFac())											  << endl
       << Form("Maxntracks: %i",getMAXNTRACKS())                                                                                          << endl
       << Form("dxO_n: %f",get_dxOn())                                                                                                    << endl
       << Form("dyO_n: %f",get_dyOn())                                                                                                    << endl
       << Form("dxsign_n: %f",get_dxsign())                                                                                               << endl
       << Form("dysign_n: %f",get_dysign())                                                                                               << endl
       << Form("dxO_p: %f",get_dxOp())                                                                                                    << endl
       << Form("dyO_p: %f",get_dyOp())                                                                                                    << endl
       << Form("dxsign_p: %f",get_dxsigp())                                                                                               << endl
       << Form("dysign_p: %f",get_dysigp())                                                                                               << endl
       << Form("dx_pn: %f",get_dxpn())                                                                                                    << endl
       << Form("dx low: %f",get_dxLow())                                                                                                  << endl
       << Form("dx high: %f",get_dxHigh())                                                                                                << endl
       << Form("dy low: %f",get_dyLow())                                                                                                  << endl
       << Form("dy high: %f",get_dyHigh())                                                                                                << endl
       << Form("use Al shield: %i",getAlshield())                                                                                         << endl
       << Form("Bin Fac: %f",getBinFac())												  << endl
       << Form("W2 FitMax: %f",getW2FitMax())												  << endl
       << Form("W2 FitMax Wide: %f",getW2FitMaxWide())											  << endl
       << Form("HCal Emin: %f",getHCaleMin())												  << endl
       << Form("Thetapq Low: %f",getThetapqLow())											  << endl
       << Form("Thetapq High: %f",getThetapqHigh())   
       << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl;  

  }

  //only work on MC LD2
  void parse_config::printMCYields(){
  //TODO, will do once I have such config files to begin with
  }

  //only work on HCal Efficiency MC LH2
  //need to test 
  void parse_config::printMCHCalEff(){
  cout << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl
       << Form("Proton File: %s",getProtonFileName().Data())										  << endl
       << Form("Neutron File: %s",getNeutronFileName().Data())                                                                            << endl
       << Form("Experiment: %s",(getExp()).Data())                                                                                        << endl
       << Form("Kinematic: %s",(getKin()).Data())                                                                                         << endl
       << Form("Kinematic File Name: %s",getKinFileName().Data())                                                                         << endl
       << Form("Proton Thresh Fac: %f",getProtonThreshFac())										  << endl
       << Form("Neutron Thresh Fac: %f",getNeutronThreshFac())                                                                            << endl
       << Form("Pmin: %f",getPmin())													  << endl
       << Form("Pmax: %f",getPmax())	   	                                                                                          << endl
       << Form("Emin: %f",getEmin())   	                                                                                           	  << endl
       << Form("Emax: %f",getEmax())                                                                                             	  << endl
       << Form("Num Bin: %f",getNumBin())												  << endl
       << Form("dx low: %f",get_dxLow())                                                                                                  << endl
       << Form("dx high: %f",get_dxHigh())                                                                                                << endl
       << Form("dy low: %f",get_dyLow())                                                                                                  << endl
       << Form("dy high: %f",get_dyHigh())                                                                                                << endl
       << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl;
  }
  