#include "../include/utility.h"
#include <iostream>
#include <filesystem>
#include <string>
#include "TMath.h"
#include <regex>


//Author Ezekiel Wertz
//Implementation for useful task functions

namespace utility{
  //function to determine if what is found was a number or a string
  bool check_number(const char *myChar){
   bool mybool=false;
  	for(int i=0; i< strlen(myChar); i++){
		//So we ever find not a digit. Immediately return mybool=false. So we have a string
		//cout << mybool << endl;
		if(isdigit(myChar[i])==false){
		return mybool;
		}
  	}
  //If we get through the for loop we never found not a digit. Make it true and return
  mybool=true;
  return mybool;
  }//end check_number

  //convert from degrees to radians
  double DegToRad(double myNum){
  return (myNum * (TMath::Pi()/180.0));
  }
  //convert from radians to degrees
  double RadToDeg(double datRad){
  return ((datRad * 180.0)/(TMath::Pi()));
  }

  //convert an int to a TString. Just a helper function
  TString intToTString(int datInt){
  stringstream ss;
  ss << datInt;
  std::string datInt_temp;
  ss >> datInt_temp;
  TString datInt_string(datInt_temp.c_str());
  //cout << datInt_string << endl;
  return datInt_string;
  }

  //This is a helper function that will tell you the default directory for all analysis output
  TString getOutputDir(){
  return output_directory;
  }

  //Helper function to make output file name for parsed data root file
  TString makeOutputFileNameParse(TString exp,TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/Zeke_parsed_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function to make output file name for parse mc root file
  TString makeOutputFileName_MCParse(TString exp, TString Kin, int SBS_field){
  TString outfile = Form("%s/Zeke_MC_pn_parsed_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),Kin.Data(),SBS_field);
  return outfile;
  }


  //Helper function to make output file name for yield or ratio information
  TString makeYieldReportFileName(TString exp,TString pass,  TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/yieldRep_Zeke_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(), pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }
  
  //Helper function to make output file name for HCal Efficiency from data
  TString makeOutputFileNameHCalEffParse(TString exp,TString pass, TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/Zeke_HCalEff_parsed_%s_%s_%s_%s_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function to make output file name for HCal Efficiency from MC
  TString makeOutputFileName_HCalEffMC(TString exp, TString Kin){
  TString outfile = Form("%s/Zeke_HCalEff_MC_%s_%s",(getOutputDir()).Data(),exp.Data(),Kin.Data());
 
  return outfile;
  }

  //TODO: Will need a function to make file name for compare HCal Efficiency data with MC

  //Helper function to make output file for HCal Efficiency Bootstrap
  TString makeOutputFileName_HCalBS(TString exp,TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event){
  TString outfile = Form("%s/Zeke_HCalEff_bootstrap_%s_%s_%s_%s_%i_event_%i_%i.root",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,first_event,last_event);
  //cout << outfile << endl;
  return outfile;
  }
  //Helper function for making file name for HCal Efficiency report
  TString makeReportFileName_HCalEff(TString exp,TString pass,  TString Kin, int SBS_field,TString target){
  TString outfile = Form("%s/efficiencyRep_Zeke_%s_%s_%s_%s_%i.txt",(getOutputDir()).Data(),exp.Data(), pass.Data(),Kin.Data(),target.Data(),SBS_field);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function for making file name for HCal Efficiency report with boot strap
  TString makeReportFileName_HCalEffBS(TString exp,TString pass, TString Kin, int SBS_field,TString target,int first_event, int last_event){
  TString outfile = Form("%s/efficiencyRep_Zeke_%s_%s_%s_%s_%i_event_%i_%i.txt",(getOutputDir()).Data(),exp.Data(),pass.Data(),Kin.Data(),target.Data(),SBS_field,first_event,last_event);
  //cout << outfile << endl;
  return outfile;
  }

  //Helper function for output file name for scale field study for MC
  TString makeOutputFileName_MCsf(TString exp, TString Kin, int SBS_field,double sf,TString target){
  TString outfile = Form("%s/MC_Zeke_%s_%s_%s_%i_sf%.1f.root",(getOutputDir()).Data(),exp.Data(),Kin.Data(),target.Data(),SBS_field,sf);
  //cout << outfile << endl;
   return outfile;
  }

  //Helper function to find MC histogram files. Right now this only supports John Boyd simulation files. But in future might support others.
  vector<string> findHistFiles(TString replay_type,TString histDirectory,TString partialName){
  vector<string> myFiles;
  const string& histDir_str = histDirectory.Data();  
  const string& partialName_str = partialName.Data();
	if(replay_type == "jboyd"){
		//populate the vector with all matching hist files
		for(const auto& entry : filesystem::directory_iterator(histDir_str)){
			//check that the file is a regular file and that the partial name is found somewhere in the filename
			if(entry.is_regular_file() && (entry.path().filename().string().find(partialName_str) != string::npos) ){
			//if we find one that should be a match put it on the vector
			myFiles.push_back(entry.path().string());
			}
		}
	}else{
	cout << "Error: MC file type not supported. Replay Type: " << replay_type << " This needs fixing!" << endl;
	}
  return myFiles;
  }
 
  //Helper function to match all the MC root files with the already found MC hist files. This function is not ideal, it modifies both of the vectors by the time the function is done. Though it does not have a return type
  void matchMCFiles(TString replay_type,vector<string>& histFiles,vector<string>& rootFiles, TString rootDirectory){
	//vector to track file indices that don't match
	vector<int> unmatchedIndices;
	int matches = 0;  
	const string& rootDir_str = rootDirectory.Data();
	if(replay_type == "jboyd"){
		//loop over the values in the hist file vector
		for(size_t i = 0; i < histFiles.size(); ++i){
		cout << "Looking for hist to root dir  matches, " << i << "/" << histFiles.size() << " total matches " << matches << "\r";
       		cout.flush();
		filesystem::path histPath = histFiles[i];
		string hist_filename = histPath.filename().string();

		//modfiy the filename for the expect root file name
		string expect_root_filename = "replayed_digitized_" + hist_filename;
		expect_root_filename = expect_root_filename.substr(0, expect_root_filename.find_last_of('.')) + ".root";

		//boolean to track if we find match
		bool matchFound = false;

			//iterate over the root file directory and see if there is a match
			for(const auto& entry : filesystem::directory_iterator(rootDir_str)){
				//check if there is a match
				if(entry.is_regular_file() && (entry.path().filename() == expect_root_filename)){
				//if we did put the file name on the root file vector and update other parameters
				rootFiles.push_back(entry.path().string());
				matchFound = true;
				matches++;
				break; //Since we found the match we can stop searching for this iteration
				}//conditional
			}//inner for loop
		
			if(!matchFound){
			//if match not found indicate that index needs removal
			unmatchedIndices.push_back(i);
			}
		}// end outer for loop

		//Remove any unmatched indices
		//Need to iterate in reverse because of implementation of vectors and to get correct indices
		for(auto j = unmatchedIndices.rbegin(); j != unmatchedIndices.rend() ; ++j){

		histFiles.erase(histFiles.begin() + *j);
		}

	}else{
        cout << "Error: MC file type not supported during matching. Replay Type: " << replay_type << " This needs fixing!" << endl;
        }
  //There is no return here so this is a little confusing but when the function terminates the root file vector should be populate and the non matching hist vector elements should be removed. So this function modifies both vectors to have matching information, since we are sending references to the original vector.
  }//end matching function

  //Function that searches through two vectors of strings and removes entries without a matching pair
  //since we are sending references to vectors it should modifiy the original reference
  void syncJobNumbers(vector<string>& proton_vec,vector<string>& neutron_vec){
  regex jobNumberRegex("job([0-9]+)");
  std::smatch match;
  std::unordered_set<string> jobNumbers_pro, jobNumbers_neu;

  	//Get job numbers for proton files
  	for(const auto& str: proton_vec){
  		//search for files with the job numbers
		if(std::regex_search(str,match,jobNumberRegex)){
		//store references if we find one
		jobNumbers_pro.insert(match[0]);
		}
  	}
  	//Get job numbers for nuetron files
  	for(const auto& str: neutron_vec){
        	//search for files with the job numbers
        	if(std::regex_search(str,match,jobNumberRegex)){
        	//store references if we find one
        	jobNumbers_neu.insert(match[0]);
        	}
  	}
	//remove unpaired entries from proton vec
	//this codes a little odd, unfamiliar library.
	//I suspect it is doing something like searching from the begin to end of the vector and removing the string. But the string is made on the fly and is some how checking if one string found in the first vector is also in the second vector
	proton_vec.erase(std::remove_if(proton_vec.begin(),proton_vec.end(),[&](const string& str){
			if(std::regex_search(str,match,jobNumberRegex) && jobNumbers_neu.find(match[0]) == jobNumbers_neu.end() ){
			return true;
			}
		return false;
		}),proton_vec.end());
	//remove unpaired entries from neutron vec
        neutron_vec.erase(std::remove_if(neutron_vec.begin(),neutron_vec.end(),[&](const string& str){
                        if(std::regex_search(str,match,jobNumberRegex) && jobNumbers_pro.find(match[0]) == jobNumbers_pro.end() ){
                        return true;
                        }
                return false;
                }),neutron_vec.end());

  }// end sync job function

  //function to parse simc hist files replayed by jboyd. Need to verify that this is still valid.
  double searchSimcHistFile(const TString &pattern, const TString &filename){
	ifstream inputFile(filename.Data());//Try to open the file with the corresponding name of interest

  	//make sure the file exists. If it does not return
	if(inputFile.fail()){
	cerr << "Error: There was a problem with the input hist file!" << endl;
	cerr << "Filename: " << filename << endl;
	return -99; //This value does not matter as long as it is nonsensical for the other parameters.
	}

  TString myLine;
  double foundValue = -1; //Default to a value that would not make sense physically for any parameter
  TString previousLine; //Keep track of the previous line for exceptions
  	//search through each line of the file for the designated pattern
	while(myLine.ReadLine(inputFile)){
		//If we found the pattern in the line
		if(myLine.Contains(pattern)){
		//handle the case of luminosity or other patterns with scientific notation
			if(myLine.EndsWith("ub^-1") || myLine.EndsWith("geV^2")){
			//Extract the numerical value
			TString delim("=");
			TObjArray* tokens = myLine.Tokenize(delim);
				//check that we found enough elements
				if(tokens->GetEntries() >= 2){
				TObjString* value = dynamic_cast<TObjString*>(tokens->At(1)); 
					if(value){
					TString valueString = value->String();
					//convert to double
					foundValue = valueString.Atof();
					delete tokens; //cleanup
					return foundValue;
					}
				}
			delete tokens; //cleanup if for some reason we don't find a value this way
			}else{
			//For other nonspecific patterns
			TString valueString;
			TString delim("=");
			TObjArray* tokens = myLine.Tokenize(delim);
				if(tokens->GetEntries() >= 2){
				TObjString* value = dynamic_cast<TObjString*>(tokens->At(1));
					if(value){
					valueString = value->String();
					}
				}
				if(valueString.IsFloat()){
				//convert to double
				foundValue = valueString.Atof();
                                delete tokens; //cleanup
                                return foundValue;
				}
			delete tokens; //cleanup if for some reason we don't find a value this way
			}
		}//end conditional on pattern line check
	previousLine = myLine;
	}//end of while loop
  //Handle luminosity and other execptions when an appropriate line is not found
  	if((pattern == "luminosity" && previousLine.EndsWith("ub^-1")) || (pattern == "genvol")){
	TString delim("=");
	TObjArray* tokens = myLine.Tokenize(delim);
		if(tokens->GetEntries() >= 2){
                TObjString* value = dynamic_cast<TObjString*>(tokens->At(1));
			if(value){
                        TString valueString = value->String();
                       	foundValue = valueString.Atof();
			}

		}
	delete tokens; //cleanup
	}
  inputFile.close();

  return foundValue;
  }//end of hist search function
  


} //end namespace
