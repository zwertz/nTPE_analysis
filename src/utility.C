#include "../include/utility.h"
#include <iostream>
#include <string>
#include "TMath.h"

//Author Ezekiel Wertz
//Implementation for useful task functions

namespace utility{

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

} //end namespace
