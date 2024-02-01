#ifndef BLINDFACTOR_H
#define BLINDFACTOR_H
#include "TRandom.h"
#include "TMath.h"
#include <ctime>
#include <chrono>
#include <iostream>
#include <sstream>
#include <map>
#include <array>
#include "analysis_utility_functions.h"
map<char,int> myMap;
char letters[] = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',' '};

int datsize = sizeof(letters);

void makeMap(){
for(int i = 0;i < datsize ;i++){
int mynum = i % 10;
myMap[letters[i]]= mynum;
//cout << mynum << endl;
}
}
long int wordToInt(TString myword){

makeMap();

TString myint_string;

	for(int j=0;j<myword.Length();j++){
	char datLet = myword[j];
	//cout << datLet << endl;
	int datInt = myMap[datLet];
	
	//cout << datInt << endl;
	myint_string.Append(intToTString(datInt));
	//cout << myint_string[j] << endl;	
	}
//cout << myint_string << endl;
long int myint;
stringstream stream(myint_string.Data());
stream >> myint;

//cout << myint << endl;
return myint;
}
class blind_factor{
private:
//Define useful variable
double myblind_fac;

public:
blind_factor(TString datWord){
long int datNum = wordToInt(datWord);
//get the current time in seconds
std::time_t time = std::time(0);
//set the random number generator seed to the time you just found
gRandom->SetSeed(datNum);
//get the random number from the generator
double rand_num = gRandom->Uniform(0.88,1.12);
myblind_fac = rand_num;
}

//return a blind number. Will be blinded at yields
double blind(double mynum){
double blindNum = mynum*myblind_fac;
return blindNum;
}
//return a number that removes the blind factor
//should only be used to diagnostic check script and to actually remove number when physics analysis complete
double unblind(double datnum){
double unblindNum = datnum / myblind_fac;
return unblindNum;

}
//get the blind factor itself
//also for diagnositc check
double getBlindFac(){

return myblind_fac;
}
};

#endif
