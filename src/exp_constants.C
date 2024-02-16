//exp_constants.C
//Author: Ezekiel Wertz
//The companion implementation for the header file
//Really this will just hold getter functions for any parameters that are meant to be constants but are either pass or kinematic dependent and not in config files

#include "../include/exp_constants.h"
namespace exp_constants{
  //m, height of the center of hcal above beam (m)
  const double getHCalOffset(TString myKin, TString pass){
  double hcal_offset;

  //First segment by pass. Later we may need to segment by kinematic
  
  	if(pass == "pass0" || pass == "pass1"){
	hcal_offset = -0.2897;
	}else if(pass == "pass2"){
	hcal_offset = 0.0;
	}else{
 	//We should never get here, cause then we have kinematic for data we dont have
 	cout << "Error: No hcal_offset for this pass: " << pass <<" Figure it out nerd!" << endl;
	}//end conditional
  return hcal_offset;
  }//end function

}//end namespace
