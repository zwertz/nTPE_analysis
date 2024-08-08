#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <TF1.h>
#include "TStopwatch.h"

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/ewertz/nTPE_analysis/include/physics_constants.h"
#include "/w/halla-scshelf2102/sbs/ewertz/nTPE_analysis/include/calc_errors.h"

//World data files:
	#include "./jboyd_data_points.h"
	//ntpe
#include "./ntpe/christy_ntpe_gmp.h"
#include "./ntpe/ntpe_proposal_fig2.h"
#include "./ntpe/ntpe_proposal_fig2_LT_PT_overlay.h"
#include "./ntpe/ntpe_proposal_fig2_V2.h"
#include "./ntpe/ntpe_RS.h"
#include "./ntpe/rosenbluth_sep_ntpe.h"

	//gmn
#include "./gmn/gmn_proposal_data.h"

	//gmp
#include "./gmp/gmp_pr12_07_108.h"

	//parameterizations
#include "./parameterizations/ye_arrington.h"
#include "./parameterizations/christy_ffr_param.h"
#include "./parameterizations/kelly_param.h"
#include "./parameterizations/bosted_param.h"
#include "./parameterizations/bradford_param.h"

Int_t n_eff_proton_rad;

TGraph *TG = new TGraph();

void plot_world_data(){

}
