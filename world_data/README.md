JB's World Data Plotting Mechanism

Developed by John Boyd (jab7bp@virginia.edu) circa October 2023

Primary use: to re-create world data plots as seen in papers and presentations.

Secondary use: to add personal data points to world data plots and to create my own/new versions of world data

Usage:

The main function to run this script is:

	plot_world_data.C

This looks for various scripts in the sub-directories of the main directory (world_data/)

This also pulls in a header file called jboyd_data_points which includes JB's experimental/calculated values for values such as GMn, nTPE, Rosenbluth Slope, etc.

Each sub-directory contains a header file which itself contains a plotting script for the given (named) world data plot.

The directories named after experiments (gmn, gmp, ntpe, etc.) contain files for the data/plots in the primary documents for the related experimental documents (proposals, etc.). 

The directory 'parameterizations' contains various scripts/files for plotting many different paramterizations for things like: GMn, GMp, GEn, GEp, GEn/GMn, GEp/GMp, etc. The specific parameterization method is in the naming of the file and possibly in the "intro" of that file. 

For example, parameterizations/kelly_param.h will contain the information and schema for plotting the Kelly parameterizations for the various form factors as seen in the G4SBS machinery.

To run:

From /world_data directory:

root -l ./plot_world_data.C

from here you can type: 

	plot_

and then press tab to populate a possible list of function names starting with 'plot_'

For example, for Kelly Paramterization type:

	plot_kelly_FFs(0) 

to plot all FFs. Or, use value1: 1 for GEp, 2 for GMp, 3 for GEn, 4 for GMn.

So, for GMn you can type:

	plot_kelly_FFs(4)

You can also add my personal data points to these world data plots by adding arguments to this:

	plot_kelly_FFs(4, 1)

This will plot the Kelly Parameterization for GMn along with jboyd's experimental values for GMn for SBS8 and SBS9 as defined i jboyd_data_points.h


