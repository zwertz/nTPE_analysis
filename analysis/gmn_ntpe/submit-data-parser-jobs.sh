#!/bin/bash

# swif2 job submit for data_elastic_parse.C for gmn analysis

#Usage
#./submit-data-parse-jobs.sh <string: config file name> 

file=$1

work_flow='elastic_data_parse'

swif2 create $work_flow

script='/w/halla-scshelf2102/sbs/ewertz/nTPE_analysis/analysis/gmn_ntpe/run-data-parser.sh'

cd /w/halla-scshelf2102/sbs/ewertz/nTPE_analysis/analysis/gmn_ntpe

echo -e "\n Submitting " $script $file "\n"

swif2 add-job -workflow $work_flow -partition production -name elastic_data_parse -cores 1 -disk 10GB -ram 2500MB $script $file

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow




