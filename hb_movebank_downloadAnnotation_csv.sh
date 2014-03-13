#!/bin/bash

# Purpose: This shell script uses Rcurl to retrieve annotation requests from movebank.
# The arguments required include: ak_lut=csv lut file generated from hb-submit_movebank_request.sh; 
# url=url where files can be retrieved - currently hardcoded; outdir=output directory for downloaded files. 
#
# Author: Tina Cormier
# Date: March 12, 2014
# Status: Runs
#
if [ "$#" -ne 2 ]; then
    echo "usage: hb-submit_movebank_request.sh <url> <xy> <xml> <un> <pw> <outdir>"
    echo
    echo "example: hb-submit_movebank_request.sh http://www.bioinfo.mpg.de/orn-gateway/request-annotation-xml.jsp /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/test_calliope.csv /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/NDVI_EVI.xml tcormier nohg3ITh /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests/"
    exit
fi

ak_lut=${1}
outdir=${2}

###############################
# #hardcoded for testing:
ak_lut=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests_accessKeys_lut.csv
outdir=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/downloads/
###############################

#cat the file and check for available downloads

cat ${ak_lut} | while read line
do
	#if the accesskey is listed as "requested," check to see if it's available yet
	ak_status=$(echo $line | cut -d',' -f5)
	if [ ${ak_status} = requested ]
	then
		echo "status = requested: checking download availability now. . ."
		#get accesskey
		ak=$(echo $line | cut -d',' -f1)
		url=http://www.bioinfo.mpg.de/orn-gateway/status.jsp?access-key=${ak}
		#get status page
		stPage=$(curl ${url} | grep "Status" | awk -F ';' '{print $2}')
		if [ stPage == "Available" ]
		then
			track=$(basename $(echo ${line} | cut -d',' -f3) .csv)
			vars=$(basename $(echo ${line} | cut -d',' -f4) .xml)
			today=$(echo ${line} | cut -d',' -f2)
			outfile=${outdir}${track}_${vars}_${today}_${ak}.csv
			dl_url=http://www.bioinfo.mpg.de/orn-gateway/download.jsp?access-key=${ak}
			curl -o ${outfile} ${dl_url}
		else
			echo "${ak} still pending. . .check back later"
		fi #end status available if
	fi #end status requested if
