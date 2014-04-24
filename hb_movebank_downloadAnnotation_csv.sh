#!/bin/bash

# Purpose: This shell script uses Rcurl to retrieve annotation requests from movebank.
# The arguments required include: ak_lut=csv lut file generated from hb-submit_movebank_request.sh; 
# url=url where files can be retrieved - currently hardcoded; outdir=output directory for downloaded files. 
#
# Author: Tina Cormier
# Date: March 12, 2014
# Status: Runs
#
if [ "$#" -ne 3 ]; then
    echo "usage: hb_movebank_downloadAnnotation_csv.sh <ak_lut> <outdir> <logfile>"
    echo
    echo "example: hb_movebank_downloadAnnotation_csv.sh /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests_accessKeys_lut.csv /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/downloads/, /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/logs/hb_movebank_downloadAnnotation_csv_LOG.txt"
    exit
fi

ak_lut=${1}
outdir=${2}
log=${3}

###############################
# #hardcoded for testing:
#ak_lut=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests_accessKeys_lut.csv
#outdir=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/downloads/
###############################

#if logfile doesn't exist, create it:
if [ ! -f ${log} ]; then
	touch ${log}
fi #end logfile if

#cat the file and check for available downloads
i=1
cat ${ak_lut} | while read line
do
	#if the accesskey is listed as "requested," check to see if it's available yet
	ak_status=$(echo $line | cut -d',' -f5)
	if [ ${ak_status} = requested ]; then
		#get accesskey
		ak=$(echo $line | cut -d',' -f1)
		echo "WHRC status = requested: checking Movebank.org for download availability for AccessKey ${ak} now. . ." >> ${log}
		url=http://www.bioinfo.mpg.de/orn-gateway/status.jsp?access-key=${ak}
		#get status page
		stPage=$(curl ${url} | grep "Status" | awk -F ';' '{print $2}')
		
		#if the page is available, download it
		if [ stPage == "Available" ]; then
			track=$(basename $(echo ${line} | cut -d',' -f3) .csv)
			vars=$(basename $(echo ${line} | cut -d',' -f4) .xml)
			today=$(echo ${line} | cut -d',' -f2)
			outfile=${outdir}${track}_${vars}_${today}_${ak}.csv
			dl_url=http://www.bioinfo.mpg.de/orn-gateway/download.jsp?access-key=${ak}
			curl -o ${outfile} ${dl_url}
			
			#after download, update ${ak_lut}
			if [ -f ${outfile} ]; then
				sed -i .bak ''"$i"'s/requested/delivered/' ${ak_lut}
		#otherwise, do nothing - move to next record in $ak_lut
		elif [ stPage == "delivered" ]; then
			now=$(date +%Y%m%d-%H:%M:%S)
			echo "status = delivered: ${ak} has already been downloaded" >> ${log}
		else
			echo "${ak} still pending as of ${now}. . .check back later" >> ${log}
		fi #end stPage if
	fi #end ak_status if
	i=i+1
done #end while