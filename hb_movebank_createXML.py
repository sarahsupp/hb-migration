#!/usr/bin/python

# Purpose: This script creates xml files that will be submitted to movebank as annotation requests.
# Author: Tina A. Cormier
# Date: 03/06/2014
# Status: In progress.

import os
import xml.etree.ElementTree as ET
from xml.dom import minidom
import argparse

#Deal with passed arguments
p = argparse.ArgumentParser(prog="hb_movebank_createXML.py", description="generate annotation request xml for submission to movebank.")
p.add_argument("-vl", dest="var_labels", required=True, help="String list of variable names. These are the column names that will appear in your output table. e.g., 'Quality,NDVI,EVI'")
p.add_argument("-i", dest="int_type", required=True, help="Interpolation type. Must be nearest-neighbor, bilinear, or inverse-distance-weighted")
p.add_argument("-t", dest="type_name", required=True, help="movebank directory and sensor type (e.g., 'modis-land/MOD13Q1.005'")
p.add_argument("-vn", dest="mb_var_names", required=True, help="variable names exactly as they are listed in movebank (e.g.,'250m 16 days VI Quality,250m 16 days NDVI,250m 16 days EVI'")
p.add_argument("-o", dest="outdir", required=True, help="output directory for xml file - use trailing slash")

args = p.parse_args()

var_labels = args.var_labels
int_type = args.int_type
type_name = args.type_name
mb_var_names = args.mb_var_names
outdir = args.outdir

print var_labels
print int_type
print type_name
print mb_var_names
print outdir

############### Hardcoded Variables for Testing ###############
## variable names you want to extract (these are the names that will appear in the
## output table - can name them as you wish).
##var_labels = "lwrf,swrf,t10m,u10m,v10m"
#var_labels = "Quality,NDVI,EVI"

## Interpolation type - one of "inverse-distance-weighted", "bilinear", or "nearest-neighbor"
#int_type = "inverse-distance-weighted"

## Type Name - need directory and sensor type (e.g. "modis-land/MOD13Q1.005")
#type_name = "modis-land/MOD13Q1.005"
##type_name = "nomads.ncdc.noaa.gov/NARR"

## Variable names - must be names listed in movebank
#mb_var_names = "250m 16 days VI Quality,250m 16 days NDVI,250m 16 days EVI"
##mb_var_names = "Downward_longwave_radiation_flux_sfc,Downward_shortwave_radiation_flux_sfc,Temp._10_m_above_gnd,u_wind_10_m_above_gnd,v_wind_10_m_above_gnd"

## Set output directory for xml file
#outdir = "C:/Share/tcormier/hummingbirds/migration_study/movebank/env_request_xmls/"

##############################

#set output file
#first add trailing slash if missing
outdir = os.path.normpath(outdir) + os.sep
outfile = outdir + '_'.join(var_labels.split(",")) + ".xml"


########
#function for xml pretty printing (from: http://pymotw.com/2/xml/etree/ElementTree/create.html).
def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

########
#function to create xml from user-input variables.
def CreateXML(var_labels, int_type, type_name, mb_var_names, outfile):
    #Create xml
    root = ET.Element("AnnotationRequest")
    root.set("annotationType","track2")
    root.set("name","test_annotation")
    root.set("notificationEmail","tcormier@whrc.org")
    root.set("owner","tcormier")
    elements = ET.SubElement(root, "Elements")
    are = ET.SubElement(elements, "AnnotationRequestElement")
    are.set("variableLabel", var_labels)
    properties = ET.SubElement(are, "Properties")
    arep = ET.SubElement(properties, "AnnotationRequestElementProperty")
    arep.set("name", "interpolation-type")
    arep.set("value", int_type)
    arep2 = ET.SubElement(properties, "AnnotationRequestElementProperty")
    arep2.set("name","type-name")
    arep2.set("value", type_name)
    arep3 = ET.SubElement(properties, "AnnotationRequestElementProperty")
    arep3.set("name", "type")
    arep3.set("value","simple")
    arep4 = ET.SubElement(properties, "AnnotationRequestElementProperty")
    arep4.set("name","distance-function-spatial")
    arep4.set("value","geodetic")
    arep5 = ET.SubElement(properties, "AnnotationRequestElementProperty")
    arep5.set("name","variable-names")
    arep5.set("value", mb_var_names)
    messages = ET.SubElement(root,"Messages")

    #print prettify(root)

    output_file = open( outfile, 'w' )
    output_file.write( '<?xml version="1.0"?>' )
    output_file.write( ET.tostring( root ) )
    output_file.close()

#create xml
CreateXML(var_labels, int_type, type_name, mb_var_names, outfile)