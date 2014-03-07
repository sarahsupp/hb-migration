#!/usr/bin/python

# Purpose: This script creates xml files that will be submitted to movebank as annotation requests.
# Author: Tina A. Cormier
# Date: 03/06/2014
# Status: In progress.

import os
import xml.etree.ElementTree as ET
from xml.dom import minidom
import os

# variable names you want to extract (these are the names that will appear in the 
# output table - can name them as you wish).
var_labels = "NDVI,EVI"
# Interpolation type - one of "inverse-distance-weighted", "bilinear", or "nearest-neighbor"
int_type = "inverse-distance-weighted"
# Type Name - need directory and sensor type (e.g. "modis-land/MOD13Q1.005")
type_name = "modis-land/MOD13Q1.005"
# Variable names - must be names listed in movebank
var_names = "250m 16 days NDVI,250m 16 days EVI"
# Set output directory for xml file
outdir = "/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/xml/"

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
def CreateXML(var_labels, int_type, type_name, var_names, outfile):
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
    arep5.set("value", var_names)
    messages = ET.SubElement(root,"Messages")
    
    #print prettify(root)
    
    output_file = open( outfile, 'w' )
    output_file.write( '<?xml version="1.0"?>' )
    output_file.write( ET.tostring( root ) )
    output_file.close()
