'''
Created on Feb 24, 2012

@author: mkiyer
'''
import logging
import argparse
import os
import subprocess
import xml.etree.ElementTree as etree 

from base import check_executable, check_sam_file

RETCODE_SUCCESS = 0
RETCODE_ERROR = 1
DESCRIPTION = "Chimerascan2 chimeric transcript (gene fusion) discovery tool"
DEFAULT_OUTPUT_DIR = "chimerascan2_out"
CONFIG_XML_ROOT_NAME = "chimerascan2_config"

# library strand types
LIB_FR_UNSTRANDED = 0
LIB_FR_FIRSTSTRAND = 1
LIB_FR_SECONDSTRAND = 2

#def get_genome_orientation(r, library_type):
#    if library_type == LibraryTypes.FR_FIRSTSTRAND:
#        if r.is_read2:
#            return OrientationTags.FIVEPRIME
#        else:
#            return OrientationTags.THREEPRIME
#    elif library_type == LibraryTypes.FR_SECONDSTRAND:
#        if r.is_read1:
#            return OrientationTags.FIVEPRIME
#        else:
#            return OrientationTags.THREEPRIME
#    return OrientationTags.NONE



def indent_xml(elem, level=0):
    """ in-place XML prettyprint formatter """
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_xml(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def get_command_line_parser():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--config", dest="config_file", default=None,
                        help="Config file (XML format) with chimerascan "
                        "parameters")
    parser.add_argument("-o", "--output_dir", dest="output_dir", 
                        default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("bam_file")


class RunConfig(object):
    
    attrs = (("velveth_bin", "velveth"),
             ("velvetg_bin", "velvetg"),
             ("samtools_bin", "samtools"),
             ("config_file", None),
             ("bam_file", None),
             ("output_dir", None))
    
    def __init__(self):
        for attrname, attrval in RunConfig.attrs:
            setattr(self, attrname, attrval)

    @staticmethod
    def from_command_line():
        parser = get_command_line_parser()
        config = RunConfig()
        return config.update_command_line(parser)

    def update_command_line(self, parser):
        args = parser.parse_args()
        if (args.config_file is not None):
            if not os.path.exists(args.config_file):
                parser.error("config file %s not found" % (args.config_file))
            self.update_xml(args.config_file)
    
    def update_xml(self, xmlfile):
        tree = etree.parse(xmlfile)        
        root = tree.getroot()
        for attrname, defaultval in RunConfig.attrs:
            attrval = root.findtext(attrname)
            if (attrval is None) or (not attrval) or (attrval == "None"):
                attrval = defaultval
            setattr(self, attrname, attrval)

    def to_xml(self):
        root = etree.Element(CONFIG_XML_ROOT_NAME)
        for attrname, defaultval in RunConfig.attrs:
            val = getattr(self, attrname)
            if val is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(val)
        # indent for pretty printing
        indent_xml(root)        
        return etree.tostring(root)

    def is_valid(self):
        valid = True
        # check third-party software programs
        for attrname in ("velveth_bin", "velvetg_bin", "samtools_bin"):
            prog = getattr(self, attrname)
            if check_executable(prog):
                logging.debug("Checking for '%s' binary... found" % prog)
            else:
                logging.error("'%s' binary not found or not executable" % prog)
                valid = False
        #
        # check third-party python packages
        #
        # check for bx python library
        try:
            import bx.intervals.intersection
            logging.debug("Checking for 'bx python' library... found")
        except ImportError, e:
            logging.error("Package 'bx python' not found")
            valid = False
        # check for pysam library
        try:
            import pysam
            logging.debug("Checking for 'pysam' library... found")
        except ImportError, e:
            logging.error("Package 'pysam' not found")
            valid = False
        # check sorted abundant reads bam file
        if not check_sam_file(self.bam_file, isbam=True):
            logging.error("Input file %s missing/corrupt" % (self.bam_file))
            valid = False
        return valid
