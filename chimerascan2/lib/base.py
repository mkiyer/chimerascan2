'''
Created on Feb 24, 2012

@author: mkiyer
'''
import os
import subprocess

def check_executable(filename):
    # check that samtools binary exists
    devnullfh = open(os.devnull, 'w')        
    try:
        p = subprocess.Popen([filename], stdout=devnullfh, stderr=devnullfh)
        p.kill()
    except OSError:
        return False
    devnullfh.close()
    return True

def check_sam_file(filename, isbam=False):
    is_valid = True
    if not file_exists_and_nz_size(filename):
        is_valid = False
    else:
        import pysam
        try:
            fmt = "rb" if isbam else "r"
            samfh = pysam.Samfile(filename, fmt)
            samfh.close()   
        except:
            is_valid = False
    return is_valid

def up_to_date(outfile, infile, nzsize=True):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if nzsize and (os.path.getsize(outfile) == 0):
        return False
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

def file_exists_and_nz_size(filename):
    if filename is None:
        return False
    if not os.path.exists(filename):
        return False
    if os.path.getsize(filename) == 0:
        return False
    return True

def parse_bool(s):    
    return True if s[0].lower() == "t" else False

def parse_string_none(s):
    return None if s == "None" else s
