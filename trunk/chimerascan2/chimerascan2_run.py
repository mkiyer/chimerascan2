'''
Created on Feb 24, 2012

@author: mkiyer
'''
import sys
import os
import logging

import pysam
from bx.intervals.cluster import ClusterTree

from chimerascan2.lib import config
from chimerascan2.lib.sam import get_aligned_intervals

def cluster_chrom(riter):
    locus_start = 0
    locus_end = 0
    strand_cluster_trees = {"+": ClusterTree(0,1),
                            "-": ClusterTree(0,1)}
    for i,r in enumerate(riter):
        if r.is_unmapped or r.is_qcfail:
            continue
        # this read marks the beginning of a new locus
        if r.pos > locus_end:
            
            strand_cluster_trees = {"+": ClusterTree(0,1),
                                    "-": ClusterTree(0,1)}
            

        intervals = get_aligned_intervals(r)
        strand = r.opt("XS")
        
        if intervals[0][0] > locus_end:
        


def cluster_bam_file(bamfh):

    for rname in bamfh.references:
        cluster_chrom(bamfh.fetch(rname))


def run(config):
    bamfh = pysam.Samfile(config.bam_file, "rb")
    bamfh.close()
    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    config = config.RunConfig.from_command_line()
    run(config)
    return config.RETCODE_SUCCESS

if __name__ == '__main__':
    sys.exit(main())