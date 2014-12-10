#!/usr/bin/env python
'''=================================================================================================
Infer RNA-seq experiment design from SAM/BAM file. This module will determine if the RNA-seq
experiment is:
1) pair-end or single-end
2) if experiment is strand-specific, how reads were stranded.
 * For pair-end RNA-seq, there are two different ways to strand reads:
  i) 1++,1--,2+-,2-+
     read1 mapped to '+' strand indicates parental gene on '+' strand
     read1 mapped to '-' strand indicates parental gene on '-' strand
     read2 mapped to '+' strand indicates parental gene on '-' strand
     read2 mapped to '-' strand indicates parental gene on '+' strand
  ii) 1+-,1-+,2++,2--
     read1 mapped to '+' strand indicates parental gene on '-' strand
     read1 mapped to '-' strand indicates parental gene on '+' strand
     read2 mapped to '+' strand indicates parental gene on '+' strand
     read2 mapped to '-' strand indicates parental gene on '-' strand
 * For single-end RNA-seq, there are two different ways to strand reads:
  i) ++,--
     read mapped to '+' strand indicates parental gene on '+' strand
     read mapped to '-' strand indicates parental gene on '-' strand
  ii) +-,-+
     read mapped to '+' strand indicates parental gene on '-' strand
     read mapped to '-' strand indicates parental gene on '+' strand
================================================================================================='''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
    print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
    sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
from time import strftime

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012. All rights reserved."
__credits__ = []
__license__ = "GPL"
__version__="2.3.5"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
    '''print progress into stderr and log file'''
    mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    LOG=open('class.log','a')
    print >>sys.stderr,mesg
    print >>LOG,mesg


def main():
    usage="%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Input alignment file in SAM or BAM format")
    parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed fomat.")
    parser.add_option("-s","--sample-size",action="store",type="int",dest="sample_size",default=200000, help="Number of reads sampled from SAM/BAM file. default=%default")
    (options,args)=parser.parse_args()

    if not (options.input_file and options.refgene_bed):
        parser.print_help()
        print >>sys.stderr, '\n\n' + __doc__
        sys.exit(0)
    for f in (options.input_file,options.refgene_bed):
        if not os.path.exists(f):
            print >>sys.stderr, '\n\n' + f + " does NOT exists" + '\n'
            sys.exit(0)
    if options.sample_size <1000:
        print >>sys.stderr, "Warn: Sample Size too small to give a accurate estimation"
    obj = SAM.ParseBAM(options.input_file)
    (protocol,sp1,sp2,other)=obj.configure_experiment(refbed=options.refgene_bed, sample_size = options.sample_size)
    if other <0: other=0.0
    if protocol == "PairEnd":
        if sp1 > 0.8:
            print "MATE1_SENSE"
        elif sp2 > 0.8:
            print "MATE2_SENSE"
        else:
            print "NONE"
    elif protocol == "SingleEnd":
        if sp1 > 0.8:
            print "SENSE"
        elif sp2 > 0.8:
            print "ANTISENSE"
        else:
            print "NONE"
    else:
        print "Unknown Data type"
    #print mesg

if __name__ == '__main__':
    main()
