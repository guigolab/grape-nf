#!/usr/bin/env python

from __future__ import print_function, division
import numpy as np
import argparse
from collections import defaultdict,namedtuple
import gzip
from math import log
import sys
import logging

class Contig:
  """Hold contig information and compute simple scores"""
  def __init__(self,start,stop,score):
    """Chr start/stop and total (pooled) score"""
    self.start = start
    self.stop = stop
    self.score = score

  def __len__(self): return int(self.stop-self.start)
  def bpkm(self,total):
    """Bases per kilobase (length) per million mapped bases (library depth)"""
    return self.score / (self.stop-self.start) * 1000 / total * 1000000

  def bedscore(self, total):
    """Convert bpkm to bedscore filed for UCSC"""
    return min(1000, int(round( 100 * log(100 * self.bpkm(total) + 1))))

  def checkAntisense(self, maxAnti):
    """check for possible shadow antisense artifact"""
    for scr,ascr in zip(self.scores,self.antiscores):
      if ascr * maxAnti < scr: return True
    return False

def load(inp,chrLengths):
  """Load bedgraph file by chromosome"""
  cchr = None
  skip = set()
  for l in inp:
    line = l.split()
    chr = line[0]
    if chr in skip: continue
    if cchr != chr:
      if cchr: yield (cchr,signal)
      if chr in chrLengths:
        signal = np.zeros(chrLengths[chr])
        cchr = chr
      else:
        logging.warn('Chromosome %s not in genome file. Skipping', chr)
        skip.add(chr)
        cchr = None
        continue

    start = int(line[1])
    stop = int(line[2])
    score = float(line[3])
    signal[start:stop] = score
  if cchr: yield (cchr,signal)

def rawcontigs(score, maxGap, minDepth):
  """Call contigs on pooled signal track"""
  pos = np.where(score >= minDepth)[0] #All positions along chr which signal
  dist = np.where(pos[1:]-pos[:-1] > maxGap+1)[0] #All positions followed by a long gap
  if dist.any():
    start = pos[0] #Start at beginning of covered region
    stop = pos[dist[0]]+1 #Until first long gap
    yield (start,stop,score[start:stop].sum()) #Yield first contig

    for i in xrange(1,len(dist)): #Generate contigs
      start = pos[dist[i-1]+1]
      stop = pos[dist[i]] + 1
      yield (start,stop,score[start:stop].sum())
    start = pos[dist[-1]+1] #Last contig
    stop = pos[-1]+1
    yield (start,stop,score[start:stop].sum())

def scoredContigs(scores, minSig, minCov, maxGap, minDepth, antiscores = []):
  """Filter contigs and score against individual tracks"""
  signal = np.sum(scores,0) if scores else 0

  for c in rawcontigs(signal, maxGap, minDepth):
    contig = Contig(*c)
    seq = slice(contig.start,contig.stop)
    if contig.score > minSig and (signal[seq]>0).mean > minCov: #Check contig for score, coverage
      contig.scores = [scr[seq].sum() for scr in scores]
      contig.antiscores = [scr[seq].sum() for scr in antiscores]
      yield contig

def readChrFile(chrFile):
  """Load chrNames and chrLengths from tab-seperated text file"""
  chrs = {}
  for l in chrFile:
    chr, length = l.rstrip().split()[:2]
    chrs[chr] = int(length)
  return chrs

def loadSignal(chrFile, fileP, fileM=None, gz=True):
  """Open bedgraph files"""
  logging.debug(chrFile)
  logging.debug("+:%s\t-:%s",fileP,fileM)
  chrLengths = readChrFile(chrFile)
  inp = {} #File streams for all bedgraph inputs
  if fileM:
    if not len(fileP) == len(fileM): raise ValueError("Unequal number of + and - strand files")
    inp['+'] = [gzip.GzipFile(fileobj=f) if gz else f for f in fileP]
    inp['-'] = [gzip.GzipFile(fileobj=f) if gz else f for f in fileM]
  else:
    inp['u'] = [gzip.GzipFile(fileobj=f) if gz else f for f in fileP]

  files = defaultdict(list) #Generators loading signal for file (one chr at a time)
  signal = defaultdict(dict) #Dictionary containing the signal for all chromosomes and strands
  for strand in inp:
    for file in inp[strand]:
      logging.info("%s\t%s",strand,file.name)
      files[strand].append(load(file,chrLengths)) #Generator reading bedgraph 'file'

  for strand in inp:
    for file in files[strand]:
      logging.debug("Loading next %s ",strand)
      for chr, sig in file:
        if not chr in signal:
          signal[chr] = defaultdict(list)
        signal[chr][strand].append(sig)
          
  for chr in signal:
    if len(signal[chr].keys()) != len(inp.keys()):
      for strand in inp:
        if not strand in signal[chr]:
          signal[chr][strand].append(np.zeros(chrLengths[chr]))
    total = sum([sig.sum() for strand in signal[chr] for sig in signal[chr][strand]])
    yield chr, signal[chr], total

def cmdopts():
  parser = argparse.ArgumentParser()
  parser.add_argument("--gz", action='store_true', help='Input files are gzipped')
  parser.add_argument("--sortOut", action='store_true', help='Output sorted by chromosome, start, stop', default=False)
  parser.add_argument("--chrFile", type=argparse.FileType('r'), help='File with chromosome names\t length', required=True)
  parser.add_argument("--fileP", type=argparse.FileType('r'), nargs='+', help='+ strand or unstranded bedgraphs', required=True)
  parser.add_argument("--fileM", type=argparse.FileType('r'), nargs='*', help='- strand bedgraphs')
  parser.add_argument("--minSig", type=int, default=100, help='minimum pooled score')
  parser.add_argument("--maxAnti", type=float, default=0.1, help='Max Sense / antisense ratio (filter antisense artifacts)')
  parser.add_argument("--minDepth", type=int, default=1, help='Read depth threshold for contigs')
  parser.add_argument("--minCov", type=int, default=0.5, help='Minimal coverage within a contig (<1 due to allowed gaps)')
  parser.add_argument("--maxGap", type=int, default=10, help='Maximum gap length')
  parser.add_argument('-v', "--verbose", dest='verb', action='count',default=0)
  parser.add_argument("--logfile")

  opts = parser.parse_args()
  return parser

def main():
  stranded = bool(opts.fileM)

  id = 1 #Contig id (name)
  contigs = []
  for chr, signal, total in loadSignal(opts.chrFile, opts.fileP, opts.fileM, opts.gz):
    for strand in signal:
      sense = signal[strand]
      antisenseStrand = '-' if strand == '+' else '+'
      anti = signal[antisenseStrand] if stranded else []
      for contig in scoredContigs(sense, opts.minSig, opts.minCov, opts.maxGap, opts.minDepth, anti):
        contig.chr = chr
        contig.strand = strand
        contig.id = id
        id += 1
        if not stranded or contig.checkAntisense(opts.maxAnti):
          contigs.append(contig)
  if opts.sortOut:
    contigs = sorted(contigs, key=lambda c: (c.chr,c.start,c.stop))
  for contig in contigs:
    print(contig.chr,contig.start,contig.stop,contig.id,contig.bedscore(total),
        contig.strand,format(contig.bpkm(total),'.3g'),*(contig.scores+contig.antiscores))

if __name__ == '__main__':
  opts = cmdopts().parse_args()
  logLevel=[logging.WARNING, logging.INFO, logging.DEBUG, 5,1][opts.verb]
  logging.basicConfig(level=logLevel, filename=opts.logfile, format="%(levelname)s:%(filename)s:%(lineno)d %(message)s")
  main()
