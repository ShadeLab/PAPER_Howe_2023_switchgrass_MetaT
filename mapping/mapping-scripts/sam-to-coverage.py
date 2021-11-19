import sys
import screed
import pysam
import numpy

#takes bowtie maps and calculates coverage for contigs.
#output = contig id, base pair coverage median, contig length, read coverage
  
samfile = pysam.Samfile(sys.argv[1], 'rb')
reference_file = sys.argv[2]

def get_coverage(list):
    covered_bp = len([i for i, e in enumerate(list) if e != 0])
    total_bp = len(list)
    return covered_bp, total_bp


d = {}
for record in screed.open(reference_file):
    d[record.name.split(' ')[0]] = len(record.sequence)

for key in d.iterkeys():
    list_of_coverage = [0] * d[key]
    for pileupcolumn in samfile.pileup(key):
        list_of_coverage[pileupcolumn.pos] = pileupcolumn.n 
    covered, total = get_coverage(list_of_coverage)
    print key, numpy.median(list_of_coverage), numpy.average(list_of_coverage), min(list_of_coverage), max(list_of_coverage), covered, total, float(covered)/total 
    
