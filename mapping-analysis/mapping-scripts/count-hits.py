import sys

fp = open(sys.argv[1], 'w')
files = sys.argv[2:]
total_files = len(files)
d = {}

for n, file in enumerate(files):
    for line in open(file):
        line = line.split(' ') 
        contig = line[0]
        max_reads = float(line[1])
        if d.has_key(contig):
            d[contig][n] = max_reads
        else:
            d[contig] = [0]*total_files
            d[contig][n] = max_reads

for x in files:
    fp.write('\t%s' % x)
fp.write('\n')

for key in d.keys():
    fp.write('%s\t' % key)
    for i in range(total_files):
        fp.write('%d\t' % d[key][i])
    fp.write('\n')

        

