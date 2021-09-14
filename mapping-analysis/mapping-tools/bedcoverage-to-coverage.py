import sys, numpy

d = {}
for line in open(sys.argv[1]):
    dat = line.rstrip().split('\t')
    contig_id = dat[0]
    cov_bp = int(dat[-1])
    if d.has_key(contig_id):
        d[contig_id].append(cov_bp)
    else:
        d[contig_id] = [cov_bp]

sorted_keys = sorted(d.keys())

for key in sorted_keys:
    l = d[key]
    summy = 0
    median = numpy.median(l)
    avg = numpy.average(l)
    minl = min(l)
    maxl = max(l)
    for n, x in enumerate(l):
        if x > 0:
            summy += 1
    total = n + 1
    cov_rat = float(summy)/total
    print key, median, avg, minl, maxl, summy, total, cov_rat
