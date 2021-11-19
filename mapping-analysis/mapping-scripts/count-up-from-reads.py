import sys

d_gene = {}

for f in sys.argv[1:]:
    for line in open(f):
        mg_id = f.split('.')[0]
        dat = line.rstrip().split('\t')
        gene = dat[1]
        if d_gene.has_key(gene):
            d_gene[gene][mg_id] = d_gene[gene].get(mg_id,0) + 1
        else:
            d_gene[gene] = {}
            d_gene[gene][mg_id] = 1

fp = open('summary-count.tsv', 'w')

sorted_samples = sys.argv[1:]

for x in sorted_samples:
    fp.write('\t%s' % x.split('.')[0])

fp.write('\n')

for gene in d_gene:
    fp.write('%s\t' % gene)
    for x in sorted_samples:
        x1 = x.split('.')[0]
        if d_gene[gene].has_key(x1):
            fp.write('%s\t' % d_gene[gene][x1])
        else:
            fp.write('0\t')
    fp.write('\n')




