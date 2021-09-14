import screed, sys

for record in screed.open(sys.argv[1]):
    l = [record.name.rstrip(), str(int(1)), str(len(record.sequence))]
    print '\t'.join(l)
      
