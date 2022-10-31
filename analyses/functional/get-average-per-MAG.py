import sys, numpy

d = {}

for n, line in enumerate(open(sys.argv[1])):
    if n == 0:
        list_of_metags = line.rstrip().split(',')
    else:
        dat = line.rstrip().split(',')
        id = dat[0]
        mag_id = id.split('_')
        mag_id2 = mag_id[2] + mag_id[3].split('.')[1]
        if mag_id2 in d:
            for i,x in enumerate(dat[1:]):
                d[mag_id2][i].append(float(x))
        else:
            d[mag_id2] = [0]*len(list_of_metags)
            for i, x in enumerate(dat[1:]):
                d[mag_id2][i] =[float(x)]


new_data = {}
for key in d:
    new_data[key] = []
    for x in d[key]:
        new_data[key].append(numpy.mean(x))
print('\t' + '\t'.join(list_of_metags))
for x in new_data:
    print(x + '\t' + '\t'.join(str(item) for item in new_data[x]))
