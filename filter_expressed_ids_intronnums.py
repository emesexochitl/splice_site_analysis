#!/usr/bin/python

import datetime
import operator

now = datetime.datetime.now()

print now

ref="../references/Arabidopsis_thaliana.TAIR10.34.total_intronnums.txt"
reflist=[]
with open(ref) as inputfile:
    for line in inputfile:
       reflist.append(line.strip().split('\t'))

ref1=[item[0] for item in reflist]
ref2=[item[1] for item in reflist]
ref3=[item[2] for item in reflist]
ref4=[item[3] for item in reflist]

print reflist[:10]
expr="counts_salt_mock_TPM_vals_expressed.txt"
exprlist=[]
with open(expr) as inputfile:
    for line in inputfile:
       exprlist.append(line.strip().split(' '))

expr1=[item[0] for item in exprlist]

print exprlist[:10]

output="Arabidopsis_thaliana.TAIR10.34_salt_mock_isoform_intron_expressed.txt"
myfile = open(output, "w")

gene = None
isoform = None
intnum = None
ss_id = None
hit = False

print "Merging has been started..."

for i in xrange(0, len(reflist)):
    gene = None
    isoform = None
    intnum = None
    ss_id = None
    hit = False

    for j in xrange(0, len(exprlist)):
        if ref1[i] == expr1[j]:
            hit = True
            gene = ref1[i]
            isoform = ref2[i]
            intnum = ref3[i]
            ss_id = ref4[i]
        else: continue

    if hit == True:        
        #print gene, isoform, intnum, ss_id
        myfile.write("%s\t%s\t%s\t%s\n" % (gene, isoform, intnum, ss_id))
myfile.close()

now = datetime.datetime.now()
print now
