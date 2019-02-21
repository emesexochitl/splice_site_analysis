#!/usr/bin/python
# Prepare separate llist of upregulated introm numbers, because genes with 0 intron don't have splice sites and therefore dropping out of the SS analysis and also I have to merge upregulated genes' SS and introm number data.

import datetime
import operator

now = datetime.datetime.now()

print now

# It's been generated.

ref="../references/Arabidopsis_thaliana.TAIR10.34.total_intronnums_mod.txt"
reflist=[]
with open(ref) as inputfile:
    for line in inputfile:
       reflist.append(line.strip().split('\t'))

ref1=[item[0] for item in reflist]
ref2=[item[1] for item in reflist]
ref3=[item[2] for item in reflist]
ref4=[item[3] for item in reflist]
ref5=[item[4] for item in reflist]
ref6=[item[5] for item in reflist]

print reflist[:10]

# Here expression values are coming from the DESeq output, with padj < 0.1 and log2foldchange > 1.

expr="deseq_Mock_150_mM_fc1up_01.txt"
exprlist=[]
with open(expr) as inputfile:
    for line in inputfile:
       exprlist.append(line.strip().split(' '))

expr1=[item[1] for item in exprlist]
expr2=[item[6] for item in exprlist]
expr3=[item[7] for item in exprlist]

print expr1[:3]
print expr2[:3]
print expr3[:3]

output="Arabidopsis_thaliana.TAIR10.34_salt_150mM_isoform_intron_upreg.txt"
myfile = open(output, "w")

gene = None
isoform = None
intnum = None
ss_id = None
ori = None
annot = None
pval = None
padj = None
hit = False

print "Merging has been started..."

for i in xrange(0, len(reflist)):
    gene = None
    isoform = None
    intnum = None
    ss_id = None
    ori = None
    annot = None
    pval = None
    padj = None
    hit = False

    for j in xrange(0, len(exprlist)):
        if ref1[i] == expr1[j]:
            hit = True
            gene = ref1[i]
            isoform = ref2[i]
            intnum = ref3[i]
            ss_id = ref4[i]
            ori = ref5[i]
            annot = ref6[i]
            pval = expr2[j]
            padj = expr3[j]
        else: continue

    if hit == True:        
        #print gene, isoform, intnum, ss_id
        myfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, isoform, intnum, ss_id, ori, annot, pval, padj))
myfile.close()

now = datetime.datetime.now()
print now
