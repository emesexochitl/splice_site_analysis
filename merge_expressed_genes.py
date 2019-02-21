#!/usr/bin/python
# Get the details of stress-induced upregulated genes and their splice sites.

import datetime
import operator

now = datetime.datetime.now()

print now

# Use the extended SS isoform refereence


sample = "50mM"

ref="../references/8_ss/Arabidopsis_thaliana.TAIR10.34_isoform_SS_background_8.txt"
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
ref7=[item[6] for item in reflist]
ref8=[item[7] for item in reflist]
ref9=[item[8] for item in reflist]
ref10=[item[9] for item in reflist]

#print reflist[:10]

# Here expression values are coming from the R TPM output, filtered to the RPKM ==1 TPM value ( usually around ~ 1.6).

expr="counts_salt_%s_TPM_vals_expressed.txt" % sample
exprlist=[]
with open(expr) as inputfile:
    for line in inputfile:
       exprlist.append(line.strip().split(' '))
#print exprlist[:10]

expr1=[item[0] for item in exprlist]
expr2=[item[3] for item in exprlist]
expr3=[item[4] for item in exprlist]

print expr1[:10]
print expr2[:10]
print expr3[:10]

output="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_%s_expressed.txt" % sample
myfile = open(output, "w")

gene_id = None
isoform = None
ss_id = None
chr = None
start = None
end = None
ori = None
annot = None
donor = None
acceptor = None
tpmval = None
log2tpm = None
hit = False

print "Merging has been started..."

for i in xrange(0, len(reflist)):
    gene_id = None
    isoform = None
    ss_id = None
    chr = None
    start = None
    end = None
    ori = None
    annot = None
    donor = None
    acceptor = None
    tpmval = None
    log2tpm = None
    hit = False

    for j in xrange(0, len(exprlist)):
        if ref1[i] == expr1[j]:
            hit = True
            gene_id = ref1[i]
            isoform = ref2[i]
            ss_id = ref3[i]
            chr = ref4[i]
            start = ref5[i]
            end = ref6[i]
            ori = ref7[i]
            annot = ref8[i]
            donor = ref9[i]
            acceptor = ref10[i]
            tpmval = expr2[j]
            log2tpm = expr3[j]

        else: continue

    if hit == True:        
        #print gene_id, isoform, ss_id, chr, start, end, ori, annot, donor, acceptor
        myfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_id, isoform, ss_id, chr, start, end, ori, annot, donor, acceptor, tpmval, log2tpm))
myfile.close()

now = datetime.datetime.now()
print now
