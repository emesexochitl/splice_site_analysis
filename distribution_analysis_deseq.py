
#!/usr/bin/python

import datetime
import operator

now = datetime.datetime.now()

print now

# Reference is going to be the detailed table of upregulated genes (merge_genes_deseq.py). Now I add the number of introns to it, and separate the output according to the represented intron numbers and prepare outputs for WebLogo visualization. Later it could be modified according to other parameters, such as strand orientation.
# It will continue the intronless hits as well, but with missing values, bebause the intron number list is the base of the loop.

sample = "300mM"

ref="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_%s_upreg.txt" % sample
#ref="test_expressed.txt"
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
ref11=[item[10] for item in reflist]
ref12=[item[11] for item in reflist]

#print reflist[:10]

# Intron number is coming from filter_deseq_ids_intronnums.py, because of intronless genes.

intnum="Arabidopsis_thaliana.TAIR10.34_salt_%s_isoform_intron_upreg.txt" % sample
#intnum="test_I.txt"
intnumlist=[]
with open(intnum) as inputfile:
    for line in inputfile:
       intnumlist.append(line.strip().split('\t'))

int1=[item[0] for item in intnumlist]
int2=[item[1] for item in intnumlist]
int3=[item[2] for item in intnumlist]
int4=[int(item[3]) for item in intnumlist]
int5=[item[4] for item in intnumlist]
int6=[item[5] for item in intnumlist]
int7=[item[6] for item in intnumlist]
int8=[item[7] for item in intnumlist]

#print intnumlist[:10]

output="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_%s_upreg_intnums.txt" % sample
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
intronnum = None
pval = None
padj = None
hit = False
counts = dict()
intcounts = dict()

print "Merging has been started..."

for i in xrange(0, len(intnumlist)):
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
    pval = None
    padj = None
    hit = False

    for j in xrange(0, len(reflist)):

        if int3[i] == ref3[j]:
            hit = True
            gene_id = ref1[j]
            isoform = ref2[j]
            ss_id = ref3[j]
            chr = ref4[j]
            start = ref5[j]
            end = ref6[j]
            ori = ref7[j]
            annot = ref8[j]
            donor = ref9[j]
            acceptor = ref10[j]
            pval = ref11[j]
            padj = ref12[j]
            intronnum = int(int4[i])

            intcounts[isoform] = intronnum
            break

        elif int4[i] == 0:
            hit = True
            gene_id = int1[i]
            isoform = int2[i]
            ss_id = int3[i]
            chr = int4[i]
            start =  "Na" 
            end =  "Na"  
            ori = int5[i]  
            annot =  int6[i]  
            donor =  "NN"  
            acceptor = "NN"
            pval = int7[i]
            padj = int8[i]

            intronnum = int(int4[i])
            intcounts[isoform] = intronnum

        else: continue

    if hit == True:        
        #print gene_id, isoform, ss_id, chr, start, end, ori, annot, donor, acceptor, intronnum
        myfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_id, isoform, ss_id, chr, start, end, ori, annot, donor, acceptor, pval, padj, intronnum))

myfile.close()

# No. of isoforms!

for k,v in sorted(intcounts.iteritems(),key=operator.itemgetter(0)):
    counts[v] = counts.get(v, 0) + 1

print counts

# Collect SS sites for logos, in separate files.
  
input="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_%s_upreg_intnums.txt" % sample
inputlist=[]
with open(input) as inputfile:
    for line in inputfile:
       inputlist.append(line.strip().split('\t'))

in1=[item[0] for item in inputlist]
in2=[item[1] for item in inputlist]
in3=[item[2] for item in inputlist]
in4=[item[3] for item in inputlist]
in5=[item[4] for item in inputlist]
in6=[item[5] for item in inputlist]
in7=[item[6] for item in inputlist]
in8=[item[7] for item in inputlist]
in9=[item[8] for item in inputlist]
in10=[item[9] for item in inputlist]
in11=[item[10] for item in inputlist]
in12=[item[11] for item in inputlist]
in13=[item[12] for item in inputlist]

#print inputlist[:10]
#print in11[:10]

for k, v in counts.items():
   #print k

   specoutput="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_%s_upreg_intnums_%s.txt" % (sample, str(k))
   myfile = open(specoutput, "w")

   intnumspec = None
   ss_spec = None


# to make it string-specific, I have to make one more block like this, but with an additional condition regarding in7.

   for l in xrange(0, len(inputlist)):
       intnumspec = None
       ss_spec = None
       if in13[l] == str(k):
           intnumspec =  in3[l] + "_" + in13[l]
           ss_spec =  in9[l] + in10[l]
           myfile.write("%s\t%s\n" % (intnumspec, ss_spec))
           #print intnumspec, ss_spec 
           
       else: continue   

   myfile.close()

#print sorted(counts.iteritems(),key=operator.itemgetter(int(0)))
#print sorted(counts.iteritems(),key=operator.itemgetter(1))

#print intcounts
now = datetime.datetime.now()
print now

#for k,v in sorted(intcounts.iteritems(),key=operator.itemgetter(0)):
#    counts[v] = counts.get(v, 0) + 1

#print counts   
