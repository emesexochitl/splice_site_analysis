
#!/usr/bin/python

import datetime
import operator

now = datetime.datetime.now()

print now

ref="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_mock_expressed.txt"
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

#print reflist[:10]

intnum="Arabidopsis_thaliana.TAIR10.34_salt_mock_isoform_intron_expressed.txt"
#intnum="test_I.txt"
intnumlist=[]
with open(intnum) as inputfile:
    for line in inputfile:
       intnumlist.append(line.strip().split('\t'))

int1=[item[0] for item in intnumlist]
int2=[item[1] for item in intnumlist]
int3=[item[2] for item in intnumlist]
int4=[int(item[3]) for item in intnumlist]

#print intnumlist[:10]

output="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_mock_expressed_intnums.txt"
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
            intronnum = int(int4[i])

            intcounts[isoform] = intronnum
            break

        elif int4[i] == 0:
            hit = True
            gene_id = int1[i]
            isoform = int2[i]
            ss_id = int3[i]
            chr = "Na"
            start =  "Na" 
            end =  "Na"  
            ori = "Na"  
            annot =  "Na"  
            donor =  "NN"  
            acceptor = "NN"
            intronnum = int(int4[i])

        else: continue

    if hit == True:        
        #print gene_id, isoform, ss_id, chr, start, end, ori, annot, donor, acceptor, intronnum
        myfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_id, isoform, ss_id, chr, start, end, ori, annot, donor, acceptor, intronnum))

myfile.close()

# No. of isoforms!

for k,v in sorted(intcounts.iteritems(),key=operator.itemgetter(0)):
    counts[v] = counts.get(v, 0) + 1

print counts

# Collect SS sites for logos, in separate files.
  
input="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_mock_expressed_intnums.txt"
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

#print inputlist[:10]
#print in11[:10]

for k, v in counts.items():
   #print k

   specoutput="Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_mock_expressed_intnums_%s.txt" % str(k)
   myfile = open(specoutput, "w")

   intnumspec = None
   ss_spec = None


# to make it string-specific, I have to make one more block like this, but with an additional condition regarding in7.

   for l in xrange(0, len(inputlist)):
       intnumspec = None
       ss_spec = None
       if in11[l] == str(k):
           intnumspec =  in3[l] + "_" + in11[l]
           ss_spec =  in9[l] + in10[l]
           myfile.write("%s\t%s\n" % (intnumspec, ss_spec))
           print intnumspec, ss_spec 
           
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
