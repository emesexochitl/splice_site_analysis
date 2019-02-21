#!/usr/bin/python

import datetime
import operator
import math
import scipy
from scipy import stats
from numpy import *
import os
from Bio import SeqIO

now = datetime.datetime.now()

#print now

class iso_set():

    def __init__(self,iso_id,exon_num):

        # Collect known informations for GTF based file.

        self.iso_id = iso_id[:9]
        self.chr = [x[1] for x in intron_list if x[0] == iso_id][0]
        self.intron_lefts = list(set([int(x[2]) for x in intron_list if x[0] == iso_id]))
        self.intron_rights = list(set([int(x[3]) for x in intron_list if x[0] == iso_id]))
        #self.intron_lefts = [int(x[2]) for x in intron_list if x[0] == iso_id]
        #self.intron_rights = [int(x[3]) for x in intron_list if x[0] == iso_id]

        self.ori = [x[4] for x in intron_list if x[0] == iso_id][0]
         
        self.linear_introns()
        self.extract_splice_sites()




    def linear_introns(self):

        self.chr_hit = chr_list.index(self.chr)

        for i in xrange(0, len(intron_list)):

            if self.ori == "+":

                self.upstream_starts = [x-4 for x in self.intron_lefts] ### ok
                self.upstream_ends = [x+20 for x in self.intron_lefts] ### ok
                self.downstream_starts = [x-21 for x in self.intron_rights] ### ok
                self.downstream_ends = [x+3 for x in self.intron_rights] ### ok

            if self.ori == "-":

                self.upstream_ends = [x+3 for x in self.intron_rights] ### ok
                self.upstream_starts = [x-21 for x in self.intron_rights] ### ok
                self.downstream_ends = [x+20 for x in self.intron_lefts] ### ok
                self.downstream_starts = [x-4 for x in self.intron_lefts] ### ok


    def extract_splice_sites(self):

        self.chr_hit = chr_list.index(self.chr)

        for i in xrange(0,len(self.upstream_starts)):

            if self.ori == "+" and len(self.upstream_starts) > 0:
                with  open(out_upstream, 'a') as f:
                    f.write(">%s_%d\t%s\n" % (self.iso_id,i+1,(records[self.chr_hit].seq[int(self.upstream_starts[i]):int(self.upstream_ends[i])])))

            if self.ori == "-" and len(self.upstream_starts) > 0:
                with  open(out_upstream, 'a') as f:
                    f.write(">%s_%d\t%s\n" % (self.iso_id,i+1,(records[self.chr_hit].seq[int(self.upstream_starts[i]):int(self.upstream_ends[i])]).reverse_complement()))

            else:break

        for i in xrange(0,len(self.downstream_starts)):

            if self.ori == "+" and len(self.downstream_starts) > 0:
                with  open(out_downstream, 'a') as f:
                    f.write(">%s_%d\t%s\n" % (self.iso_id,i+1,(records[self.chr_hit].seq[int(self.downstream_starts[i]):int(self.downstream_ends[i])])))

            if self.ori == "-" and len(self.downstream_starts) > 0:
                with  open(out_downstream, 'a') as f:
                    f.write(">%s_%d\t%s\n" % (self.iso_id,i+1,(records[self.chr_hit].seq[int(self.downstream_starts[i]):int(self.downstream_ends[i])]).reverse_complement()))

            else: break

try:
    os.remove(out_upstream)
except:
    pass

try:
    os.remove(out_downstream)
except:
    pass



out_upstream = "linear_exp_upstream_nonred2.clustalw"
out_downstream = "linear_exp_downstream_nonred2.clustalw"

### gtf or something

ref_file = "merge_RnaseR_polyA_loose_forplots.txt"
ref_list = []
with open(ref_file) as inputfile:
    for line in inputfile:
        ref_list.append(line.strip().split('\t'))


### need expressed ids!!!

intron_file = "Arabidopsis_thaliana.TAIR10.34_isoform_intron_unique.saf"
intron_list = []
with open(intron_file) as inputfile:
    for line in inputfile:
        intron_list.append(line.strip().split('\t'))

intron_id =[item[0] for item in intron_list]


exp_file = "average_exp1_exon_len_num2.txt"
#exp_file="test.txt"
exp_list = []
with open(exp_file) as inputfile:
    for line in inputfile:
        exp_list.append(line.strip().split('\t'))

exp_id =[item[0][:9] for item in exp_list if item[0][:9] in intron_id]


### circ instead of ref.

circ_file = "merge_RnaseR_polyA_loose_forplots.txt"
circ_list = []
with open(circ_file) as inputfile:
    for line in inputfile:
        circ_list.append(line.strip().split('\t'))

records = list(SeqIO.parse("../../references/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa","fasta"))

chr_list = ["1","2","3","4","5","Mt","Pt"]

#intron_id =[item[0] for item in intron_list]

circ_id =[item[0][:9] for item in circ_list]


#ids = [x for x in intron_id if x not in circ_id]
ids = list(set([x for x in exp_id if x not in circ_id]))


print exp_id[:5], circ_id[:5], ids[:5]
print len(intron_id), len(circ_id), len(exp_id), len(ids)

counts = dict()

transformtuple = tuple(ids)

for p in transformtuple:
   counts[p] = counts.get(p, 0) + 1


for k,v in counts.items():
    print k, v

    iso_id = iso_set(k,v)
    print iso_id.__dict__
