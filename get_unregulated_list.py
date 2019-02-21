
#!/usr/bin/python

import datetime
import operator

now = datetime.datetime.now()

print now

expr="counts_salt_mock_TPM_vals_expressed.txt"
exprlist=[]
with open(expr) as inputfile:
    for line in inputfile:
       exprlist.append(line.strip().split(' '))

expr1=[item[0] for item in exprlist]
expr2=[float(item[3]) for item in exprlist]
print expr1[:3]
print expr2[:3]

upreg="deseq_Mock_50mM_upreglist.txt"
upreglist=[]
with open(upreg) as inputfile:
    for line in inputfile:
       upreglist.append(line.strip().split('\t'))

#upreg=[item[0] for item in upreglist]
print upreglist[:3]
output="Mock_50mM_expr_noninduced.txt"
myfile = open(output, "w")

gene_id = None
tpm_val =  None

for i in xrange(0, len(exprlist)):
    gene_id = None
    tpm_val = None

    for j in xrange(0, len(upreglist)):

        if expr1[i] == upreglist[j]:
            break
        elif expr1[i] != upreglist[j] and expr2[i] > 1.65:
            gene_id = expr1[i]
            tpm_val = expr2[i]
        else: continue

    print gene_id, tpm_val
    myfile.write("%s\t%.3f\n" % (gene_id, tpm_val))

myfile.close()


