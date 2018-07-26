import os
import cplex
import sys
import csv

# change the path below to the adp's output for a sample
os.chdir(sys.argv[1])
# os.chdir("C:\\Users\\parha\\Downloads\\Compressed\\SRR2034333")

#----------------------------------------------------------------------
#               Creating the matrices used in LP
#----------------------------------------------------------------------

propFiles = []
for (dirpath, dirnames, filenames) in os.walk(os.getcwd()):
    for fn in filenames:
        if '.csv' in fn:
            propFiles.append(fn)
    break
# print propFiles, '\n'

cols = 1
dims = []
Y = []
A = []

#---------------------- Reading the output of adp ----------------------

for pf in propFiles:
    with open(pf, 'rb') as csvFile:
        alleleReader = csv.reader(csvFile)
        locus_alleles = 0
        for allele_prop in alleleReader:
            locus_alleles += 1
            # print allele_prop[1]
            Y.append(float(allele_prop[1]))
        dims.append(locus_alleles)
        cols *= locus_alleles
# print cols, dims

#--------------------- Generating the Coeff. matrix --------------------

remainder = cols
for d in dims:
    for r in range(d):
        row = []
        for repeat in range(cols / (remainder / d)):
            for i in range(remainder / d):
                if repeat % d == r:
                    row.append(1)
                else: row.append(0)
        A.append(row)
    remainder /= d
# print Y
# for row in A:
#     print row

#----------------------------------------------------------------------
#                     Solving the CS using CPLEX
#----------------------------------------------------------------------

