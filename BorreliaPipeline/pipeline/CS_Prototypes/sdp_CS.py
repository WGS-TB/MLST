import os
import cplex
from cplex.exceptions import CplexError
import sys
import csv

# change the path to the adp's output for a sample
# os.chdir(sys.argv[1])
os.chdir("./sample_data")
print
#----------------------------------------------------------------------
#               Creating the matrices used in LP
#----------------------------------------------------------------------

print "----------------------------------------------------------------------\n               Creating the matrices " \
      "used in LP\n---------------------------------------------------------------------- "
print

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

#---------------------- Reading the output of adp ---------------------


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

#--------------------- Generating the Coeff. matrix -------------------

remainder = cols
for d in dims:
    for r in range(d):
        row = []
        for repeat in range(cols / (remainder / d)):
            for i in range(remainder / d):
                if repeat % d == r:
                    row.append(1)
                else:
                    row.append(0)
        A.append(row)
    remainder /= d
print "    Y:\n", Y, "\n    A"
for row in A:
    print row
print
print

#----------------------------------------------------------------------
#                     Solving the CS using CPLEX
#----------------------------------------------------------------------

print "----------------------------------------------------------------------\n                     Solving the CS " \
      "using CPLEX\n---------------------------------------------------------------------- "
print

#----------------------- Preparing the Inputs -------------------------

epsilon = 0.0001 # this can be changed later to determine accuracy
obj = [1.0 for x in range(cols)]
ub = [1.0 for x in range(cols)]
cnames = ['x' + str(x) for x in range(cols)]
rnames = ['c' + str(x) for x in range(len(A) * 2)]

rhs = []
rhs.extend(Y)
for y in Y:
    rhs.append(-y)
for i in range(len(rhs)):
    rhs[i] += epsilon
constraints = []
for row in A:
    constraints.append([range(cols), row])
for row in A:
    constraints.append([range(cols), [-x for x in row]])

senses = ""
for i in range(len(constraints)):
    senses += "L"

print "    Constraints:"
for r in range(len(constraints)):
    print constraints[r][1], "\t\t\t", rhs[r]
print "\n"

#-------------------- Using cplex to solve the LP ----------------------

try:
    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.variables.add(obj=obj, ub=ub, names=cnames)
    prob.linear_constraints.add(lin_expr=constraints, senses=senses, rhs=rhs, names=rnames)

    prob.solve()

    print
    print "Solution status:", prob.solution.get_status(), prob.solution.status[prob.solution.get_status()]
    print "Solution value  = ", prob.solution.get_objective_value()
    for i in range(prob.variables.get_num()):
        print prob.variables.get_names()[i], '\t=\t', prob.solution.get_values()[i]


except CplexError, exc:
    print exc
    exit(-1)
