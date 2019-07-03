from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cplex
from cplex.exceptions import CplexError
from pipeline_functions import *
from scipy.spatial.distance import hamming
from spgl1 import spg_bpdn

matplotlib.use('Agg')
plt.style.use('ggplot')

# ----------------------------------------------------------------------
#               Creating the matrices used in LP
# ----------------------------------------------------------------------


def defineProblem(cols, dims):
    print "----------------------------------------------------------------------\n               Creating the " \
          "matrices " \
          "used in LP\n---------------------------------------------------------------------- "

    # --------------------- Generating the Coeff. matrix -------------------

    A = []
    remainder = cols
    for d in dims:
        for r in range(d):
            row = []
            for repeat in range(int(cols / (remainder / d))):
                for i in range(int(remainder / d)):
                    if repeat % d == r:
                        row.append(1)
                    else:
                        row.append(0)
            A.append(row)
        remainder /= d
    return A


# ----------------------------------------------------------------------
#                     Solving the CS using CPLEX
# ----------------------------------------------------------------------

def spgl1_solve(A, b, epsilon, obj):
    W = [w + 1.01 for w in obj]
    x, o1, o2, o3 = spg_bpdn(np.array(A), np.array(b), epsilon, weights=np.array(W), verbosity=1)
    return x, sum(x), np.linalg.norm(np.subtract(np.dot(A, x), np.array(b)), 1)


def solve_cs_qcp(A, b, cols, epsilon, obj):
    ub = [1.0 for x in range(cols)]
    lb = [0.0 for x in range(cols)]

    #---------Constructing matrices of the quadratic constraint---------
    A_t = map(list, zip(*A))
    m = len(A)
    n = len(A_t)
    cnames = ['x' + str(x) for x in range(cols)]

    # Q : quadratic expression, which in our problem is A_t*A
    A_t_A = [[],[],[]] # matrix in index-index-value format used by CPLEX
    for i in range(n):
	row = []
	for j in range(n):
	    val = 0
	    for k in range(m):
		val += A_t[i][k] * A[k][j]
	    row.append(val)
	    A_t_A[0].append(i)
	    A_t_A[1].append(j)
	    A_t_A[2].append(val)

    # L : linear expression, which in our problem is 2(A_t*b)_t
    A_t_b = [[],[]] # vector in index-value format used by CPLEX
    for i in range(n):
	val = 0
	for j in range(m):
	    val += -2 * A_t[i][j] * b[j]
	A_t_b[0].append(i)
	A_t_b[1].append(val)

    # RHS : righthand side, which in our problem is epsilon - ||b||^2
    b_norm2 = 0
    for i in range(len(b)):
	b_norm2 += b[i] * b[i]
    RHS = epsilon - b_norm2

    try:
        prob = cplex.Cplex()
        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.variables.add(obj=obj, ub=ub, lb=lb, names=cnames)
        prob.quadratic_constraints.add(name="error",
                                lin_expr=A_t_b,
                                quad_expr=A_t_A,
                                rhs=RHS,
                                sense="L")

        prob.solve()

        print "Solution status:", prob.solution.get_status(), prob.solution.status[prob.solution.get_status()]

        x = prob.solution.get_values()
	print prob.solution.get_objective_value()
        error = np.sum(np.subtract(np.dot(A, x), b))
        return x, sum(x), error

    except CplexError, exc:
        print exc
        exit(-1)


def solve_cs(A, b, cols, epsilon=0.01, obj=None, qcp=False):
    print "----------------------------------------------------------------------\n                     Solving the CS " \
          "using CPLEX\n---------------------------------------------------------------------- "
    print

    # ----------------------- Preparing the Inputs -------------------------

    if obj == None:
	obj = [1.0 for x in range(cols)]

    if qcp:
	return solve_cs_qcp(A, b, cols, epsilon, obj)

    ub = [1.0 for x in range(cols)]
    lb = [0.0 for x in range(cols)]
    cnames = ['x' + str(x) for x in range(cols)]
    rnames = ['c' + str(x) for x in range(len(A) * 2)]

    rhs = []
    rhs.extend(b)
    for y in b:
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

    # -------------------- Using cplex to solve the LP ----------------------

    try:
        prob = cplex.Cplex()
        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.variables.add(obj=obj, ub=ub, lb=lb, names=cnames)
        prob.linear_constraints.add(lin_expr=constraints, senses=senses, rhs=rhs, names=rnames)

        prob.solve()

        print "Solution status:", prob.solution.get_status(), prob.solution.status[prob.solution.get_status()]

        x = prob.solution.get_values()
	print prob.solution.get_objective_value()
        error = np.sum(np.subtract(np.dot(A, x), b))
        return x, sum(x), error

    except CplexError, exc:
        print exc
        exit(-1)


'''
Predict strains using CS and output in csv file
Input:
    dataPath, absolute path to directory containing samples' alleles and proportions
    pathToDistMat, path to directory containing edit distances matrix for each gene
    refStrains, path to strain_ref.txt
    outputPath, path to output csv file
    loci, list of locus
    objectiveOption, "all" means all objective components and "noPropAndErr" means omitting proportion and error terms
    globalILP_option, "all" if running on all samples
'''


def strainSolver(dataPath, refStrains, outputPath, objectiveOption, globalILP_option='all', timelimit=600, gap=8,
                 loci=["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"], pathToDistMat=None, eps=0.01):
    # ------------------------- Data handling ----------------------------
    # Parameters
    propFormat = 1  # proportion in percentage or fraction
    numLoci = len(loci)

    # read data for samples and reference
    data, numSamples = readData(dataPath, loci, globalILP_option)
    newNameToOriName = dict()
    namingIndex = 1
    for i in sorted(data.keys()):
        newNameToOriName["s{}".format(namingIndex)] = i
        data["s{}".format(namingIndex)] = data.pop(i)
        namingIndex += 1
    reference = pd.read_csv(refStrains, sep="\t", usecols=range(1, numLoci + 1))
    lociNames = list(reference.columns.values)
    numReference = reference.shape[0]
    allSamples = data.keys()

    # check proportions sum to 100
    checkProp(data, propFormat)

    # round the proportions to 3 decimal places
    data = roundProp(data)

    # As reference only contains numbers as entries, add gene name to the variants for better identification
    for name in lociNames:
        reference["%s" % name] = name + "_" + reference["%s" % name].astype(str)

    # Get proportions of variants at different locus for each sample
    varAndProp = returnVarAndProportions(data)

    # Get the combinations at all loci across all samples
    strains, numOfComb = returnCombinationsAndNumComb(data, numLoci, loci)
    uniqueStrains = strains.drop_duplicates(loci)
    uniqueStrains = (uniqueStrains[loci]).reset_index(drop=True)
    uniqueStrains["ST"] = uniqueStrains.index.values + 1  # assign indices for strains or each unique combinations
    strains = strains.merge(uniqueStrains, indicator=True, how="left")  # assign the index to the combinations(as strain data frame contains duplicated rows)
    strains = strains.drop("_merge", 1)


    # -------------------- Creating the weights matrix ----------------------
    alleleSet = []	# set of alleles present in sample
    for l in loci:
        alleleSet.append((l, varAndProp.loc[varAndProp['Locus'] == l]['Variant'].tolist()))
    noc = 1		# noc: number of combinations
    for g in alleleSet:
	noc *= len(g[1])

    # The matrix below is the hamming distance between each pair of a possible strain in the sample
    # and a known strain. Its rows correspond to sample strains and its columns to known strains.
    distMat = [[8 for x in range(len(reference.values.tolist()))] for x in range(noc)]

    # distmat is calculated using an undirect approach
    for i in range(len(reference.values.tolist())):
	rem = noc
        for g in alleleSet:
	    refV = reference.iloc[i][g[0]]
	    for j in range(len(g[1])):
		v = g[1][j]
		if refV == v:
		    for repeat in range(int(noc/rem)):
			offset = repeat * rem
			for index in range(int(rem/len(g[1]))):
			    distMat[index + j * int(rem / len(g[1])) + offset][i] -= 1
	    rem = int(rem / len(g[1]))

    weights = [-1 if min(i) == 0 else min(i) for i in distMat]

    # For each variants, get a mapping of which strains it maps to
    varSampToST = mapVarAndSampleToStrain(strains, loci, allSamples)

    if len(allSamples) == 1:
        print "\n\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n"

        # b is the measurements matrix in cs
        b = []
        # dims is the number of variants in each locus
        dims = []

        for p in varAndProp['Proportion']:
            if p < 1.0:
	        b.append(float(p))


	# below we calculate dims, which is the number of alleles per each locus
        loci = varAndProp['Locus'].tolist()
        prev = ''
        for l in loci:
	    if prev == l:
	        continue
	    prev = l
	    if len(varAndProp.index[varAndProp['Locus'] == l]) > 1:
                dims.append(len(varAndProp.index[varAndProp['Locus'] == l]))

        # A is the coefficient matrix for the cs problem
        samp, numOfC = numOfComb.popitem()
	numOfComb[samp] = numOfC
	# print samp, numOfC
        A = defineProblem(numOfC, dims)

        # Now we solve the cs problem
        # strainProps, sumOfProps, error = solve_cs(A, b, numOfC, eps, weights, True)
        strainProps, sumOfProps, error = spgl1_solve(A, b, eps, weights)
        print '~~~', sumOfProps, error

        strainNums = np.nonzero(strainProps)[0]
        strainNumsIndexes = [i + 1 for i in strainNums]

        # output fromatting
        output = strains.merge(reference, indicator=True, how="left")
        output["_merge"].replace(to_replace="both", value="Existing", inplace=True)
        output["_merge"].replace(to_replace="left_only", value="New", inplace=True)
        output = output.rename(columns = {"_merge":"New/Existing"})
        out = output.drop_duplicates(loci)
        retainCol = ["ST", "New/Existing"]
        out = out[retainCol].reset_index(drop=True)
        out = pd.merge(out, strains.iloc[strainNums], on='ST')
        out['Proportions'] = [strainProps[i] for i in strainNums]
	out_cols = out.columns.tolist()
	out_cols = out_cols[0:-2] + out_cols[-1:] + out_cols[-2:-1]
	out = out[out_cols]
        print out
	out.drop('Sample', axis=1, inplace=True)
        out.to_csv("{0}/{1}_strainsAndProportions.csv".format(outputPath, newNameToOriName[samp]))

        print "\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n\n"
	return out

    else:
	print "error: this simulation is for single sample mode"
	return
