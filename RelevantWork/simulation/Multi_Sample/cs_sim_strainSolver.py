from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from cplex.exceptions import CplexError
from pipeline_functions import *
from scipy.spatial.distance import hamming
from spgl1 import spgl1, spg_mmv, spgSetParms

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
            for repeat in xrange(int(cols / (remainder / d))):
                for i in xrange(int(remainder / d)):
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


def solve_cs(A, b, cols, epsilon=0.01, obj=None):
    print "----------------------------------------------------------------------\n                     Solving the CS " \
          "using CPLEX\n---------------------------------------------------------------------- "
    print

    # ----------------------- Preparing the Inputs -------------------------

    if obj == None:
	obj = [1.0 for x in range(cols)]
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

    # Proportions should add up to 1
    '''rhs.append(1 + 10 * epsilon)
    rhs.append(-1 + 10 * epsilon)
    constraints.append([range(cols), [1 for i in row]])
    constraints.append([range(cols), [-1 for i in row]])'''

    #print
    #for row in constraints:
	#print row, '     ', rhs[constraints.index(row)]
    #print

    senses = ""
    for i in range(len(constraints)):
        senses += "L"

    # print "    Constraints:"
    # for r in range(len(constraints)):
    #     print constraints[r][1], "\t\t\t", rhs[r]
    # print "\n"

    # -------------------- Using cplex to solve the LP ----------------------

    try:
        prob = cplex.Cplex()
        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.variables.add(obj=obj, ub=ub, lb=lb, names=cnames)
        prob.linear_constraints.add(lin_expr=constraints, senses=senses, rhs=rhs, names=rnames)

        prob.solve()

        # print
        print "Solution status:", prob.solution.get_status(), prob.solution.status[prob.solution.get_status()]
        # print "Solution value  = ", prob.solution.get_objective_value()
        # for i in range(prob.variables.get_num()):
        #     print prob.variables.get_names()[i], '\t=\t', prob.solution.get_values()[i]

        x = prob.solution.get_values()
	print prob.solution.get_objective_value()
        error = np.sum(numpy.subtract(numpy.dot(A, x), b))
        return x, sum(x), error

    except CplexError, exc:
        print exc
        exit(-1)

def createMeasureMatrix(strains, loci, variants):
	# todo: create measurment matrix from strains and variants
	# strains is a dataframe, variants is a list
	return

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
    print allSamples

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

    # For each variants, get a mapping of which strains it maps to
    varSampToST = mapVarAndSampleToStrain(strains, loci, allSamples)

    if len(allSamples) == 1:
        print "\n\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n"

	alleleSet = []
    	for l in loci:
            alleleSet.append((l, varAndProp.loc[varAndProp['Locus'] == l]['Variant'].tolist()))
    	noc = 1
    	for g in alleleSet:
	    noc *= len(g[1])
   
    	distMat = [[8 for x in range(len(reference.values.tolist()))] for x in range(noc)]

    	# Calculating weights for minimization
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

        # b is the measurements matrix in cs
        b = []
        # dims is the number of variants in each locus
        dims = []

        for p in varAndProp['Proportion']:
            if p < 1.0:
	        b.append(float(p))
	    

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
	print samp
        A = defineProblem(numOfC, dims)

        # Now we solve the cs problem
        strainProps, sumOfProps, error = solve_cs(A, b, numOfC, eps, weights)
        print '~~~', sumOfProps, error

        strainNums = np.nonzero(strainProps)[0]

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

	print "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n\n"

	print strains
	
	# Solving for strains present in several samples
	# some data processing
	vpdf = varAndProp.drop_duplicates(['Variant'])
	loci = vpdf['Locus'].tolist()
	from collections import OrderedDict
	loci = list(OrderedDict.fromkeys(loci))

	strSet = set(strains.drop_duplicates(loci)['ST'].tolist())
	chosenStrs = list()
	for s in strSet:
	    sampSet = set(strains.loc[strains['ST'] == s].drop_duplicates('Sample')['Sample'].tolist())
	    if len(sampSet) > 1:
		chosenStrs.append(s - 1)
	sharedStrs = strains.iloc[chosenStrs].reset_index(drop=True)
	print 'SHARED:\n', sharedStrs
	chosenVars = list()
	dims = []
	variants = []
	for l in loci:
	    dims.append(len(set(sharedStrs[l].tolist())))
	    variants.append(list(set(sharedStrs[l].tolist())))
	A = createMeasureMatrix(sharedStrs, loci, variants)
	
	
        sampAlleles = [] # List for having each samples data

	nrows = len((varAndProp.drop_duplicates(['Variant'])).values.tolist()) # Number of rows of A and rows of B
	nsamp = 0 # Number of columns of X and columns of B

	for samp, numOfC in numOfComb.iteritems():
	    sampAlleles.append(varAndProp.loc[varAndProp['Sample'] == samp])
	    nsamp += 1
	
	# Defining the problem
	B = [[0 for x in range(nsamp)] for x in range(nrows)]
	
	
	dims = []
        for l in loci:
	    # if len(vpdf.index[vpdf['Locus'] == l]) > 1:
            dims.append(len(vpdf.index[vpdf['Locus'] == l]))
	ncols = 1
	for d in dims:
	    ncols *= d
	A = defineProblem(ncols, dims)

	for ind in range(nsamp):
	    for indV in range(len(vpdf['Variant'].tolist())):
		if vpdf['Variant'].tolist()[indV] in sampAlleles[ind]['Variant'].tolist():
		    B[indV][ind] = sampAlleles[ind].loc[sampAlleles[ind]['Variant'] == vpdf['Variant'].tolist()[indV]]['Proportion'].values.tolist()[0]
	
	# finding all combinations across samples
	vlist = list()
	for l in loci:
	    vlist.append(sorted(list(set(varAndProp.loc[varAndProp['Locus'] == l]['Variant'].tolist()))))
	combination = itertools.product(*vlist)
	combinationIndices = [list(comb) for comb in itertools.product(*[range(len(var)) for var in vlist])]
	allCombs = list()
	for strain,strainIndex in itertools.izip(combination,combinationIndices):
            temp = list(strain)
            allCombs.append(temp)
	allCombs = pd.DataFrame(allCombs, columns=(loci))
	allCombs["ST"] = allCombs.index.values + 1

	    
	# Solving the problem using spgl1	
	print "----------------------------------------------------------------------\n           Solving the Joint Sparsity " \
          "using SPGL1\n---------------------------------------------------------------------- "
        print
	opts = spgSetParms({'verbosity': 1})
	X, _, _, _ = spg_mmv(np.array(A), np.array(B), eps, opts)
	for Row in X.transpose():
	    print Row
	print
	
	# output formatting
	output = allCombs.merge(reference, indicator=True, how="left")
        output["_merge"].replace(to_replace="both", value="Existing", inplace=True)
        output["_merge"].replace(to_replace="left_only", value="New", inplace=True)
        output = output.rename(columns = {"_merge":"New/Existing"})
        out = output.drop_duplicates(loci)
        retainCol = ["ST", "New/Existing"]
	out = out[retainCol].reset_index(drop=True)
	samplecnt = 1
	outputdict = dict()
	for sampStrains in X.transpose():
	    persampout = out.copy(deep=True)
	    strainNums = np.nonzero(sampStrains)[0]
	    persampout = pd.merge(persampout, allCombs.iloc[strainNums], on='ST')
            persampout['Proportions'] = [sampStrains[i] for i in strainNums]
	    persampout['Sample'] = 's' + str(samplecnt)
	    samplecnt += 1
	    print persampout
	    outputdict['s' + str(samplecnt)] = persampout

	print "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n\n"
	return outputdict


