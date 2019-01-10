from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cplex
from cplex.exceptions import CplexError
from pipeline_functions import *
from utility_functions import *
from scipy.spatial.distance import hamming
from spgl1 import spgl1, spg_mmv, spgSetParms, spg_bpdn

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
        print "\nERROR: cplex could not solve the single sample problem with changed proportions of previous step\n"
	return [0 for _ in range(len(b))], 0, np.inf


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
	print "\nERROR: cplex could not solve the single sample problem with changed proportions of previous step\n"
	return [0 for _ in range(len(b))], 0, np.inf

'''
Creating the measurment matrix (A) via inspecting the strains
Input:
    strains, the strains to be inspected
    loci, list of locus
    variant, list of lists of variants in each locus
    dims, list of lists of number of variants in each locus
'''
def createMeasureMatrix(strains, loci, variants, dims):
    AT = [] 				# the matrix created here is the transpose of the matrix we want
    for index, strain in strains.iterrows():
	strainCol = []
	for li in range(len(loci)): 	#li is locus index
	    l = loci[li]
	    for vi in range(dims[li]): 	#vi is variant index
		v = variants[li][vi]
		if strain[l] == v:
		    strainCol.append(1)
		else:
		    strainCol.append(0)
	AT.append(strainCol)	
    A = (np.array(AT)).transpose()
    return A

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

def single_sample_solver(loci, varAndProp, reference, strains, eps, re=False, br=None, vm=None):
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
    bcheck = []
    for p in varAndProp['Proportion']:
        if p < 1.0:
	    b.append(float(p))
	bcheck.append(p)
    
    if re:
	ind1 = 0
	ind2 = 0
	for v in varAndProp['Variant']:
	    if v in vm and bcheck[ind2] != 1.0:
		b[ind1] -= br[vm.index(v)]
		ind1 += 1
	    ind2 += 1
	
    # dims is the number of variants in each locus
    dims = []
    prev = ''
    for l in loci:
	if prev == l:
	    continue
	prev = l
	if len(varAndProp.index[varAndProp['Locus'] == l]) > 1:
            dims.append(len(varAndProp.index[varAndProp['Locus'] == l]))

    # A is the coefficient matrix for the cs problem
    A = defineProblem(noc, dims)

    print 'Proportions:'
    for p in b:
	print p
    print 'varAndProp input'
    print varAndProp
    print np.array(A).shape

    # Now we solve the cs problem
    strainProps, sumOfProps, error = solve_cs(A, b, noc, eps, weights, True)
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
    print "\n**SINGLE SAMPLE OUTPUT:"
    print out
    print "\n-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n\n"
    return out

	

def strainSolver(dataPath, refStrains, outputPath, objectiveOption, globalILP_option='all', timelimit=600, gap=8,
                 loci=["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"], pathToDistMat=None, eps=0.05):
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

    # For each variants, get a mapping of which strains it maps to
    varSampToST = mapVarAndSampleToStrain(strains, loci, allSamples)


    print "\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n\n"
	
    # some data processing
    vpdf = varAndProp.drop_duplicates(['Variant'])
    loci = vpdf['Locus'].tolist()
    from collections import OrderedDict
    loci = list(OrderedDict.fromkeys(loci))
    sampAlleles = [] 	# List for having each samples data
    nsamp = 0 		# Number of columns of X and columns of B
    sampNames = []
    for samp, numOfC in numOfComb.iteritems():
	sampNames.append(int(samp[1:]))
	nsamp += 1
    sampNames = sorted(sampNames)
    print sampNames
    for samp in sampNames:
	sampAlleles.append(varAndProp.loc[varAndProp['Sample'] == 's' + str(samp)])

    for sa in sampAlleles:
        print sa
    print '\n'



    # ~~~~~~~~~~~~  Solving for strains present in several samples  ~~~~~~~~~~~~~~

    # Below we detect the shared strains between all samples
    strSet = set(strains.drop_duplicates(loci)['ST'].tolist())
    chosenStrs = list()
    sampsets = list()
    for s in strSet:
	# those strains which are present in more than two samples are detected by 
	# counting the number of samples they appear in
	sampSet = set(strains.loc[strains['ST'] == s].drop_duplicates('Sample')['Sample'].tolist())
	if len(sampSet) > 1:
	    chosenStrs.append(s - 1)
	    sampsets.append(sampSet)
    sharedStrs = strains.iloc[chosenStrs].reset_index(drop=True)
    sharedStrs['Sample'] = sampsets
    sharedStrs = sharedStrs.drop_duplicates(loci)
    print '**SHARED STRAINS'
    print sharedStrs, '\n'

    # detecting the samples which have a shared strain
    samplesSharing = []
    for ss in sampsets:
	for s in ss:
	    samplesSharing.append(int(s[1:]))
    samplesSharing = sorted(list(set(samplesSharing)))
    nsampS = len(samplesSharing)

    # dictionary to hold output dataframes for each sample
    outputdict = dict()
    # matrix of residues of shared calculation allele proportions
    BR = np.array([[0]])
    variantsMixed = []

    # In some iterations it's possible not to have any shared strains
    if len(chosenStrs) > 0:
        chosenVars = list()
        dims = []
        variants = [[],[]]
        for l in loci:
            dims.append(len(set(sharedStrs[l].tolist())))
            locusvars = list(set(sharedStrs[l].tolist()))
            variants[0].append(locusvars)
	    variants[1].append((l, locusvars))
            variantsMixed.extend(locusvars)

        # A: Measurement matrix
        A = createMeasureMatrix(sharedStrs, loci, variants[0], dims)
	

        # B: Observation Matrix		
        nrows = len(variantsMixed) # Number of rows of A and rows of B
        B = [[0 for x in range(nsampS)] for x in range(nrows)]
        for ind in range(nsampS):
            for indV in range(nrows):
	        if variantsMixed[indV] in sampAlleles[samplesSharing[ind] - 1]['Variant'].tolist():
	            B[indV][ind] = sampAlleles[samplesSharing[ind] - 1].loc[sampAlleles[samplesSharing[ind] - 1]['Variant'] == variantsMixed[indV]]['Proportion'].values.tolist()[0]
	


        # Now we solve the shared strains with mmv, then the remaining strains for each sample with single sample solution
        # For solving the shared strains problem,  we choose a range  of sigmas, and use each of them as  input to get the
        # best result.

        # calculating weights for shared strains
	weights = []
	for st in sharedStrs[loci].values.tolist():
	    weights.append(min([hamming(st, rst) for rst in reference[loci].values.tolist()]) * 8 + 1)
    
        print "----------------------------------------------------------------------\n           Solving the Joint Sparsity " \
        "using SPGL1\n---------------------------------------------------------------------- "
        print

        sigmas = []
        errs = []
        min_err = 100.0
        best_sig = 0
        X1 = np.array([[0]])
	# in order to find the best solution, we solve the join sparsity problem with different values of sigma
        # to find the result with least error and then report that result.
        for i in range(11):
            sigmas.append(i * 0.05)
        for sig in sigmas:
            opts = {'weights' : np.array(weights), 'verbosity' : 1}
            X, _, _, _ = spg_mmv(A, np.array(B), sig, opts)
            errmat = np.dot(A,X) - np.array(B)
            err = np.sqrt(sum([np.linalg.norm(row) for row in errmat]))
            errs.append(err)
    	if min_err > err:
    	    min_err = err
    	    best_sig = sig
    	    X1 = X
        print '**BEST SIG AND ITS ERR ', best_sig, min_err
	print '**RESULT column=samp row=strain:'
	print pd.DataFrame(X1)
	
	# setting the residue matrix for samples that were included in this part
	BR = np.array(B) - np.dot(A,X1)

	# creating the output of shared strains. Latter output for individual samples will be appended to these.
	output = sharedStrs.merge(reference, indicator=True, how="left")
        output["_merge"].replace(to_replace="both", value="Existing", inplace=True)
        output["_merge"].replace(to_replace="left_only", value="New", inplace=True)
        output = output.rename(columns = {"_merge":"New/Existing"})
        out = output.drop_duplicates(loci)
	
	# more output formatting and adding the results to proper dataframes
        retainCol = ["ST", "New/Existing"]
        out = out[retainCol].reset_index(drop=True)
        samplecnt = 0
        for sampStrains in X1.transpose():
	    sampName = 's' + str(samplesSharing[samplecnt])
            persampout = out.copy(deep=True)
            strainNums = np.nonzero(sampStrains)[0]
	    trueStrainNums = []
	    for i in strainNums:
		if sampName in sharedStrs.iloc[i]['Sample']:
		    trueStrainNums.append(i)
            persampout = pd.merge(persampout, sharedStrs.iloc[trueStrainNums], on='ST')
            persampout['Proportions'] = [sampStrains[i] for i in trueStrainNums]
            persampout['Sample'] = sampName
	    columns = persampout.columns.tolist()
	    columns = columns[:-2] + columns[-1:] + columns[-2:-1]
	    persampout = persampout[columns]
            #print persampout
            outputdict['s' + str(samplesSharing[samplecnt])] = persampout
	    samplecnt += 1

    # ------------------ Per-sample solution using single sample solver -------------------
    # after subtraction of the results of shared strains, the remaining proportions are pas
    # -sed to the function which solves the single sample problem.
    for i in range(nsamp):
	sampName = 's' + str(i + 1)
	if i + 1 in samplesSharing:
	    persampout = single_sample_solver(loci, sampAlleles[i], reference, strains.loc[strains['Sample'] == sampName], eps, True, BR.transpose().tolist()[samplesSharing.index(i + 1)], variantsMixed)
	    outputdict[sampName] = outputdict[sampName].append(persampout)
	else:
            persampout = single_sample_solver(loci, sampAlleles[i], reference, strains.loc[strains['Sample'] == sampName], eps)
	    outputdict[sampName] = persampout
	print outputdict[sampName]
    return outputdict

