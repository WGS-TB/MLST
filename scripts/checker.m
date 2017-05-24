clear;
clc;
format long;

numLoci = 3;
startingSampleNum = 1;
numSamples = 2;
loci=['clpA';'clpX';'nifS'];

%read files output by the python file
proportion = readtable('proportion.csv', 'Delimiter',',');
strain = readtable('strain.csv');
error = readtable('error.csv');
propConstrRHS = importdata('propConstrRHS.csv');
propConstrRHS = propConstrRHS(2:end,2);

%Separate piDotComb.csv into cells
fid = fopen('piDotComb.csv');
textLine = fgets(fid); % Read first line.
lineCounter = 1;
while ischar(textLine)
	% get into array
	row = sscanf(textLine, '%s');
    row = strsplit(row, ',');
	% Alternate way where the whole array is in one cell.
    
    for i=1:length(row)
        piDotComb(lineCounter, i) = row(i);
    end
    
	% Read the next line.
    textLine = fgets(fid);
	lineCounter = lineCounter + 1;
end
fclose(fid);
piDotComb = piDotComb(2:end, 2:end);

%getting the possible subsets that can explain the data
strainSubset = strainSubsets(strain, proportion, numSamples, loci);

fval_min = Inf;
x_min = [];
var_min = [];
subset_min = [];
sol = [];
num_variable = [];
%Here goes the solving
for problem=1:size(strainSubset,1)
    mySubset = strainSubset(problem);
    
    %error sum to 0 constraint
    [ Aeq_err, beq_err, err, lb_err, ub_err ] = errSumTo0(error, numSamples, numLoci);
    
    %proportion decision variables sum to 1
    [ Aeq_prop, beq_prop, prop, lb_prop, ub_prop ] = propSumTo1( mySubset, proportion, numSamples );
    
    %piDotCombination constraint
    [ Aeq_piDotComb, beq_piDotComb ] = piDotCombConstr( prop, piDotComb, propConstrRHS );
    
    %error bound constraint
    [ A_errB1, b_errB1, A_errB2, b_errB2 ] = errorBoundConstr(mySubset, error, proportion );
    
    variables = [err; prop];
    numVariables = length(variables);
    
    for i=1:numVariables
        eval([variables{i},' = ', num2str(i),';']); 
    end
    
    %set lower bound of variables
    lb = zeros(size(variables));
    exist_err = ismember(variables, err);
    index_err = find(exist_err==1);
    exist_prop = ismember(variables, prop);
    index_prop = find(exist_prop==1);
    lb([index_err; index_prop]) = [lb_err; lb_prop];
    
    %set upper bound of variables
    ub=zeros(size(variables));
    ub([index_err; index_prop]) = [ub_err; ub_prop];
    
    %error sum to 0 and prop sum to 1
    Aeq = blkdiag(Aeq_err, Aeq_prop);
    beq = [beq_err; beq_prop];
    
    %piDotComb constraint
    Aeq = [Aeq; zeros(size(err,1)) Aeq_piDotComb];
    beq = [beq; beq_piDotComb];
    
    %error bound constraint
    A = [A_errB1; A_errB2];
    b = [b_errB1; b_errB2];
    
    %objective function
    weightOfProp = proportion{ismember(proportion{:,'DecisionVariable'}, prop), 'Weights'};
    f = [ones(size(err,1),1)' weightOfProp'];
    strainBool = ismember(proportion{:,'ST'}, mySubset);
    remainingStrain = proportion(strainBool, :);
    [uniqueStrain, strainIndex] = unique(remainingStrain{:, 'ST'});
    weightOfStrain = remainingStrain{strainIndex, 'Weights'};
    totalStrainUsage = sum(weightOfStrain);
    numStrains = length(mySubset);
    
    options = optimoptions('linprog','Algorithm','dual-simplex');
	[x, fval, exitflag, output] = linprog(f,A,b,Aeq,beq,lb,ub,options);
    
    if exitflag == 1 %if feasible
         nonZeroPropIndex = [];
        %get the index of non zero prop variables
        for k=size(err,1)+1:size(variables,1)
            if x(k) ~= 0
                nonZeroPropIndex = [nonZeroPropIndex; k];
            end
        end
    
        nonZeroProp = variables(nonZeroPropIndex);
        matchStrain = ismember(proportion{:,'DecisionVariable'},nonZeroProp);
        %get the unique strain these prop variables referring to
        nonZeroStrain = unique(proportion{matchStrain, 'ST'});
        
        %run update iff the number of strains in the subset is the same as
        %the strain which has non zero prop var
        if size(nonZeroStrain,1) == size(mySubset,2)
            fval = fval + totalStrainUsage;
            if fval <= fval_min
                sol = [sol; fval];
                fval_min = fval;
                x_min = [x_min; x];
                num_variable = [num_variable; size(x,1)];
                var_min = [var_min; variables];
                subset_min = [subset_min; mySubset];
            end
        end
    end
end

%Print the optimal value and values of decision variable
start=0;
for num_min=1:size(sol,1)
    numberOfVariable = num_variable(num_min);
    for i = 1:numberOfVariable
        fprintf('%12.5f \t%s\n',x_min(start+i),var_min{start+i}) 
    end
    fprintf('\n')
    start = start+i;
end

fval_min
subset_min
sol

%% ===================================== Functions Definition =================================== %%
function [ st ] = strainSubsets( strain, proportion, numSamples, loci )
%Return a container.Map where the key=subset index and value=the strains
%that the subset contains
numOfStrains = length(strain{:, 'ST'});
subsetIndices = dec2bin(1:2^numOfStrains - 1) - '0';    %0/1 representation of a subset

numOfSub = length(subsetIndices);
totalVariants = totalNumVar(proportion, loci);

st = containers.Map('KeyType','int64','ValueType','any');

track=1;    %track subset index
for i=1:numOfSub
    sub = subsetIndices(i,:);
    temp = [];

    %For those position having 1, grab the particular ST index from the data we loaded
    for j=1:length(sub)
        if sub(j) == 1
            temp = [temp strain{j, 'ST'}];  
        end
    end

    %Add the subset into consideration iff it covers all the samples and
    %all variants observed at all samples are explained
    if numSampExplain(temp, proportion) == numSamples && numVarExplain(temp, proportion, loci) >= totalVariants
        st(track) = temp;
        track = track+1;
    end
end

end


function [ numUniqueSamp ] = numSampExplain( subset, proportion )
%Return the number of samples cover by a given subset of strains
match = ismember(proportion{:,'ST'}, subset);
samples = proportion{match, 'Sample'};
uniqueSamp = unique(samples);
numUniqueSamp = length(uniqueSamp);

end


function [ varExplain ] = numVarExplain( subset, proportion, loci )
%This function returns the number of variants explainable by the given
%subset of strains(there could be duplicated variants identified at
%different samples)
samples = unique(proportion{:,'Sample'});
%get details on the subset of strains
match = ismember(proportion{:,'ST'}, subset);
%get rid of proportion decision variables which are unrelated to the subset
remaining = proportion(match, :);
varExplain = 0;

%count the number of variants in remaining portion
for i=1:length(samples)
    oneSample = ismember(remaining{:, 'Sample'}, samples{i});
    temp = remaining(oneSample, :);
    
    for j=1:size(loci,1)
        tempLoc = temp{:,loci(j,:)};
        tempLoc = unique(tempLoc);
        varExplain = varExplain + length(tempLoc);
    end
end

end


function [ totalVar ] = totalNumVar( proportion, loci )
%This function returns the number of variants in all sample(there could be
%duplications)
samples = unique(proportion{:,'Sample'});
totalVar = 0;

for i=1:length(samples)
    %get all the variants for a single sample
    oneSample = ismember(proportion{:, 'Sample'}, samples{i});
    temp = proportion(oneSample, :);
    
    %iterate over each locus
    for j=1:size(loci,1)
        tempLoc = temp{:,loci(j,:)};
        tempLoc = unique(tempLoc);
        totalVar = totalVar + length(tempLoc);
    end
end


end


function [ Aeq, beq, err, lb, ub ] = errSumTo0( error, numSamples, numLoci )
%Return matrix constraints related to the constraint where error sum to 0
Aeq = zeros(numLoci*numSamples, length(error{:,'DecisionVariable'}));

row=1;col=1;
lb=[]; ub=[];

%Fill in matrix and vectors
while 1
    Aeq(row, col) = 1;
    currSamp = error{col, 'Sample'};
    currLocus = error{col, 'Locus'};
    lb = [lb; -error{col,'Proportion'}];
    ub = [ub; 1-error{col,'Proportion'}];
    
    %stop loop after getting all error variables
    if col == length(error{:,'DecisionVariable'})
        break;
    end
    
    col = col+1;
    nextSamp = error{col, 'Sample'};
    nextLocus = error{col, 'Locus'};
    
    %if next sample or next locus isn't the same as current one, go to next
    %row
    if ~(strcmp(currSamp, nextSamp) && strcmp(currLocus, nextLocus))
        row=row+1;
    end

end

beq = zeros(size(Aeq,1), 1);
err = error{:,'DecisionVariable'};


end


function [ Aeq, beq, prop, lb, ub ] = propSumTo1( strains, proportion, numSamples )
%return matrix constraints, Aeq and beq related to the constraints where
%for each sample, the sum of the proportion decision variables =1
%prop, lb, ub are the prop decision variable names, lower bound and upper
%bound for them respectively

%get prop dec var related to this subset
remain = ismember(proportion{:,'ST'}, strains);
propVar = proportion{remain, {'Sample','DecisionVariable'}};

Aeq = zeros(numSamples, length(propVar));

row=1;col=1;

%fill in the matrix
while 1
    Aeq(row, col) = 1;
    currSamp = propVar(col,1);
    
    if col == length(propVar)
        break;
    end
    
    col = col+1;
    nextSamp = propVar(col,1);
    
    %change row when next sample is different
    if ~(strcmp(currSamp, nextSamp))
        row=row+1;
    end
end

beq = ones(numSamples,1);
lb = zeros(length(propVar),1);
ub = ones(length(propVar),1);
prop = propVar(:,2);

end


function [ Aeq, beq ] = piDotCombConstr( remainPropVar, piDotComb, propConstrRHS)
%Return matrix and vector constraints related to the preservation of
%proportions in each sample(the constraint which has the sum of pi * V in
%the paper). Here remainPropVar is the prop dec var from function
%propSumTo1

numRows = size(piDotComb,1);
for row=1:numRows
    %As each row of piDotComb has different lengths, here is to grab those
    %non-empty cells
    notEmpty = ~cellfun(@isempty, piDotComb(row,:));
    temp = piDotComb(row,notEmpty);
    %grab those related proportion decision variables
    exist = ismember(temp, remainPropVar);
    piDotCombLHS{row} = temp(exist);
end

%There may be chances that "temp" in the for loop above does not contain
%any common prop dec var
cellIsNotEmpty = ~cellfun(@isempty, piDotCombLHS);
piDotCombLHS = piDotCombLHS(cellIsNotEmpty);
piDotCombRHS = propConstrRHS(cellIsNotEmpty);

Aeq = [];
beq = piDotCombRHS;

%According to the order of the given prop dec var, fill in matrix A
%accordingly
for rowTrack=1:size(beq,1)
    temp = piDotCombLHS{rowTrack};
    exist = ismember(remainPropVar, temp);
    
    Aeq = [Aeq; double(exist)'];
end

end


function [ A_1, b_1, A_2, b_2 ] = errorBoundConstr(strains, error, proportion )
%get only relevant proportion var
remain = ismember(proportion{:,'ST'}, strains);
remainPropVar = proportion(remain, :);

%error <= sum of prop variable - proportion and error <= proportion - sum
%of prop variable
error_variable = error{:,'DecisionVariable'};
errorVarMatrix = eye(length(error_variable));

allPropVariables = remainPropVar{:,'DecisionVariable'};

propVarMat_errLessSumMinProp = zeros(length(error_variable), length(allPropVariables));
propVarMat_errLessPropMinSum = zeros(length(error_variable), length(allPropVariables));

for i=1:size(error, 1)
    samp = error{i, 'Sample'};
    variant = error{i,'Variant'};
    locus = error{i,'Locus'};
    
    %unique identifier of error variable: sample and variant at that locus
    temp = ismember(remainPropVar{:,'Sample'}, samp) & ismember(remainPropVar{:,locus}, variant);
    tempPropVar = remainPropVar{temp, 'DecisionVariable'};
    
    propVarIsIn = ismember(allPropVariables, tempPropVar);
    index = [];
    
    %get the index where the proportion variable is related to the error
    for j=1:length(propVarIsIn)
        if propVarIsIn(j) ==1
            index = [index j];
        end
    end
    
    propVarMat_errLessSumMinProp(i, index) = -1;
    propVarMat_errLessPropMinSum(i, index) = 1;
end

A_1 = [errorVarMatrix propVarMat_errLessSumMinProp];
A_2 = [errorVarMatrix propVarMat_errLessPropMinSum];

b_1 = -error{:, 'Proportion'};
b_2 = error{:, 'Proportion'};

end

