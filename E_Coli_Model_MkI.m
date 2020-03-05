%% Preamble
% Program: E_Coli_Model_MkI.m
% Author: Jonathan Myers
% Date: December 8, 2019
% Purpose: Determine gene knockout in E. Coli metabolisim model.
% Arguments: None.
% Loads: Ec_iJO1366.mat
% Calls: None.
% Returns: None.

%% Setup
close all
clear
clc

%% Single Cell Without Knockout
load Ec_iJO1366
mod.sense = '=';
mod.modelsense = 'max';
mod.A = model.S;
mod.obj = model.c;
mod.rhs = model.b;
mod.ub = model.ub;
mod.lb = model.lb;
param.OutputFlag = 0;
param.Method = 1;
result = gurobi(mod,param);
factor = result.objval;

%% Determine Single Gene-Reaction Knockouts
removeReactionGeneList = cell(length(model.genes),1);                       % Initialize reactions removed for gene index.
for numGene = 1:length(model.genes)
    rxnList = find(model.rxnGeneMat(:,numGene)==1);                         % Reactions that use gene index enzyme.
    if (~isempty(rxnList))
        x = true(size(model.genes));                                        % Build good genes matrix.
        x(numGene) = false;                                                 % Knock out gene index.
        removeReactionGeneList{numGene} = [];                                   % Initialize reactions to remove.
        for i = 1:length(rxnList)
            if (~eval(model.rules{rxnList(i)}))                             % Looks at x vector contents.
                removeReactionGeneList{numGene} = union(removeReactionGeneList{numGene},rxnList(i));
            end
        end
    end
end

%% Single Gene Knockout Reactions in Stoichiometry
biomassGeneration = zeros(length(model.genes),3);
biomassGeneration(:,1) = 1:length(model.genes);

parfor i = transpose(biomassGeneration(:,1))
    mod = [];
    param = [];
    rxnS = model.S;
    for j = 1:length(removeReactionGeneList(i))
        rxnS(:,removeReactionGeneList{i}) = 0;
    end
    mod.sense = '=';
    mod.modelsense = 'max';
    mod.A = rxnS;
    mod.obj = model.c;
    mod.rhs = model.b;
    mod.ub = model.ub;
    mod.lb = model.lb;
    param.OutputFlag = 0;
    param.Method = 1;
    result = gurobi(mod,param);
    biomassGeneration(i,2) = result.objval;
end
biomassGeneration(:,3) = biomassGeneration(:,2) >= 0.5;

%% Determine Double Gene-Reaction Knockouts
% Best Optimization Strategy, but too time intensive.
% nonLethalGene = biomassGeneration(biomassGeneration(:,3) == 1,1:2);
% knockoutNumber = 2;
% genePair = combnk(nonLethalGene(:,1),knockoutNumber);
% genePairNumber = size(genePair,1);
% removeReactionGeneList = cell(genePairNumber,2);           % Initialize reactions removed for gene index.
% for i = 1:size(removeReactionGeneList,1)
%     removeReactionGeneList{i,2} = genePair(i,:);
% end

totalNonLethalGene = biomassGeneration(biomassGeneration(:,3) == 1,1:2);    % Array of all nonlethal genes.
nonLethalGeneIndices = totalNonLethalGene(:,2) == max(totalNonLethalGene(:,2)); % Array of best biomass flux.
nonLethalGene = num2cell(totalNonLethalGene(nonLethalGeneIndices,:));       % Cell array of best nonlethal genes.
knockoutNumber = 2;     % Number of knockouts.

genePair = zeros((size(model.genes,1))*size(nonLethalGene,1),knockoutNumber); % Array of all genes knockouts.
for i = 0:size(nonLethalGene,1)-1
    for j = 1:size(model.genes,1)
        genePair(j+(i*(size(model.genes,1))),1) = nonLethalGene{i+1,1};
    end
end
for i = 2:knockoutNumber
    for j = 0:size(nonLethalGene,1)-1
        for k = 1:size(model.genes,1)
            genePair(k+(j*(size(model.genes,1))),i) = k;
        end
    end
end

genePairNumber = size(genePair,1);
removeReactionGeneList = cell(genePairNumber,2);           % Initialize reactions removed for gene index.
for i = 1:size(removeReactionGeneList,1)
    removeReactionGeneList{i,2} = genePair(i,:);
end

for knockoutIncrement = 1:genePairNumber
    geneVector = genePair(knockoutIncrement,:);
    rxnList = [];
        for i=1:length(geneVector)
            rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1)); % Reactions that use gene index enzymes.
        end
    if (~isempty(rxnList))
        x = true(size(model.genes,1));                                      % Build total genes matrix.
        x(geneVector) = false;                                              % Knock out gene indices.
        removeReactionGeneList{knockoutIncrement,1} = [];                     % Initialize reactions to remove.
        for i = 1:length(rxnList)
            if (~eval(model.rules{rxnList(i)}))                             % Use determine if x contents 
                removeReactionGeneList{knockoutIncrement,1} = union(removeReactionGeneList{knockoutIncrement,1},rxnList(i));
            end
        end
    end
end

%% Double Gene Knockout Reactions in Stoichiometry
biomassGeneration = zeros(genePairNumber,3);
biomassGeneration(:,1) = 1:genePairNumber;

parfor i = transpose(biomassGeneration(:,1))
    mod = [];
    param = [];
    rxnS = model.S;
    for j = 1:length(removeReactionGeneList(i))
        rxnS(:,removeReactionGeneList{i}) = 0;
    end
    mod.sense = '=';
    mod.modelsense = 'max';
    mod.A = rxnS;
    mod.obj = model.c;
    mod.rhs = model.b;
    mod.ub = model.ub;
    mod.lb = model.lb;
    param.OutputFlag = 0;
    param.Method = 1;
    result = gurobi(mod,param);
    biomassGeneration(i,2) = result.objval;
end
biomassGeneration(:,3) = biomassGeneration(:,2) >= 0.5;

%% Multiple Gene-Reaction Knockout and Stoichiometry Optimization


% M02 End Program