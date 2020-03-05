%% Preamble
% Program: E_Coli_Model_Single_Knockout.m
% Author: Jonathan Myers
% Date: December 8, 2019
% Purpose: Determine gene knockout in E. Coli metabolisim model.
% Arguments: None.
% Loads: Ec_iJO1366.mat
% Calls: None.
% Returns: None.

%% Setup
close all
%clear
%clc

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
singleRemoveReactionGeneList = cell(length(model.genes),1);                       % Initialize reactions removed for gene index.
for numGene = 1:length(model.genes)
    rxnList = find(model.rxnGeneMat(:,numGene)==1);                         % Reactions that use gene index enzyme.
    if (~isempty(rxnList))
        x = true(size(model.genes));                                        % Build good genes matrix.
        x(numGene) = false;                                                 % Knock out gene index.
        singleRemoveReactionGeneList{numGene} = [];                                   % Initialize reactions to remove.
        for i = 1:length(rxnList)
            if (~eval(model.rules{rxnList(i)}))                             % Looks at x vector contents.
                singleRemoveReactionGeneList{numGene} = union(singleRemoveReactionGeneList{numGene},rxnList(i));
            end
        end
    end
end

%% Single Gene Knockout Reactions in Stoichiometry
singleBiomassGeneration = zeros(length(model.genes),3);
singleBiomassGeneration(:,1) = 1:length(model.genes);

parfor i = transpose(singleBiomassGeneration(:,1))
    mod = [];
    param = [];
    rxnS = model.S;
    for j = 1:length(singleRemoveReactionGeneList(i))
        rxnS(:,singleRemoveReactionGeneList{i}) = 0;
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
    singleBiomassGeneration(i,2) = result.objval;
end
singleBiomassGeneration(:,3) = singleBiomassGeneration(:,2) >= 0.5;

% M02 End Program