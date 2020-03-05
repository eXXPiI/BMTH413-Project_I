%% Preamble
% Program: E_Coli_Model_MkII.m
% Author: Jonathan Myers
% Date: December 15, 2019
% Purpose: Determine gene knockout in E. Coli metabolisim model.
% Arguments: None.
% Loads: Ec_iJO1366.mat
% Calls: None.
% Returns: None.

%% Setup
close all
clear
clc
load Ec_iJO1366

%% Single Gene-Reaction Knockout and Stoichiometry Optimization
genomeSize = size(model.genes,1);
viableKnockoutGenes = transpose(1:genomeSize);
lenViableKnockoutGenes = size(viableKnockoutGenes,1);
trackingKnockout = cell(lenViableKnockoutGenes,5);
trackingKnockout(:,1) = num2cell(1:lenViableKnockoutGenes);
for numGene = 1:genomeSize
    rxnList = find(model.rxnGeneMat(:,numGene)==1);                         % Reactions that use gene index enzyme.
    if (~isempty(rxnList))
        trackingKnockout{numGene,3} = numGene;
        x = true(genomeSize);                                        % Build good genes matrix.
        x(numGene) = false; %#ok<NASGU>                                     % Knock out gene index.
        trackingKnockout{numGene,2} = [];                                   % Initialize reactions to remove.
        for i = 1:length(rxnList)
            if (~eval(model.rules{rxnList(i)}))                             % Looks at x vector contents.
                trackingKnockout{numGene,2} = union(trackingKnockout{numGene,2},rxnList(i));
            end
        end
    end
end

biomassGeneration = zeros(lenViableKnockoutGenes,1);
parfor i = transpose(cell2mat(trackingKnockout(:,1)))
    mod = [];
    param = [];
    rxnS = model.S; 
    for j = 1:size(trackingKnockout{i,2})
        rxnS(:,trackingKnockout{i,2}) = 0;
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
    biomassGeneration(i) = result.objval;
end
trackingKnockout(:,4) = num2cell(biomassGeneration(:));
trackingKnockout(:,5) = num2cell(cell2mat(trackingKnockout(:,4)) >= 0.5);
viableSingleKnockoutGenes = cell2mat(trackingKnockout(cell2mat(trackingKnockout(:,5)) == 1,3));
viableSingleKnockoutGenesIndices = cell2mat(trackingKnockout(:,4)) == max(cell2mat(trackingKnockout(:,4))); % Array of best biomass flux.
viableKnockoutGenes = num2cell(viableSingleKnockoutGenes(viableSingleKnockoutGenesIndices,:));       % Cell array of best nonlethal genes.
lenViableSingleKnockoutGenes = size(viableSingleKnockoutGenes,1);
clear i x numGene rxnList genomeSize biomassGeneration

%% Multiple Loop
while ~isempty(viableKnockoutGenes)
    [trackingKnockout,viableKnockoutGenes,viableSingleKnockoutGenes] = GeneKnockoutTreeIncrement(viableKnockoutGenes,viableSingleKnockoutGenes,model);
end

% M02 End Program