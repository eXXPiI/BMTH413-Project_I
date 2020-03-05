%% Preamble
% Program: GeneKnockoutTree.m
% Author: Jonathan Myers
% Date: December 8, 2019
% Purpose: Determine gene knockout using tree algorithm.
% Arguments: biomassGeneration, knockoutNumber, and model.
% Loads: None.
% Calls: None.
% Returns: None.

%% Function
function trackingKnockout = GeneKnockoutTree(knockoutNumber,model)
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
    rxnS = model.S; %#ok<PFBNS>
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
clear i x numGene rxnList biomassGeneration
if knockoutNumber == 1
    return
end

for iteration = 2:knockoutNumber
    lenViableKnockoutGenes = size(viableKnockoutGenes,1);
    geneSetNumber = lenViableSingleKnockoutGenes*lenViableKnockoutGenes;
    trackingKnockout = cell(geneSetNumber,5);
    trackingKnockout(:,1) = num2cell(1:geneSetNumber);
    geneSet = cell(geneSetNumber,iteration); % Array of all gene knockouts.
    for i = 0:lenViableKnockoutGenes-1
        for j = 1:lenViableSingleKnockoutGenes
            for k = 1:iteration-1
                geneSet{j+(i*lenViableSingleKnockoutGenes),k} = viableKnockoutGenes{i+1,k};
            end
        end
    end
    for i = 0:lenViableKnockoutGenes-1
        for j = 1:lenViableSingleKnockoutGenes
            geneSet{j+(i*lenViableSingleKnockoutGenes),iteration} = viableSingleKnockoutGenes(j);
        end
    end
    
    for i = 1:geneSetNumber
        trackingKnockout{i,3} = cell2mat(geneSet(i,:));
    end

    for knockoutIncrement = 1:geneSetNumber
        geneVector = cell2mat(geneSet(knockoutIncrement,:));
        rxnList = [];
        for i=1:length(geneVector)
            rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1)); % Reactions that use gene index enzymes.
        end
        if (~isempty(rxnList))
            x = true(lenViableSingleKnockoutGenes);                                      % Build total genes matrix.
            x(geneVector) = false; %#ok<NASGU>                              % Knock out gene indices.
            %trackingKnockout{knockoutIncrement,2} = [];                     % Initialize reactions to remove.
            for i = 1:length(rxnList)
                if (~eval(model.rules{rxnList(i)}))                             % Use determine if x contents 
                    trackingKnockout{knockoutIncrement,2} = union(trackingKnockout{knockoutIncrement,2},rxnList(i));
                end
            end
        end
    end
    
    biomassGeneration = zeros(geneSetNumber,1);
    parfor i = transpose(cell2mat(trackingKnockout(:,1)))
        mod = [];
        param = [];
        rxnS = model.S; %#ok<PFBNS>
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
    viableKnockoutGenes = cell2mat(trackingKnockout(cell2mat(trackingKnockout(:,5)) == 1,3));
    viableKnockoutGenesIndices = cell2mat(trackingKnockout(:,4)) == max(cell2mat(trackingKnockout(:,4))); % Array of best biomass flux.
    if sum(viableKnockoutGenesIndices) == 0
        break
    end
    viableKnockoutGenes = num2cell(viableKnockoutGenes(viableKnockoutGenesIndices,:));       % Cell array of best nonlethal genes.
    clear i x numGene rxnList biomassGeneration
end
end