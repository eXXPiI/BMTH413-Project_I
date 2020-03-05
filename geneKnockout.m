fprintf('loading model file\n')
load Ec_iJO1366.mat;

% Construct the flux lists to remove
fprintf('Calculating fluxes list to remove\n')
removeList = cell(length(model.genes),1);
for numGene = 1:length(model.genes)
   geneVector = [numGene];
   rxnList = [];
   for i=1:length(geneVector)
      rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
   end
   rxnList = sort(rxnList);
   if (~isempty(rxnList))
      x = true(size(model.genes));
      x(geneVector) = false;
      removeList{numGene} = [];
      for i = 1:length(rxnList)
         if (~eval(model.rules{rxnList(i)})) % Looks at x vector.
            removeList{numGene} = union(removeList{numGene},rxnList(i));
         end
      end
   end
end
