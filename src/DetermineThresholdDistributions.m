function alldists = DetermineThresholdDistributions(geneExpressionReference, bootstrapSamples, precision, sortby, takeMaxThresholds)
%DetermineThresholdDistributions
% INPUT:
% geneExpressionReference: m by n matrix containing m genes with n samples
%                          each. This is the reference distribution for
%                          each gene
% bootstrapSamples: how many bootstrap samples should be drawn to
%                   approximate the parameter distributions. (natural
%                   number)
% precision: the precision of the thresholds. A value of 10^(-x) indicates
%            a precision up to the x-th decimal place.
% takeMaxThresholds: boolean value indicating whether to choose the maximum
% (TRUE) or minimum (FALSE) observed gene expression value if no extreme
% point is detected for the lower threshold. Minimum value is recommended!
% OUTPUT:
% cell array that contains for each gene two distributions, one for the
% lower thresholds and one for the upper thresholds. Used later on to
% determine p-values and q-values of query genes.

[num_genes,~] = size(geneExpressionReference);

FittedDists = cell(num_genes,1);
ProbDists = cell(num_genes,1);

%allfitdist folder is necessary in the path
addpath('allfitdist')

%Get distribution type for each reference gene
parfor g=1:num_genes
   [~,PD] = allfitdist(geneExpressionReference(g,:),sortby);
   FittedDists{g} = PD{1}.DistributionName;
   ProbDists{g} = PD{1};
end


alldists = cell(num_genes,1);
parfor g=1:num_genes
    [dists] = DetermineDistributions(geneExpressionReference(g,:),FittedDists{g},ProbDists{g},bootstrapSamples,precision);
    if ~takeMaxThresholds
       for i=1:numel(dists{1}(:,1))
          if dists{1}(i,1) > dists{2}(i,1)
             dists{1}(i,1) = min(geneExpressionReference(g,:)); 
          end
       end
    end
    alldists{g} = dists;
end

end

