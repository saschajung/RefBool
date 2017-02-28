function [distributions] = DetermineDistributions( SingleExpressionReference, FittedDists, probdists, bootstrapSamples, precision)
%DETERMINEDISTRIBUTIONS
%   Should be only called from within DetermineThresholdDistributions

firstBound = [];
secondBound = [];
isMixtureModel = strcmpi(probdists.DistributionName,'gaussian mixture distribution');
if isMixtureModel
    g = @(x) infEx(FittedDists,x,precision,probdists.NumComponents);
else
    g = @(x) infEx(FittedDists,x,precision,probdists.NumParameters);
end
[boot_params, ~]=bootstrp(bootstrapSamples,g,SingleExpressionReference','Options',statset('UseParallel',true));
boot_params = boot_params(~any(isnan(boot_params),2),:);
if isMixtureModel
    [out,upper] = MinimizeRectangle(boot_params,probdists.DistributionName,SingleExpressionReference,precision,probdists.NumComponents,probdists.ComponentProportion);
else
    [out,upper] = MinimizeRectangle(boot_params,probdists.DistributionName,SingleExpressionReference,precision);
end
firstBound = [firstBound out];
secondBound = [secondBound upper];

distributions = {firstBound,secondBound};


end

