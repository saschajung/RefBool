function [ Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj ] = CalculatePValues( dists, expression, varargin )
%CALCULATEPVALUES
% Calculate p- and q-values for the expression values in 'expression' with
% respect to the threshold distributions in 'dists'.
% INPUT: dists: Threshold distribution cell array derived from
%               DetermineThresholdDistributions.m
%        expression: Query gene expression values as matrix
% NOTE: the ordering of genes for 'dists' and 'expression' needs to be the
% same!

[r,c] = size(expression);
Pvals_low = zeros(r,c);
Pvals_up = zeros(r,c);
Pvals_inter_low = zeros(r,c);
Pvals_inter_up = zeros(r,c);
parfor i=1:r
    xq = 0:1e-3:max(max(expression+1e-3));
    if isempty(dists{i}{1})
        Pvals_low(i,:) = NaN;
        Pvals_up(i,:) = NaN;
        continue
    end
    [~,cols] = size(dists{i}{1});
    if cols > 1
        [f_low,x_low] = ecdf(dists{i}{1}(:,1),'frequency',dists{i}{1}(:,2));
        [f_up,x_up] = ecdf(dists{i}{2}(:,1),'frequency',dists{i}{2}(:,2));
    else
        [f_low,x_low] = ecdf(dists{i}{1}(:,1));
        [f_up,x_up] = ecdf(dists{i}{2}(:,1));
    end
    x_low(1) = 0;
    x_up(1) = 0;
    if x_low(end) ~= xq(end)
        x_low(end+1) = xq(end);
        f_low(end+1) = 1;
    end
    if x_up(end) ~= xq(end)
        x_up(end+1) = xq(end);
        f_up(end+1) = 1;
    end
    if x_low(1)==x_low(2)
        x_low = [x_low(2);x_low(3:end)];
        f_low = [f_low(2);f_low(3:end)];
    end
    if x_up(1)==x_up(2)
        x_up = [x_up(2);x_up(3:end)];
        f_up = [f_up(2);f_up(3:end)];
    end
    fq_low = interp1(x_low,f_low,xq,'linear','extrap');
    fq_up = interp1(x_up,f_up,xq,'linear','extrap');
    low_pval_i = zeros(1,c);
    up_pval_i = zeros(1,c);
    inter_low_pval_i = zeros(1,c);
    inter_up_pval_i = zeros(1,c);
    for j=1:c
        inter_low_pval_i(1,j) = 1-fq_low(find(xq<=expression(i,j),1,'last'));
        inter_up_pval_i(1,j) = fq_up(find(xq>=expression(i,j),1));
        low_pval_i(1,j) = fq_low(find(xq>=expression(i,j),1));
        up_pval_i(1,j) = 1-fq_up(find(xq<=expression(i,j),1,'last'));
    end
    Pvals_low(i,:) = low_pval_i;
    Pvals_up(i,:) = up_pval_i;
    Pvals_inter_low(i,:) = inter_low_pval_i;
    Pvals_inter_up(i,:) = inter_up_pval_i;
end

%Calculate adjusted P-vals based on Benjamini-Hochberg correction
tmp = reshape(Pvals_low,[c*r,1]);
adjp = mafdr(tmp,'BHFDR',true);
Pvals_low_adj = reshape(adjp(1:c*r,1),[r c]);
tmp = reshape(Pvals_up,[c*r,1]);
adjp = mafdr(tmp,'BHFDR',true);
Pvals_up_adj = reshape(adjp(1:c*r,1),[r c]);
tmp = reshape(Pvals_inter_low,[c*r,1]);
adjp = mafdr(tmp,'BHFDR',true);
Pvals_inter_low_adj = reshape(adjp(1:c*r,1),[r c]);
tmp = reshape(Pvals_inter_up,[c*r,1]);
adjp = mafdr(tmp,'BHFDR',true);
Pvals_inter_up_adj = reshape(adjp(1:c*r,1),[r c]);
end

