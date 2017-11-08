%%% Example script for booleanization with RefBool
% prerequisites:
%
% 1. get the RefBool code
% 
% $git clone https://github.com/saschajung/RefBool.git
%
% 2. Download the TPM background distribution from:
%
% https://webdav-r3lab.uni.lu/public/cbg/RefBool_ReferenceDistributions/
%
% You can use:
%
% $wget -r -nH --cut-dirs=2 --no-parent https://webdav-r3lab.uni.lu/public/cbg/RefBool_ReferenceDistributions/
%
% 3. Set the path in matlab to have access to both the source and the background files

%% 
% Initialize background and load file to booleanize
background_genes = readtable('GeneOrderReference_CoTFCRF.txt', ...,
    'delimiter','\t', 'ReadVariableNames',false );

gNames = background_genes.(1);

load Thresholds_CoTFcRF.mat;

filename = 'E008_EXPRESSION.txt';
Expr = readtable(filename, 'delimiter', '\t');

%%
% Align expression to the distributions
indices = -ones(size(gNames));
for j = 1:length(indices)
    expind = find(ismember(Expr.gene_id,gNames{j}));
    if ~isempty(expind)
        indices(j)=expind;
    end
end

% Removing non-existing indices
mydists = dists;
toDel = (indices==-1);

mydists(toDel)=[];

indices(toDel)=[];

% Organize vector for expression
expression = Expr{indices,2};
Gene_ID = Expr{indices,1};

%%
% Do the actual booleanization, and get the p-values
[ Pvals_low,Pvals_up,Pvals_inter_low,Pvals_inter_up,Pvals_low_adj,Pvals_up_adj,Pvals_inter_low_adj,Pvals_inter_up_adj ] = CalculatePValues( mydists, expression);

%%
% Generate the output

% booleanizing everything
bool_value_all = Pvals_inter_low_adj < 0.95;

% only lower and upper 10% according to the adjusted p-values
bool_value = nan(size(Pvals_inter_low_adj));

bool_value (Pvals_low_adj < 0.1) = 0;
bool_value (Pvals_up_adj < 0.1) = 1;

output = table(Gene_ID, expression, bool_value_all, bool_value);
background_genes.Properties.VariableNames = {'Gene_ID', 'HGNC_symbol'};
output = innerjoin(background_genes,output);

%%
% Write file
[pathstr,name,ext] = fileparts(which(filename));
writetable(output,[pathstr filesep name '_bool' ext], 'Delimiter','\t');
