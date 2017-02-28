function np = infEx( type, x, precision, varargin )
%INFEX 
%   Estimates the parameters for a bootstrap sample given the distribution
%   'type'. In case of generalized pareto distributions, the fit could
%   abort and NaN is returned instead indicating that this value should not
%   be taken into account.
%   Should be only called from within DetermineDistributions
warning('off','all')
if strcmpi(type,'Generalized Pareto')
    y = min(x(:))-0.1:precision:min(x(:));
    try
        pd = fitdist(x,'Generalized Pareto','theta',y(find(y<min(x(:)),1,'last')),'Options',statset('MaxIter',1000,'MaxFunEvals',5000));
    catch ME
        np = NaN*zeros(3,1);
        return
    end
else if strcmpi(type,'Birnbaum-Saunders')
        try
        pd = fitdist(x,'birnbaumsaunders','Options',statset('MaxIter',1000,'MaxFunEvals',5000));
        catch err
            np = NaN*zeros(1,varargin{1});
            return;
        end
    else if strcmpi(type,'Inverse Gaussian')
            try
            pd = fitdist(x,'inversegaussian','Options',statset('MaxIter',1000,'MaxFunEvals',5000));
            catch err
                np = NaN*zeros(1,varargin{1});
                return;
            end
        else if strcmpi(type,'Log-Logistic')
                try
                pd = fitdist(x,'loglogistic','Options',statset('MaxIter',1000,'MaxFunEvals',5000));
                catch err
                    np = NaN*zeros(1,varargin{1});
                    return;
                end
            else if strcmpi(type,'t Location-Scale')
                    try
                    pd = fitdist(x,'tlocationscale','Options',statset('MaxIter',1000,'MaxFunEvals',5000));
                    catch err
                        np = NaN*zeros(1,varargin{1});
                        return;
                    end
                else if strcmpi(type,'gaussian mixture distribution')
                        try
                            pd = fitgmdist(x,varargin{1});
                        catch err
                            np = NaN*zeros(1,varargin{1}*3);
                            return
                        end
                    else
                        try
                        pd = fitdist(x,type,'Options',statset('MaxIter',1000,'MaxFunEvals',5000));
                        catch err
                           np = NaN*zeros(1,varargin{1});
                           return
                        end
                    end
                end
            end
        end
    end
end

if strcmpi(type,'gaussian mixture distribution')
    np = [pd.mu',reshape(pd.Sigma,[1,varargin{1}]),pd.ComponentProportion];
else
    np = pd.ParameterValues;
end

end

