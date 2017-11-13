function [out,upper] = MinimizeRectangle( bootstrapParams, type, data, precision, varargin)
%MINIMIZERECTANGLE
% Maximizes the rectangles C_1 and C_2 as described in the manuscript.
% Should only be called from within DetermineDistributions

maxVals = 10000;

PDs = cell(numel(bootstrapParams(:,1)),1);

if any(strcmpi(type,{'beta'}))
    for i=1:numel(bootstrapParams(:,1))
       PDs{i} =  makedist(type,'a',bootstrapParams(i,1),'b',bootstrapParams(i,2));
    end
else if any(strcmpi(type,{'birnbaumsaunders','Birnbaum-Saunders'}))
         for i=1:numel(bootstrapParams(:,1))
            PDs{i} =  makedist('birnbaumsaunders','beta',bootstrapParams(i,1),'gamma',bootstrapParams(i,2));
         end
    else if strcmpi(type,'burr')
             for i=1:numel(bootstrapParams(:,1))
                PDs{i} =  makedist(type,'alpha',bootstrapParams(i,1),'c',bootstrapParams(i,2),'k',bootstrapParams(i,3));
             end            
        else if strcmpi(type,'exponential')
                 for i=1:numel(bootstrapParams(:,1))
                    PDs{i} =  makedist(type,'mu',bootstrapParams(i,1));
                 end                
            else if strcmpi(type,'extreme value')
                     for i=1:numel(bootstrapParams(:,1))
                        PDs{i} =  makedist(type,'mu',bootstrapParams(i,1),'sigma',bootstrapParams(i,2));
                     end                    
                else if strcmpi(type,'gamma')
                         for i=1:numel(bootstrapParams(:,1))
                            PDs{i} =  makedist(type,'a',bootstrapParams(i,1),'b',bootstrapParams(i,2));
                         end                        
                    else if strcmpi(type,'generalized extreme value')
                             for i=1:numel(bootstrapParams(:,1))
                                PDs{i} =  makedist(type,'k',bootstrapParams(i,1),'sigma',bootstrapParams(i,2),'mu',bootstrapParams(i,3));
                             end
                        else if strcmpi(type,'generalized pareto')
                                 for i=1:numel(bootstrapParams(:,1))
                                    PDs{i} =  makedist(type,'k',bootstrapParams(i,1),'sigma',bootstrapParams(i,2),'theta',bootstrapParams(i,3));
                                 end                                
                            else if any(strcmpi(type,{'inversegaussian','Inverse Gaussian'}))
                                     for i=1:numel(bootstrapParams(:,1))
                                        PDs{i} =  makedist('inversegaussian','mu',bootstrapParams(i,1),'lambda',bootstrapParams(i,2));
                                     end                                    
                                else if strcmpi(type,'logistic')
                                         for i=1:numel(bootstrapParams(:,1))
                                            PDs{i} =  makedist(type,'mu',bootstrapParams(i,1),'sigma',bootstrapParams(i,2));
                                         end                                        
                                    else if any(strcmpi(type,{'loglogistic','Log-Logistic'}))
                                             for i=1:numel(bootstrapParams(:,1))
                                                PDs{i} =  makedist('loglogistic','mu',bootstrapParams(i,1),'sigma',bootstrapParams(i,2));
                                             end                                            
                                        else if strcmpi(type,'lognormal')
                                                 for i=1:numel(bootstrapParams(:,1))
                                                    PDs{i} =  makedist(type,'mu',bootstrapParams(i,1),'sigma',bootstrapParams(i,2));
                                                 end                                                
                                            else if strcmpi(type,'nakagami')
                                                     for i=1:numel(bootstrapParams(:,1))
                                                        PDs{i} =  makedist(type,'mu',bootstrapParams(i,1),'omega',bootstrapParams(i,2));
                                                     end                                                    
                                                else if strcmpi(type,'normal')
                                                         for i=1:numel(bootstrapParams(:,1))
                                                            PDs{i} =  makedist(type,'mu',bootstrapParams(i,1),'sigma',bootstrapParams(i,2));
                                                         end                                                        
                                                    else if strcmpi(type,'rayleigh')
                                                             for i=1:numel(bootstrapParams(:,1))
                                                                PDs{i} =  makedist(type,'b',bootstrapParams(i,1));
                                                             end                                                            
                                                        else if strcmpi(type,'rician')
                                                                 for i=1:numel(bootstrapParams(:,1))
                                                                    PDs{i} =  makedist(type,'s',bootstrapParams(i,1),'sigma',bootstrapParams(i,2));
                                                                 end                                                                
                                                            else if any(strcmpi(type,{'tlocationscale','t Location-Scale'}))
                                                                     for i=1:numel(bootstrapParams(:,1))
                                                                        PDs{i} =  makedist('tlocationscale','mu',bootstrapParams(i,1),'sigma',bootstrapParams(i,2),'nu',bootstrapParams(i,3));
                                                                     end                                                                    
                                                                else if strcmpi(type,'weibull')
                                                                         for i=1:numel(bootstrapParams(:,1))
                                                                            PDs{i} =  makedist(type,'a',bootstrapParams(i,1),'b',bootstrapParams(i,2));
                                                                         end
                                                                    else if strcmpi(type,'gaussian mixture distribution')
                                                                            for i=1:numel(bootstrapParams(:,1))
                                                                                PDs{i} = gmdistribution(bootstrapParams(i,1:varargin{1})',reshape(bootstrapParams(i,varargin{1}+1:2*varargin{1}),[1,1,varargin{1}]),bootstrapParams(i,2*varargin{1}+1:end));
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

matrix_size = ceil((max(data)/precision)/maxVals);
x_last = 0;
best_low = zeros(numel(bootstrapParams(:,1)),2);
best_up = zeros(numel(bootstrapParams(:,1)),2);

if strcmpi(type,'gaussian mixture distribution')
    %Sort bootstrap samples by mean
    for i=1:numel(bootstrapParams(:,1))
       [~,I] = sort(bootstrapParams(i,1:varargin{1}));
       I = [I,I+varargin{1},I+2*varargin{1}];
       bootstrapParams(i,:) = bootstrapParams(i,I);
    end
    out = [];
    upper = [];
    for comp=1:varargin{1}
        x_last = 0;
        best_low = zeros(numel(bootstrapParams(:,1)),2);
        best_up = zeros(numel(bootstrapParams(:,1)),2);
        for j=1:matrix_size
            %Determine area to explore, i.e. from x_last to x_cur such that only
            %maxVals evaluations are needed given the precision.
            x_cur = x_last + ((min(j*maxVals,max(data)/precision) - ((j-1)*maxVals))*precision);

            x = repmat(x_last:precision:x_cur,numel(bootstrapParams(:,1)),1);
            y = zeros(numel(bootstrapParams(:,1)),numel(x(1,:)));

            for i=1:numel(bootstrapParams(:,1))
                y(i,:) = cdf('normal',x(i,:),bootstrapParams(i,comp),bootstrapParams(i,varargin{1}+comp));
            end

            [m, c] = max(x-x.*y,[],2);
            c = c + ((j-1)*maxVals);
            [m1, c1] = max(y.*(max(data)-x),[],2);
            c1 = c1 + ((j-1)*maxVals);

            if j==1
                best_low = [m c];
                best_up = [m1 c1];
            else
                [m_best, c_best] = max([best_low(:,1) m],[],2);
                [m1_best, c1_best] = max([best_up(:,1) m1],[],2);
                tmp = [best_low(:,2) c];
                best_low = [m_best diag(tmp(:,c_best))];
                tmp = [best_up(:,2) c1];
                best_up = [m1_best diag(tmp(:,c1_best))];
            end

            x_last = x_cur;
        end
        
        frequency = zeros(numel(bootstrapParams(:,1)),1);
        for i = 1:numel(bootstrapParams(:,1))
            frequency(i) = round(bootstrapParams(i,2*varargin{1}+comp)*1000);
        end
        out = [out;[(best_low(:,2)-1)*precision,frequency]];
        upper = [upper;[(best_up(:,2)-1)*precision,frequency]];
    end
else
    for j=1:matrix_size
        %Determine area to explore, i.e. from x_last to x_cur such that only
        %maxVals evaluations are needed given the precision.
        x_cur = x_last + ((min(j*maxVals,max(data)/precision) - ((j-1)*maxVals))*precision);

        x = repmat(x_last:precision:x_cur,numel(bootstrapParams(:,1)),1);
        y = zeros(numel(bootstrapParams(:,1)),numel(x(1,:)));

        for i=1:numel(bootstrapParams(:,1))
            y(i,:) = cdf(PDs{i},x(i,:));
        end

        [m, c] = max(x-x.*y,[],2);
        c = c + ((j-1)*maxVals);
        [m1, c1] = max(y.*(max(data)-x),[],2);
        c1 = c1 + ((j-1)*maxVals);

        if j==1
            best_low = [m c];
            best_up = [m1 c1];
        else
            [m_best, c_best] = max([best_low(:,1) m],[],2);
            [m1_best, c1_best] = max([best_up(:,1) m1],[],2);
            tmp = [best_low(:,2) c];
            best_low = [m_best diag(tmp(:,c_best))];
            tmp = [best_up(:,2) c1];
            best_up = [m1_best diag(tmp(:,c1_best))];
        end

        x_last = x_cur;
    end

    out = (best_low(:,2)-1)*precision;
    upper = (best_up(:,2)-1)*precision;
end
end

