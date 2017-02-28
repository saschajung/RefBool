function [ DiscretizedMatrix ] = Booleanize_Bypval( qvals_low, qvals_up, qvals_inter_low, qvals_inter_up, threshold, threshold_inter )
%Booleanize_Bypval
% Given pvals and qvals (obtained from CalculatePValues.m) and a threshold
% for the significance, a matrix with the discretized values is returned.

[r,c] = size(qvals_low);
DiscretizedMatrix = zeros(r,c);
k = 0;
for i=1:r
   for j=1:c
      if qvals_low(i,j) <= threshold && qvals_up(i,j) <= threshold %Should not happen
          DiscretizedMatrix(i,j) = NaN;
          k = k+1;
      else if qvals_low(i,j) <= threshold
          DiscretizedMatrix(i,j) = 0;
          else if qvals_up(i,j) <= threshold
                  DiscretizedMatrix(i,j) = 1;
              else if qvals_inter_low(i,j) >= 1-threshold_inter && qvals_inter_up(i,j) <= threshold_inter
                  DiscretizedMatrix(i,j) = 0.5;
                  else
                      DiscretizedMatrix(i,j) = NaN;
                  end
              end
          end
      end
        
   end
end
display(sprintf('%d values had significant q-values for active and inactive expression!',k))
end

