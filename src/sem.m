function y  = sem(x,dim)
% Hannah Payne 9/30/13

if nargin==1    
  y = nanstd(x)./sqrt(sum(~isnan(x)));
else
  y = nanstd(x,0,dim)./sqrt(sum(~isnan(x),dim));
end


