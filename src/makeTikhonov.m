function T = makeTikhonov(n, type, mask)
% n can either be the matrix length, or for 0th order, can be a vector of
% weights
% mask: mask out transitions

if ~exist('type','var')
    type = 0;
end

if length(n)==1
    weights = ones(n,1);
else
    weights = n;
    n = length(weights);
end

% Make tikhonov matrix
% Sum of weights
if type == 0  
        T = diag(weights);  
        
    % First derivative
elseif type == 1
    T = diag(ones(n,1)) - diag(ones(n-1, 1),1);
    T(end,:) = 0;
    
    if exist('mask','var')
        % CHECK THIS: should just be mask?      T(mask,:) = 0;
        T(mask,:) = 0; % Delete the nth points due to overlap btwn conds
        %        temp = find(mask)-1; temp(temp<1) = [];
        %        T(temp,:) = 0;
    end
    T = bsxfun(@times, T, weights(:));
    
    % Second derivative
elseif type == 2
    T =  diag(ones(n,1)) - .5*diag(ones(n-1, 1),1) - .5*diag(ones(n-1, 1),-1);
    T([1 end],:) = 0;
    
    if exist('mask','var')
        T(mask,:) = 0; % Delete the nth points due to overlap btwn conds
        temp = find(mask)-1; temp(temp<1) = [];
        T(temp,:) = 0; % Delete the n-1th points due to overlap btwn conds
    end
    
    T = bsxfun(@times, T, weights(:));

end
