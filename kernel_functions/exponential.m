function M = exponential(coord, e)
%
%   Exponential Kernel matrix, K(x,y) = exp(-eps * ||x-y||^2)
%
if isnumeric(coord)
    N = size(coord, 1);
    dist = pdist(coord);
    M = squareform(exp(- e * dist.^2)) + eye(N);
elseif iscell(coord)&&length(coord)==2      
    M = pdist2(coord{1}, coord{2});
    M = exp(- e * M.^2);
end
end