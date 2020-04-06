function M = laplacian(coord, e)
%
%   Exponential Kernel matrix, K(x,y) = exp(-e * ||x-y||)
%
if isnumeric(coord)
    N = size(coord, 1);
    dist = pdist(coord);
    M = squareform(exp(- e * dist)) + eye(N);
elseif iscell(coord)&&length(coord)==2      
    M = pdist2(coord{1}, coord{2});
    M = exp(- e * M);
end
end