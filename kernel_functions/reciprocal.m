function M = reciprocal(coord, alpha)
%
%   Reciprocal Kernel matrix, K(x,y) = 1 / ||x-y||^alpha
%
if isnumeric(coord)
    N = size(coord, 1);
    dist = pdist(coord);
    M = squareform(1./dist) + eye(N);
elseif iscell(coord)&&length(coord)==2      
    M = pdist2(coord{1}, coord{2});
    M = 1./(M.^alpha);
end
end