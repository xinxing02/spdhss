function M = polynomialkernel(coord, r, d)
%
%   Polynomial Kernel F(x,y) = (x'*y + r) ^ d
%

if (nargin < 3)
    d = 2;
end

if isnumeric(coord)
    M = (coord * coord' + r).^d;
elseif iscell(coord)&&length(coord)==2      
    M = (coord{1} * coord{2}' + r).^d;
end