function Q = randSVD(A, r, p)
%
%   Random SVD of the matrix A, 
%
%   Presently, only give the orthogonal columns out as the projector
%

if (r > min(size(A)))
    r = min(size(A));
end

if (nargin < 3)
    p = 5;
end

randvec = rand(size(A,2), r+p);
ax = A * randvec;
[Q0, ~] = qr(ax, 0);
[Q1, ~] = qr(A' * Q0, 0);
[Q, ~, ~] = qr(A * Q1);
Q = Q(:,1:r);
end