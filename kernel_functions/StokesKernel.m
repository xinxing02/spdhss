function D = StokesKernel(coord, a, eta)
% The RPY tensor estimation of the covariance matrix D
% Var:
%   coord, the coordinates of the particles' present location
%          N*3 matrix
%   a, the radius of the particle
%   D0, k_bT/6pna
% Output:
%   D, the dense 3N*3N covariance matrix

N = size(coord,1);
D0 = 1/(6*pi*a*eta);
D = D0*eye(3*N);

for i = 2 : N
    for j = 1 : (i-1)
        %fill in the D_ij, a 3*3 matrix
        r = coord(i,:) - coord(j,:);
        rind = ((i-1)*3+1):(3*i);
        cind = ((j-1)*3+1):(3*j);
        l = norm(r);
%         if (l >= 2*a)
%             D(rind, cind) = D0*3/4*a/l*(eye(3) + 1/l/l*(r'*r));
%         else
%             D(rind, cind) = D0*( (1-9/32*l/a)*eye(3) + 3/32/l/a*(r'*r) );
%         end

        D(rind, cind) = D0*3/4*a/l*(eye(3) + 1/l/l*(r'*r));
        D(cind, rind) = D(rind, cind);
    end
end

end
