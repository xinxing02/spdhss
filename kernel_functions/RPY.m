function D = RPY(coord, a, eta)
% The RPY tensor estimation of the covariance matrix D
% Var:
%   coord, the coordinates of the particles' present location
%          N*3 matrix
%   a, the radius of the particle
%   D0, k_bT/6pna
% Output:
%   D, the dense 3N*3N covariance matrix

if isnumeric(coord)
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
            if (l >= 2*a)
                D(rind, cind) = D0*3/4*a/l*( (1+2/3*a*a/l/l) * eye(3) + 1/l/l*(1-2*a*a/l/l)*(r'*r));
            else
                D(rind, cind) = D0*( (1-9/32*l/a)*eye(3) + 3/32/l/a*(r'*r) );
            end
            D(cind, rind) = D(rind, cind);
        end
    end
else if iscell(coord)&&length(coord)==2        
        N = size(coord{1},1);
        M = size(coord{2},1);
        D0 = 1/(6*pi*a*eta);
        D = zeros(3*N,3*M);
        for i = 1 : N
            for j = 1 : M
                %fill in the D_ij, a 3*3 matrix
                r = coord{1}(i,:) - coord{2}(j,:);
                rind = ((i-1)*3+1):(3*i);
                cind = ((j-1)*3+1):(3*j);
                l = norm(r);
                if (l >= 2*a)
                    D(rind, cind) = D0*3/4*a/l*( (1+2/3*a*a/l/l) * eye(3) + 1/l/l*(1-2*a*a/l/l)*(r'*r));
                else
                    D(rind, cind) = D0*( (1-9/32*l/a)*eye(3) + 3/32/l/a*(r'*r) );
                end
            end
        end
    else
    D = [];
    end
end

end
