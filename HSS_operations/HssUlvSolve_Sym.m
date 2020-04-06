function u = HssUlvSolve_Sym(Q, Lr, Lc, Idx, rhs, htree, flag)
%
%   Use the Generalized Cholesky Factor to multiply a vector. 
%   A = Q*Lr*Lc*Q'
%
%   if flag == 1, u satisfy (Lc*Q')*u = rhs;
%   if flag == -1, u satisfy (Q*Lr)*u = rhs;
%   if flag == 0, u satisfy A * u = rhs;

    
    level = htree.level;
    nlevel = length(level);    
    u = rhs(:);
        
    if (nargin < 7)
        flag = 0;
    end
    
    %   Check Input
    if (length(Q) * length(Lr) * length(Lc) * length(Idx) == 0)
        disp('Invalid Input for ULV factors');
        return
    end
    
    %   three cases
    switch flag
        case 1
            for i = 1:nlevel
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Q{node} * (Lc{node}\u(Idx{node}));
                end
            end
        case -1
            for i = nlevel : -1 : 1
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Lr{node}\(Q{node}'*u(Idx{node}));
                end
            end
        otherwise
            for i = nlevel : -1 : 1
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Lr{node}\(Q{node}'*u(Idx{node}));
                end
            end
            for i = 1 : nlevel
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Q{node} * (Lc{node}\u(Idx{node}));
                end
            end

    end
end