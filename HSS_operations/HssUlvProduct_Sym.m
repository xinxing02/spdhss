function u = HssUlvProduct_Sym(Q, Lr, Lc, Idx, vec, htree, flag)
%
%   Use the Generalized Cholesky Factor to multiply a vector. 
%
%   With A = L * U
%   
%   trans_flag == 1 => 
%       u = L * vec
%   trans_flag == -1 =>
%       u = U * vec

    if (nargin < 7)
        flag = 0;
    end
    
    level = htree.level;
    nlevel = length(level);    
    u = vec(:);
       
    %Check Input
    if (length(Q) * length(Lr) * length(Lc) * length(Idx) == 0)
        disp('Invalid Input for ULV factors');
        return
    end
    
    switch flag
        case 1
            for i = 1 : nlevel
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Q{node}*(Lr{node}*u(Idx{node}));
                end
            end
        case -1
            for i = nlevel : -1 : 1
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Lc{node}*(Q{node}'*u(Idx{node}));
                end
            end
        otherwise
            for i = nlevel : -1 : 1
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Lc{node}*(Q{node}'*u(Idx{node}));
                end
            end
            for i = 1 : nlevel
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Q{node}*(Lr{node}*u(Idx{node}));
                end
            end
    end
end