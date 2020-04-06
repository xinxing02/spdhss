function u = HssCholProduct(Q, L, Idx, vec, htree, flag)
%
%   Use the HSS Cholesky Factor to multiply a vector. 
%
%   With A = S * S'
%   
%   trans_flag == 1 => 
%       u = S * vec
%   trans_flag == -1 =>
%       u = S' * vec
%   trans_flag ~= -1 =>
%       u = S * S' * vec

    if (nargin < 6)
        flag = 0;
    end
    
    level = htree.level;
    nlevel = length(level);    
    u = vec(:);
       
    %Check Input
    if (length(Q) * length(L) * length(Idx) == 0)
        disp('Invalid Input for ULV factors');
        return
    end
    
    switch flag
        case 1
            for i = 1 : nlevel
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Q{node}*(L{node}*u(Idx{node}));
                end
            end
        case -1
            for i = nlevel : -1 : 1
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = L{node}'*(Q{node}'*u(Idx{node}));
                end
            end
        otherwise
            for i = nlevel : -1 : 1
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = L{node}'*(Q{node}'*u(Idx{node}));
                end
            end
            for i = 1 : nlevel
                for j = 1 : length(level{i})
                    node = level{i}(j);                        
                    u(Idx{node}) = Q{node}*(L{node}*u(Idx{node}));
                end
            end
    end
end