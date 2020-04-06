function [Q,Lr,Lc,Idx] = HssULV_Sym(D, U, B, htree, shift)
%
%   Construct the ULV decomposition of a general symmetric HSS matrix. 
%
%   Q:  Orthogonal matrix w.r.t. each node's basis
%   Lr, Lc:  LU factor w.r.t. each node's diagonal block
%   Idx: Record the index of the submatrix which Q and L performs on for
%   each node in global sense. 

    if (nargin < 7)
        shift = 0;
    end
    parent   = htree.parent;
    children = htree.children;
    level    = htree.level;
    nlevel   = length(level);
    root     = htree.root;
    mcluster = htree.mcluster;
    
    Q = cell(length(parent));         % Orthogonal matrix
    Lr = cell(length(parent));        % LU factor
    Lc = cell(length(parent));        % LU factor
    D_mid = cell(length(parent));     % Diagonal block 
    U_mid = cell(length(parent));     % Basis changed
    
    for i = nlevel : -1 : 2
        for j = 1 : length(level{i})
            node = level{i}(j);                        
            if (isnan(children(node,1))) %Leaf node
                %Size information of U_node and Diag_node
                numrank = size(U{node},2);  
                nrow = size(D{node},1);
  
                %QL decomposition of U_node
                [Q{node}, tmpU] = qr(U{node});
                U_mid{node} = tmpU(1:numrank,:);
                Q{node}(:, [1:numrank, (numrank+1):end]) = Q{node}(:,  [(numrank+1):end, 1:numrank]);
                
                %LU Decomposition over left-top corner of Diag{node}
                %and eliminate the off-diagonal block of it. 
                tmpD = Q{node}' * (D{node}+shift*eye(size(D{node}))) * Q{node};
                [tmpL, tmpU] = lu(tmpD(1:(end-numrank), 1:(end-numrank)));
                if (min(abs(diag(tmpU))) < eps('double'))
                    disp('The target matrix is near-singular');
                    Q = {};
                    Lr = {};
                    Lc = {};
                    Idx = {};
                    return 
                end
                LD12 = tmpL \ tmpD(1:(end-numrank), (end-numrank+1):end);
                DU21 = tmpD((end-numrank+1):end, 1:(end-numrank)) / tmpU;  
                
                %   factors and remaining schur complement
                Lr{node} = [tmpL, zeros(nrow-numrank, numrank); DU21, eye(numrank)];
                Lc{node} = [tmpU, LD12; zeros(numrank, nrow-numrank), eye(numrank)];
                D_mid{node} = tmpD((end-numrank+1):end,(end-numrank+1):end) - DU21 * LD12;
                Idx{node} = mcluster(node,1):mcluster(node,2);
            else %NonLeaf Node
                %Children and size information
                c1 = children(node,1);
                c2 = children(node,2);
                numrank = size(U{node},2);
                nrow = size(D_mid{c1},1) + size(D_mid{c2},1);
                
                %Build the compressed diagonal block
                Offdiag = U_mid{c1}* B{node} *U_mid{c2}';
                tmpD = [D_mid{c1}, Offdiag; Offdiag', D_mid{c2}];
                
                %Build the compressed basis  
                tmpU = [U_mid{c1}*U{node}(1:size(U_mid{c1},2),:); 
                                U_mid{c2}*U{node}((size(U_mid{c1},2)+1):end,:)];
                            
                %QL decomposition of U_node
                [Q{node}, tmpU] = qr(tmpU);
                U_mid{node} = tmpU(1:numrank,:);
                Q{node}(:, [1:numrank, (numrank+1):end]) = Q{node}(:,  [(numrank+1):end, 1:numrank]);
                
                %LU Decomposition over left-top corner of Diag{node}
                %and eliminate the off-diagonal block of it. 
                tmpD = Q{node}' * tmpD * Q{node};
                [tmpL, tmpU] = lu(tmpD(1:(end-numrank), 1:(end-numrank)));
                if (min(abs(diag(tmpU))) < eps('double'))
                    disp('The target matrix is near-singular');
                    Q = {};
                    Lr = {};
                    Lc = {};
                    Idx = {};
                    return 
                end
                LD12 = tmpL \ tmpD(1:(end-numrank), (end-numrank+1):end);
                DU21 = tmpD((end-numrank+1):end, 1:(end-numrank)) / tmpU;
                
                %   factors and remaining schur complement
                Lr{node} = [tmpL, zeros(nrow-numrank, numrank); DU21, eye(numrank)];
                Lc{node} = [tmpU, LD12; zeros(numrank, nrow-numrank), eye(numrank)];
                D_mid{node} = tmpD((end-numrank+1):end,(end-numrank+1):end) - DU21 * LD12;
                Idx{node} = [Idx{c1}((end-size(D_mid{c1},1)+1):end), ...
                                Idx{c2}((end-size(D_mid{c2},1)+1):end)];                
            end
        end
    end

    %Last Step
    c1 = children(root, 1);
    c2 = children(root, 2);
    Offdiag = U_mid{c1}*B{root}*U_mid{c2}';
    D_mid{root} = [D_mid{c1}, Offdiag; Offdiag', D_mid{c2}];
    
    [Lr{root},Lc{root}] = lu(D_mid{root});
    if (min(abs(diag(Lc{root}))) < eps('double'))
        disp('The target matrix is near-singular');
        Q = {};
        Lr = {};
        Lc = {};
        Idx = {};
        return 
    end
    Q{root} = eye(size(D_mid{root}));
    Idx{root} = [Idx{c1}((end-size(D_mid{c1},1)+1):end), ...
                    Idx{c2}((end-size(D_mid{c2},1)+1):end)];
end