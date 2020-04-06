function [Q,L,Idx] = HssCholesky(D, U, B, htree, shift)
%
%   Construct the ULV decomposition of the s.p.d. HSS matrix. 
%
%   Q:  Orthogonal matrix w.r.t. each node's basis
%   L:  Cholesky factor w.r.t. each node's diagonal block
%   Idx: Record the index of the submatrix which Q and L performs on for
%   each node in global sense. 

    if (nargin < 7)
        shift = 0;
    end
    parent   = htree.parent;
    children = htree.children;
    mcluster = htree.mcluster;
    level    = htree.level;
    nlevel   = length(level);
    root     = htree.root;
    
    Q = cell(length(parent));         % Orthogonal matrix
    L = cell(length(parent));         % Cholesky Factor
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
                
                %Cholesky Decomposition over left-top corner of Diag{node}
                %and eliminate the off-diagonal block of it. 
                tmpD = Q{node}' * (D{node}+shift*eye(size(D{node}))) * Q{node};
                [tmpL,p] = chol(tmpD(1:(end-numrank), 1:(end-numrank)), 'lower');
                if (p ~= 0)
                    disp('Matrix Not Positive Definite');
                    Q = {};
                    L = {};
                    Idx = {};
                    return 
                end
                S = tmpL\tmpD(1:(end-numrank), (end-numrank+1):end);
                
                %Cholesky factor and remainder (Schur complement)
                L{node} = [tmpL, zeros(nrow-numrank, numrank); S', eye(numrank)];
                D_mid{node} = tmpD((end-numrank+1):end,(end-numrank+1):end) - S'*S;
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
                
                %Cholesky Decomposition over left-top corner of Diag{node}
                %and eliminate the off-diagonal block of it. 
                tmpD = Q{node}' * tmpD * Q{node};
                [tmpL,p] = chol(tmpD(1:(end-numrank), 1:(end-numrank)), 'lower');
                if (p ~= 0)
                    disp('Matrix Not Positive Definite');
                    Q = {};
                    L = {};
                    Idx = {};
                    return 
                end
                
                S = tmpL\tmpD(1:(end-numrank), (end-numrank+1):end);                
                %Cholesky factor and remainder (Schur complement)
                L{node} = [tmpL, zeros(nrow-numrank, numrank); S', eye(numrank)];
                D_mid{node} = tmpD((end-numrank+1):end,(end-numrank+1):end) - S'*S;
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
    
    [L{root},p] = chol(D_mid{root}, 'lower');
    if (p ~= 0)
        disp('Matrix Not Positive Definite');
        Q = {};
        L = {};
        Idx = {};
        return 
    end
    Q{root} = eye(size(D_mid{root}));
    Idx{root} = [Idx{c1}((end-size(D_mid{c1},1)+1):end), ...
                    Idx{c2}((end-size(D_mid{c2},1)+1):end)];
end