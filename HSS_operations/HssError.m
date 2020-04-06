function err = HssError(A, htree, D, U, B)
%
%   This function will calculate the approximation at the off-diagonal
%   blocks. 
%   For leaf node i, 
%       err(i,1) = norm(A(clst(i),clst(i) 'fro');
%       err(i,2) = 0;
%   For non-leaf node i with son c1, c2
%       err(i,1) = norm(A(clst(c1), clst(c2)), 'fro')
%       err(i,2) = norm(A(clst(c1), clst(c2)) - U(c1) * B(i) * U(c2)',
%       'fro');
%   
    level    = htree.level;
    children = htree.children;
    root     = htree.root;
    nlevel   = length(level);
    mcluster = htree.mcluster;
    err = zeros(size(mcluster,1), 2);
    tmpU = {[]};
    
    %   non-root level
    for i = nlevel : -1 : 2
        for j = 1 : length(level{i})
            node = level{i}(j);               
            if (isnan(children(node, 1)))
                err(node,1) = norm(D{node}, 'fro');
                err(node,2) = 0;
                tmpU{node} = U{node};
            else
                c1 = children(node,1);
                c2 = children(node,2);
                err(node,1) = sqrt(2)*norm(A(mcluster(c1,1):mcluster(c1,2), mcluster(c2,1):mcluster(c2,2)), 'fro');
                err(node,2) = sqrt(2)*norm(A(mcluster(c1,1):mcluster(c1,2), mcluster(c2,1):mcluster(c2,2)) - ...
                    tmpU{c1} * B{node} * tmpU{c2}', 'fro');
                tmpU{node} = [tmpU{c1}*U{node}(1:size(tmpU{c1},2),:);
                             tmpU{c2}*U{node}((size(tmpU{c1},2)+1):end,:)];
            end
        end
    end
    
    %   root level
    node = root; 
    c1 = children(node,1);
    c2 = children(node,2);
    err(node,1) = sqrt(2)*norm(A(mcluster(c1,1):mcluster(c1,2), mcluster(c2,1):mcluster(c2,2)), 'fro');
    err(node,2) = sqrt(2)*norm(A(mcluster(c1,1):mcluster(c1,2), mcluster(c2,1):mcluster(c2,2)) - ...
                    tmpU{c1} * B{node} * tmpU{c2}', 'fro');
end
