function A = Hss2Mat(D, U, B, htree)
%
%  Reconstruct the dense matrix from the HSS matrix
%
    level    = htree.level;
    children = htree.children;
    mcluster = htree.mcluster;
    root     = htree.root;
    N = (mcluster(end,2)-mcluster(end,1) + 1);
    A = zeros(N);
    
    nlevel = length(level);
    tmpU = {[]};
    for i = nlevel: -1 : 1
        for j = 1 : length(level{i})
            node = level{i}(j);
            index = mcluster(node,1):mcluster(node,2);
            
            if (isnan(children(node,1)))
                %Diagonal Block
                A(index, index) = D{node};
                tmpU{node} = U{node};
            else %NonLeaf Node
                c1 = children(node,1);
                c2 = children(node,2);
                if (node ~= root) %not root
                    tmpU{node} = [tmpU{c1}*U{node}(1:size(tmpU{c1},2),:);
                            tmpU{c2}*U{node}((size(tmpU{c1},2)+1):end,:)];
                end
                A(mcluster(c1,1):mcluster(c1,2), mcluster(c2,1):mcluster(c2,2)) = ...
                    tmpU{c1}*B{node}*tmpU{c2}';
                A(mcluster(c2,1):mcluster(c2,2), mcluster(c1,1):mcluster(c1,2)) = ...
                    A(mcluster(c1,1):mcluster(c1,2), mcluster(c2,1):mcluster(c2,2))';
            end
        end
    end
    A = (A+A')/2;
end