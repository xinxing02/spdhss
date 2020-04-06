function u = HssMatVecProduct(D, U, B, vec, htree)
%
%   Compute A * vec with A being a HSS matrix (D,U,B)
%
    children = htree.children;
    level    = htree.level;
    nlevel   = length(level);
    root     = htree.root;
    mcluster = htree.mcluster;
    vec = vec(:);
    u = zeros(length(vec), 1);
    y = {[]};
   
    % Upward Sweep
    for i = nlevel : -1 : 1
        for j = 1 : length(level{i})
            node = level{i}(j);
            if (isnan(children(node,1))) %Leaf node
                index = mcluster(node,1):mcluster(node,2);
                u(index) = D{node} * vec(index);
            else %NonLeaf Node 
                c1 = children(node,1);
                c2 = children(node,2);
                if isnan(children(c1,1))
                    y{c1} = U{c1}' * vec(mcluster(c1,1):mcluster(c1,2));
                else
                    y{c1} = U{c1}' * y{c1};
                end
                if isnan(children(c2,1))
                    y{c2} = U{c2}' * vec(mcluster(c2,1):mcluster(c2,2));
                else
                    y{c2} = U{c2}' * y{c2};
                end
                y{node} = [y{c1};y{c2}];
                tmp = y{c1};
                y{c1} = B{node} * y{c2};
                y{c2} = B{node}' * tmp;
            end
        end
    end
    y{root} = zeros(length(y{root}), 1);
    
    %% Downward Sweep
    for i = 1 : nlevel
        for j = 1 : length(level{i})
            node = level{i}(j);           
            if (isnan(children(node,1))) %Leaf node
                index = mcluster(node,1):mcluster(node,2);
                u(index) = u(index) + y{node};
            else %NonLeaf Node
                c1 = children(node,1);
                c2 = children(node,2);
                y{c1} = y{node}(1:length(y{c1})) + y{c1};
                y{c2} = y{node}((length(y{c1})+1):end) + y{c2};
                y{c1} = U{c1} * y{c1};
                y{c2} = U{c2} * y{c2};
            end
        end
    end        
end