function store = HssStorage(U, htree, mcluster)
%
%   Storage Cost for the HSS matrix
%
    children = htree.children;
    ncluster = size(mcluster,1);
    store = 0;
    for i = 1 : (ncluster - 1)        
        if isnan(children(i,1))
            %Diagonal
            store = store + (mcluster(i,2)-mcluster(i,1)+1)^2;
            %U_i
            store = store + numel(U{i});
        else
            c1 = children(i,1);
            c2 = children(i,2);
            store = store + size(U{c1},2) * size(U{c2}, 2);
            store = store + numel(U{i});
        end
    end
end
