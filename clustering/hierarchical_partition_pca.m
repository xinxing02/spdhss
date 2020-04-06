function [parent, children, cluster, permute, n] = hierarchical_partition_pca(coord, minSize)

%   Hierarchically bipartition of a given point set.
%
%   Input: 
%       coord:    N*dim matrix
%       minSize:  minimum size for each leaf box
%   
%   Output:
%       n:        number of nodes in the tree
%       parent:   n*1 vector, the ith entry is the index of i's parent node
%       children: n*2 matrix, the ith row contains the indices of i's two children.
%       cluster:  n*2 matrix, the ith row, [start-point, end-point],
%                   contains the range of the indices of points lying
%                   in box i. 
%       permute:  n*1 vector, the permutation of the coordinates needed that 
%                      makes each cluster of points have consecutive
%                      indices.

    N = size(coord,1);
    
    %   Small cluster
    if (N <= minSize)
        parent = [0];
        children = [nan,nan];
        cluster = [1,N];
        permute = (1:N)';
        n = 1; 
        return ;
    end

    %   Bisect the particles by PCA
    [index1,index2] = bisection_pca(coord);
    
    %   Special case: coordinates coincide.
    if (isempty(index1)||isempty(index2))        
        parent = [0];
        children = [nan,nan];
        cluster = [1,N];
        permute = (1:N)';
        n = 1; 
        return ;
    end
    
    %   Recursively Defined
    [parent1,child1,cluster1,permutation1,n1] = hierarchical_partition_pca(coord(index1,:), minSize);
    [parent2,child2,cluster2,permutation2,n2] = hierarchical_partition_pca(coord(index2,:), minSize);

    %   Combine the two subtrees
    n = n1 + n2 + 1;
    permute = zeros(n, 1);
    permute(index1) = permutation1;    
    permute(index2) = length(index1) + permutation2;
    
    %   Post-Order representation of the partition tree.
    parent = [parent1(1:end-1); n; parent2(1:end-1)+n1;  n;  0];
    children = [child1; n1+child2; n1, n1+n2];
    cluster = [cluster1; cluster2 + length(index1); 1,N];
end







