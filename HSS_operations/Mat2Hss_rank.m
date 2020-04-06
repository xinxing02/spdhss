function [D, U, B] = Mat2Hss_rank(A, htree, r)
%   Standard HSS approximation 
%   For symmetric matrix A with a given rank r for the off-diagonal blocks. 
%   
%   Modification is done directly on the matrix A. 
%   Simulated Parallel Calculation

    N = size(A, 1);
    parent   = htree.parent;
    children = htree.children;
    level    = htree.level;
    root     = htree.root;
    mcluster = htree.mcluster;
    nlevel   = length(level);
    
    D = {}; U = {}; B = {}; 
    index = cell(length(parent),1);
    
    if (isscalar(r))
        r = r * ones(nlevel,1);
    end
    rowidx = true(N,1);
    rowidx_lvl = rowidx;
    for i = nlevel : -1 : 2          
        %   Find the projection matrix U
        for j = 1 : length(level{i})
            node = level{i}(j);
            ind = mcluster(node,1):mcluster(node,2);
            xorind = rowidx_lvl;
            xorind(ind) = false;           
            %   For Leaf node 
            if (isnan(children(node,1)))                
                %   Store the diagonal block at leaf level
                D{node} = A(ind, ind);
                %   Generator U and diagonal block
                U{node} = randSVD(A(ind, xorind), r(i));           
            %   For Non-Leaf Node
            else 
                c1 = children(node,1);
                c2 = children(node,2);
                ind12 = [index{c1}; index{c2}];                 
                %   Intersection Block
                B{node} = A(index{c1}, index{c2});                              
                %   Generator (Rc1;Rc2) and diagonal block. 
                U{node} = randSVD(A(ind12, xorind), r(i));                          
            end            
        end
        
        %   Compress the columns by the projection obtained above. 
        for j = 1 : length(level{i})
            node = level{i}(j);
            ind = mcluster(node,1):mcluster(node,2);
            xorind = rowidx_lvl;
            xorind(ind) = false;   
            if (isnan(children(node,1)))
                %   update row index
                index{node} = (ind(1): (ind(1)+size(U{node},2)-1))';
                rowidx( (ind(1)+size(U{node},2)):ind(end) ) = false;
                %   compress the rows and columns
                A(index{node}, xorind) = U{node}' * A(ind, xorind);
                A(xorind, index{node}) = A(index{node}, xorind)';
            else
                c1 = children(node,1);
                c2 = children(node,2);
                ind12 = [index{c1}; index{c2}];    
                %   update row index
                index{node} = ind12(1:size(U{node},2));
                rowidx( ind12((size(U{node},2)+1) : end) ) = false;
                %   compress the rows and columns
                A(index{node}, xorind) = U{node}' * A(ind12, xorind);
                A(xorind, index{node}) = A(index{node}, xorind)';                
            end
        end
        
        %   update row index
        rowidx_lvl = rowidx;
    end
    
    %   ROOT OF THE TREE
    c1 = children(root,1);
    c2 = children(root,2);
    B{root} = A(index{c1}, index{c2});
end