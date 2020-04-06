function [D, U, B] = Mat2Hss_SPD_tol(A, htree, tol)
%   SPD HSS Method 1
%
%   SPD HSS approximation of an SPD matrix A with FIXED RELATIVE ERROR THRESHOLD 
%   for off-diagonal blocks. 
%
%   Input: 
%       A:      target kernel matrix
%       htree:  partition tree for the HSS approximation 
%       tol:    relative error threshold
%               (can be an array defiing different tol in different levels)
%   
%   Output: 
%       D, U, B: triplets that represent the constructed HSS approximation.
%           D --- the diagonal blocks defined in leaf level.     
%           U --- the nested basis/transfer matrices associated with each node
%           B --- the intermediate blocks associated with each off-diagonal
%                 block
%       One off-diagonal block A_{ij} defined by clusters i and j are
%       approximated as 
%           A_{ij} \approx Ui Bij Uj^T. 
%


    %   Output Initialization
    D = {[]}; U = {[]}; B = {[]}; 
    
    %   Tree info
    N = size(A, 1);
    parent   = htree.parent;
    children = htree.children;
    level    = htree.level;
    root     = htree.root;
    nlevel   = length(level);
    mcluster = htree.mcluster;
    
    %   Auxiliary info
    tmpE = {}; 
    effect_idx_per_node = cell(length(parent),1);
    effect_rowidx = true(N,1);   
    effect_rowidx_lvl = effect_rowidx;
    
    %   Parse the input tol
    if (isscalar(tol))
        tol = tol * ones(nlevel, 1);
    end

    %   Construction goes from the leaf level to the level below the root
    for i = nlevel : -1 : 2
        for j = 1 : length(level{i})
            node = level{i}(j);
            ind = mcluster(node,1):mcluster(node,2);
            xorind = effect_rowidx_lvl;
            xorind(ind) = false;               %   the complementary of idx
            
            %   For Leaf node 
            if (isnan(children(node,1)))                
                %   Store the diagonal block at leaf level
                D{node} = A(ind, ind);
             
                %   Find the proper sub-eigenspace of D{node};
                [V, e] = eig(D{node});
                e = diag(e);
                A(ind, xorind) = V' * A(ind, xorind);  %onsite modification
                score = sum(A(ind, xorind).^2, 2);
                [score,I] = sort(score, 'descend');
                score = cumsum(score);
                rr = find(score > (1-tol(i)^2) * score(end), 1, 'first');                
                I1 = I(1:min(rr, size(V,1)));
                Ic = I(min(rr, size(V,1))+1:end);
                
                %   Generator U and diagonal block
                U{node} = V(:, I1);
                tmpE{node} = e(I1); 
                
                %   Remaining Rows' index update. 
                effect_idx_per_node{node} = ind(1)-1+I1;               
                effect_rowidx(ind(1)-1+Ic) = false;
                
            %   For Non-Leaf Node
            else 
                c1 = children(node,1);
                c2 = children(node,2);
                ind1 = effect_idx_per_node{c1};
                ind2 = effect_idx_per_node{c2};
                ind12 = [ind1; ind2];
                
                %   Intermediate block defined by two children (c1,c2)
                B{node} = A(ind1, ind2);
                
                %   Find generator Rc1, Rc2. 
                tmpD = [diag(tmpE{c1}), B{node}; B{node}', diag(tmpE{c2})];
                [V,e] = eig(tmpD);
                e = diag(e);
                A(ind12, xorind) = V'*A(ind12, xorind);
                score = sum(A(ind12, xorind).^2, 2);
                [score,I] = sort(score, 'descend');
                score = cumsum(score);
                rr = find(score > (1-tol(i)^2) * score(end), 1, 'first');               
                I1 = I(1:min(rr, size(V,1)));
                Ic = I(min(rr, size(V,1))+1:end);
                              
                %   Generator (Rc1;Rc2) and diagonal block. 
                U{node} = V(:, I1);
                tmpE{node} = e(I1); 
                
                %   Remaining Rows' index update. 
                effect_idx_per_node{node} = ind12(I1);
                effect_rowidx(ind12(Ic)) = false; 
            end            
        end
        %   eliminate the ignored rows. 
        effect_rowidx_lvl = effect_rowidx;
        
        %   Compress the corresponding block columns by the basis/transfer 
        %   matrices constructed above. 
        for j = 1 : length(level{i})
            node = level{i}(j);
            ind = mcluster(node,1):mcluster(node,2);
            xorind = effect_rowidx_lvl;
            xorind(ind) = false;
            if (isnan(children(node,1))) 
                A(xorind, effect_idx_per_node{node}) = A(xorind, ind) * U{node};
            else
                c1 = children(node,1);
                c2 = children(node,2);
                ind12 = [effect_idx_per_node{c1}; effect_idx_per_node{c2}];
                A(xorind, effect_idx_per_node{node}) = A(xorind, ind12) * U{node};
            end
        end
    end
    
    %   Intermediate block defined by the two children of the root
    c1 = children(root,1);
    c2 = children(root,2);
    ind1 = effect_idx_per_node{c1};
    ind2 = effect_idx_per_node{c2}; 
    B{root} = A(ind1, ind2);
end