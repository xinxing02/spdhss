function [D, U, B] = Mat2Hss_SPD2_tol(A, htree, tol)
%   SPD HSS Method 2
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

    %   Output Initialization
    D = {[]}; U = {[]}; B = {[]}; 
    
    %   Tree info
    parent   = htree.parent;
    children = htree.children;
    level    = htree.level;
    root     = htree.root;
    nlevel   = length(level);
    mcluster = htree.mcluster;
    
    %   Axuliary Variables
    V = {}; S = {}; Q = {}; N = {};
    tmpDiag = {}; tmpQ = {};
    
    %   Location Auxiliary arrays
    effect_idx_per_node = cell(length(parent),1);
    effect_rowidx = true(size(A,1),1);
    effect_rowidx_lvl = effect_rowidx;
    
    %   Parse the input tol
    if (isscalar(tol))
        tol = tol * ones(nlevel,1);
    end
    
    %   Construction goes from the leaf level to the level below the root
    %   In a non-leaf level, Mij is stored in Aij (in-place modification)
    
    for i = nlevel : -1 : 2   
        %   First step: scale the off-diagonal blocks in the level i. 
        %   the scaled off-diagonal block Cij = Si^{-1}AijSj^{-T} is stored in Aij 
        for j = 1 : length(level{i})
            node = level{i}(j);
            idx = mcluster(node,1):mcluster(node,2);
            xorind = effect_rowidx_lvl;
            xorind(idx) = false;    %   the complementary of idx.
          
            %   For Leaf node 
            if (isnan(children(node,1))) 
                %   Store the diagonal block at leaf level
                D{node} = A(idx, idx);
                S{node} = chol(D{node}, 'lower');

                %   Scale the off-diagonal blocks using S{node}
                tmpS = inv(S{node});
                A(idx, xorind) = tmpS * A(idx, xorind);
                A(xorind, idx) = A(idx, xorind)';
                
            %   For Non-Leaf Node    
            else
                c1 = children(node,1);
                c2 = children(node,2);
                ind1 = effect_idx_per_node{c1};
                ind2 = effect_idx_per_node{c2};
                
                %   Intermediate block defined by two children (c1,c2)
                B{node} = A(ind1, ind2); 
                
                %   SVD of the intermediate block
                [Q{c1}, tmp, Q{c2}] = svd(B{node});
                N{node} = diag(tmp);
                r1 = size(Q{c1},1);
                r2 = size(Q{c2},1);
                rr = min(r1, r2);
                
                %   Temporary matrix
                tmpDiag{node} = [(1+N{node}); (1-N{node}); ones(r1+r2-2*rr,1)];
                tmpQ{node} = [Q{c1}(:,1:rr)/sqrt(2), Q{c1}(:,1:rr)/sqrt(2), Q{c1}(:,(rr+1):r1), zeros(r1, r2-rr); 
                    Q{c2}(:,1:rr)/sqrt(2), -Q{c2}(:,1:rr)/sqrt(2), zeros(r2, r1-rr), Q{c2}(:,(rr+1):r2)];
                
                %   Scale the off-diagonal block 
                A([ind1,ind2], xorind) = diag(tmpDiag{node}.^(-1/2)) * tmpQ{node}' * A([ind1,ind2], xorind);
                A(xorind, [ind1,ind2]) = A([ind1,ind2], xorind)';            
            end 
        end
        
        %   Second step: find the optimal ``basis matrix'' V_i to compress
        %   the scaled off-diagonal blocks C_{ii^c}.
        for j = 1 : length(level{i})
            node = level{i}(j);
            idx = mcluster(node,1):mcluster(node,2);
            xorind = effect_rowidx_lvl;
            xorind(idx) = false;    %   the complementary of idx.
            
            %   For Leaf node 
            if (isnan(children(node,1))) 
                V{node} = randSVDtol(A(idx, xorind), tol(i));
                U{node} = S{node} * V{node};
                
            %   For Non-Leaf Node
            else
                c1 = children(node,1);
                c2 = children(node,2);
                ind1 = effect_idx_per_node{c1};
                ind2 = effect_idx_per_node{c2};  
                
                %   Compute the basis matrices
                V{node} = randSVDtol(A([ind1,ind2], xorind), tol(i));
                U{node} = tmpQ{node} * diag(tmpDiag{node}.^(1/2)) * V{node};     
            end            
        end
        
        %   Third step: Update the in-place Cij to the intermediate matrix 
        %               Mij = Vi'CijVj
        for j = 1 : length(level{i})
            node = level{i}(j);
            idx = mcluster(node,1):mcluster(node,2);
            xorind = effect_rowidx_lvl;
            xorind(idx) = false;    %   the complementary of idx.                      
            if (isnan(children(node,1)))
                tmp_r = size(V{node},2);
                A(idx(1):idx(1)+tmp_r-1, xorind) = V{node}' * A(idx, xorind);
                A(xorind, idx(1):idx(1)+tmp_r-1) = A(idx(1):idx(1)+tmp_r-1, xorind)';
                
                %   Effective rows after the in-place modification
                effect_idx_per_node{node} = idx(1) : idx(1)+tmp_r -1;
                effect_rowidx(idx(1)+tmp_r:1:idx(end)) = false;
            else
                c1 = children(node,1);
                c2 = children(node,2);               
                ind1 = effect_idx_per_node{c1};
                ind2 = effect_idx_per_node{c2};
                
                tmp_r = size(V{node},2);  
                A(idx(1):idx(1)+tmp_r-1, xorind) = V{node}' * A([ind1,ind2], xorind);
                A(xorind, idx(1):idx(1)+tmp_r-1) = A(idx(1):idx(1)+tmp_r-1, xorind)';
                
                %   Effective rows after the in-place modification
                effect_idx_per_node{node} = idx(1) : idx(1)+tmp_r -1;
                effect_rowidx(idx(1)+tmp_r:1:idx(end)) = false;
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