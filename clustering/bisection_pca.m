function [index1,index2] = bisection_pca(coord)
%
%   Use PCA to partition the point set into two subsets
%
    center = mean(coord, 1);
    coord = coord - repmat(center, size(coord,1), 1);
    coeff = pca(coord);
    coeff = coeff(:,1);
    dir = coord * coeff;
    index1 = find(dir > 0);
    index2 = find(dir <= 0);
end