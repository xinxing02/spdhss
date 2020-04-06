function dimrank = HssRank(U, mcluster)
%
%   Return the rank of of each compressed off-diagonal subblock 
%   with form  [dim, rank];
%
    dimrank = zeros(size(mcluster,1),2);
    for i = 1 : size(mcluster,1) - 1
        dimrank(i,:) = [mcluster(i,2)-mcluster(i,1)+1, size(U{i},2)];
    end
    dimrank(end, :) = [(mcluster(end,2)-mcluster(end,1)+1), 0];
end
