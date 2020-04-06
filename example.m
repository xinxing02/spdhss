%%  Loading of the particle configuration (uniform)
N = 4000; 
coord = rand(N, 2);

%%  Hierarchical partition of the point set 
minSize = 50;
[parent, children, cluster, perm, ntree] = hierarchical_partition_pca(coord, minSize);
coord(perm,:) = coord;
htree.root     = length(parent);             %   root index
htree.parent   = parent;                     %   parent list
htree.children = children;                   %   children list
htree.cluster  = cluster;                    %   point cluster list
htree.mcluster = cluster;                    %   matrix row/column cluster list (usually each point corresponds to one row, sometimes one point can correspond to multiple rows)
htree.level    = nodes_of_level(parent);     %   the list of nodes in each level. 
htree.leafnode = find(isnan(children(:,1))); %   leaf node list

%%  Kernel matrix formulation
kernel = @(coord)(multiquadric(coord, .5, 50));
A = kernel(coord) + 1e-5*eye(size(coord,1));

%%  SPD HSS approx. construction (defined by triplets (D,U,B) )
% r = 30;
% [D, U, B] = Mat2Hss_SPD2_rank(A, htree, r);
% [D, U, B] = Mat2Hss_SPD_rank(A, htree, r);
tol = 5e-2;
[D, U, B] = Mat2Hss_SPD2_tol(A, htree, tol);
% [D, U, B] = Mat2Hss_SPD_tol(A, htree, tol);

%%  HSS Cholesky of the HSS approx. (defined by triplets (Q, L, Idx))
[Q, L, Idx] = HssCholesky(D, U, B, htree);

%%  Check the positive definiteness of the constructed HSS approx.
hssA = Hss2Mat(D, U, B, htree);     % dense form of the HSS approx.
eighssA = sort(eig(hssA));
min(eighssA)

%%  Preconditioner for solving Ax = b (compared to block jacobi)

%   Block Jacobi preconditioner
bjA = zeros(size(A,1));
for i = 1 : length(htree.leafnode)
    idx = htree.mcluster(htree.leafnode(i),1) : htree.mcluster(htree.leafnode(i),2);
    bjA(idx, idx) = inv(A(idx, idx));
end

%   CG solving
rhs = randn(N, 1);
hss_prec = @(x)(HssCholSolve(Q, L, Idx, x, htree));
bj_prec = @(x)(bjA * x);

[~, flag1, ~, iter1, resvec1] = pcg(A, rhs, 1e-6, 600);
[~, flag2, ~, iter2, resvec2] = pcg(A, rhs, 1e-6, 600, hss_prec);
[~, flag3, ~, iter3, resvec3] = pcg(A, rhs, 1e-6, 600, bj_prec);

fig = figure('DefaultAxesFontSize',16);
semilogy(resvec1 / resvec1(1), '-r');
hold on
semilogy(resvec2 / resvec2(1), '--b');
semilogy(resvec3 / resvec3(1), '-.');
legend('Unprecond.', 'SPDHSS precond.', 'BJ precond.');

