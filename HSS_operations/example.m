%
%   Example to use the standard HSS matrix functions by exponential kernel
%%  Coordinates Input
npts = 5000;
coord = npts^(1/3) * rand(npts, 3);

%%  Partition of the system 
minSize = 100;
[parent, children, cluster, perm, ntree] = hierarchical_partition_pca(coord, minSize);
coord(perm, :) = coord;
htree.parent   = parent;
htree.children = children; 
htree.cluster  = cluster;
htree.mcluster = cluster;
htree.level    = nodes_of_level(parent); 
htree.leafnode = find(isnan(children(:,1)));
htree.root     = find(~(parent>0));

%%  Kernel Function & Coordinates
kernel = @(x)(multiquadric(x, 1/2, 0.5));
A = kernel(coord);

%%  Standard HSS construction
% r = 50;
% [D, U, B] = Mat2Hss_rank(A, htree, r);
tau = 1e-2;
[D, U, B] = Mat2Hss_tol(A, htree, tau);


%%  Matrix Reconstruction from HSS matrix
hssA = Hss2Mat(D, U, B, htree);
norm(A - hssA, 'fro') / norm(A, 'fro')


%%  HSS matrix vector multiplication test
vec = rand(size(A,1),1);
hss_Ax = HssMatVecProduct(D, U, B, vec, htree);
norm(A * vec - hss_Ax, 'fro') / norm(A * vec, 'fro')


%%  Rank for each off-diagonal block
rr = HssRank(U, htree.mcluster);


%%  ULV decomposition of HSS matrix
[Q,Lr,Lc,Idx] = HssULV_Sym(D, U, B, htree);
hss_ulv_Ax = HssUlvProduct_Sym(Q, Lr, Lc, Idx, vec, htree);


%%  Solving by ULV decomposition
rhs = HssMatVecProduct(D, U, B, vec, htree);
u = HssUlvSolve_Sym(Q, Lr, Lc, Idx, rhs, htree);
norm(u - vec, 'fro') / norm(vec)


%%  Preconditioner for solving Ax = b (mostly not working)
%   Block Jacobi preconditioner
bjA = zeros(size(A,1));
for i = 1 : length(htree.leafnode)
    idx = htree.mcluster(htree.leafnode(i),1) : htree.mcluster(htree.leafnode(i),2);
    bjA(idx, idx) = inv(A(idx, idx));
end

%   GMRES solve
rhs = randn(size(A,1), 1);
hss_prec = @(x)(HssUlvSolve_Sym(Q, Lr, Lc, Idx, x, htree));
bj_prec = @(x)(bjA * x);

[~, flag1, ~, iter1, resvec1] = gmres(A, rhs, [], 1e-6, size(A,1));
[~, flag2, ~, iter2, resvec2] = gmres(A, rhs, [], 1e-6, size(A,1), hss_prec);
[~, flag3, ~, iter3, resvec3] = gmres(A, rhs, [], 1e-6, size(A,1), bj_prec);

fig = figure('DefaultAxesFontSize',16);
semilogy(resvec1 / resvec1(1), '-r');
hold on
semilogy(resvec2 / resvec2(1), '--b');
semilogy(resvec3 / resvec3(1), '-.');
legend('Unprecond.', 'SPDHSS precond.', 'BJ precond.');