function [Q, Y0] = randSVDtol(A, tol)
%
%   Random SVD of the matrix A with specified relative threshold 
%
%   Find Q s.t. |(I-QQ')A|_F <= tol * |A|
%
%   Presently, only give the orthogonal columns out as the projector
%
if nargout == 1
    [m,n] = size(A);
    r = 5;
    p = 2;
    Q = [];
    W = 1/sqrt(n)*randn(n, r+p);
    W = A * W;
    Y = W;

    while (norm(Y,'fro')/norm(W,'fro') > tol)
        [Q1, R1, ~] = qr(Y, 0);
        r1 = min([r, find(abs(diag(R1)) > 1e-16, 1, 'last')]);        
        Q = [Q, Q1(:,1:r1)];
        [Q, ~, ~] = qr(Q, 0);
        if (size(Q,2) >= m)
            Q = eye(m);
            break;
        end
        W = 1/sqrt(n)*randn(n,r+p);
        W = A * W;
        Y = W - Q * (Q' * W);
    end
end

if nargout == 2
    [m,n] = size(A);
    r = 5;
    p = 2;
    Q = [];
    W = 1/sqrt(n)*randn(n, r+p);
    W = A * W;
    Y = W;
    Y0 = W;
    while (norm(Y,'fro')/norm(W,'fro') > tol)
        [Q1, R1, ~] = qr(Y, 0);
        r1 = min([r, find(abs(diag(R1)) > 1e-16, 1, 'last')]);
        Q = [Q, Q1(:,1:r1)];
        [Q, ~, ~] = qr(Q, 0);
        if (size(Q,2) >= m)
            Q = eye(m);
            break;
        end
        W = 1/sqrt(n)*randn(n,r+p);
        W = A * W;
        Y0 = [Y0, W];
        Y = W - Q * (Q' * W);
    end
end

end