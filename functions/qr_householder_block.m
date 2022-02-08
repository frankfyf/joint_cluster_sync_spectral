%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The blockwise column-pivoted QR factorization by householder
% transformation
%
% Inputs:
%      X: input matrix
%      d: the dimension of the block
% Outputs: We output three matrices, Q, R and Pi that satisfy X * Pi = Q * R
%      Q: the orthogonal matrix 
%      R: the upper-triangular matrix
%      Pi: the permutation matrix that indicates the pivot selection

function [ Q, R, Pi ] = qr_householder_block(X, d)

% Parameters
[K, n] = size(X);
K = K/d;
n = n/d;

% Initialization
Pi = eye(n*d);
Q = eye(K*d);
R = X;

% QR factorization
for t = 1:min(K,n)
    
    % Determine the pivot block column
    col_l2 = zeros(1,n-t+1);
    for i = t:n
        col_l2(i-t+1) = norm(R((t-1)*d+1:end, (i-1)*d+1:i*d), 'fro');
    end
    [tmp_max, tmp_id] = max(col_l2);
    if tmp_max == 0
        break
    end
    if tmp_id ~= 1
        tmp_id = t+tmp_id-1;
        R(:, [(t-1)*d+1:t*d, (tmp_id-1)*d+1:tmp_id*d]) = R(:, [(tmp_id-1)*d+1:tmp_id*d, (t-1)*d+1:t*d]);
        Pi(:, [(t-1)*d+1:t*d, (tmp_id-1)*d+1:tmp_id*d]) = Pi(:, [(tmp_id-1)*d+1:tmp_id*d, (t-1)*d+1:t*d]);
    end
    
    % Householder transformation 
    for j = 1:d
        pivot = (t-1)*d+j;
        u = R(pivot:end, pivot);
        u(1,1) = u(1,1) + sign(R(pivot,pivot))*norm(R(pivot:end,pivot),2);
        v = u./norm(u,2);
        Q_l = eye(K*d-pivot+1) - 2*v*(v');
        R(pivot:end,pivot:end) = Q_l*R(pivot:end,pivot:end);
        Q(:,pivot:end) = Q(:,pivot:end)*Q_l;
    end
end

end

