%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The spectral method for joint clustering and synchornizaion 
%
% Inputs:
%      A: The observation matrix of size nd * nd
%      K: The number of clusters
%      d: the dimension of orthogonal transformation in SO(d)
% Outputs:
%      id: identified cluster membership 
%      V_r: the concatenation of recovered orthogonal transforms of size
%           nd * d

function [id, V_r] = sync_spectral(A, K, d)

n = size(A,1)/d;

% Compute the top Kd eigenvector 
[Evec, ~] = eigs(A, K*d);

% Blockwise column-pivoted QR factorization by householder transforamtion 
[~, R, P] = qr_householder_block(Evec', d);
R = R*(P');

% Determine the cluster label
R_blocknorm = zeros(K,n);
for i = 1:n
    for k = 1:K
        R_blocknorm(k,i) = norm(R((k-1)*d+1:k*d, (i-1)*d+1:i*d), 'fro');
    end
end
[~, id] = max(R_blocknorm, [], 1);

% Synchronization within each cluster
V_r = cell(n,1);
for i = 1:n
    C_i = id(i);
    [U,~, V] = svd(R((C_i-1)*d+1:C_i*d, (i-1)*d+1:i*d));
    V_r{i} = U*V';
end
V_r = cell2mat(V_r);

end

