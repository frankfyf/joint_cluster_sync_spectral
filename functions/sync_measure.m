%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the synchronization error
%
% Inputs:
%      V: The concatenation of the true orthorgonal transforms of each node, 
%         which is of size nd by d
%      V_r: The concatenation of the identified orthogonal transforms
%      id_true: the true cluster memberships, the value should range from 1
%               to K, where K is the number of clusters
% Outputs:
%      error_r: the synchronization error

function [ error_r ] = sync_measure( V, V_r, id_true)

%%% Parameters
n = numel(id_true); % data size
d = size(V,1)/n; % the dimension of the orthogonal transformation
K = max(id_true); % the number of clusters

%%% Compute the error for each cluster
error_r = zeros(1, K);
V = mat2cell(V, d*ones(1, n), d); 
V_r = mat2cell(V_r, d*ones(1, n), d);

for k = 1:K
    
    % Get the transforms of the nodes in the k-th cluster
    V_k = V(id_true == k, :);
    V_r_k = V_r(id_true == k, :);
    
    % Find the optimal alignment
    [U_tmp, ~, V_tmp] = svd(cell2mat(V_r_k)'*cell2mat(V_k));
    opt_ortho = U_tmp*V_tmp';
    
    % Compute the synchronization error of each node
    error_r_k = zeros(1, numel(V_k));
    for i = 1:numel(V_k)
        error_r_k(i) = norm(V_k{i} - V_r_k{i}*opt_ortho, 'fro');
    end
    error_r(k) = max(error_r_k);
end

% Take the logarithm of the maximum error 
error_r = log(max(error_r)/sqrt(d));

end

