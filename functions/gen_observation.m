%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the observation matrix A. We assume the first m_1 nodes form the first cluster,
% the following m_2 nodes form the second one, and so on...
%
% Inputs:
%      m_list: The list of cluster sizes, should be of length K, where K is
%              the number of clusters
%      p: The probability of connection within clusters
%      q: The probability of connection across different clusters
%      d: the dimension of orthogonal transformation in SO(d)
% Outputs:
%      A: The observation matrix, of size nd by nd
%      V: The concatenation of the true orthorgonal transforms of each node, is of size nd by d
%      A_true: The ground truth of A
%      id_true: The true cluster membership

function [ A, V, A_true, id_true] = gen_observation(m_list, p, q, d)

% Parameters
n = sum(m_list); % The total number of data points
K = numel(m_list); % The number of blocks

% Generate the true cluster memberships
id_true = [];
for k = 1:K
    id_true = [id_true, k * ones(1, m_list(k))];
end

% Generate the orthogonal group elements
V = repmat(eye(d), n, 1); % We assume O_i = I_d

% Generate the observations matrix
A = cell(K,K); % We treat A as a K*K block matrix.  
A_true = zeros(n*d, n*d);
count = 0;
for i = 1:K
    for j = i:K
        if j == i % The diagonal block
            
            tmp = V(count+1:count+d*m_list(i),:)*V(count+1:count+d*m_list(i),:)';
            A_true(count+1:count+d*m_list(i), count+1:count+d*m_list(i)) = tmp;
            for l1 = 1:m_list(i)
                for l2 = l1+1:m_list(i)
                    if rand(1) > p
                        tmp((l1-1)*d+1:l1*d, (l2-1)*d+1:l2*d) = zeros(d,d);
                    end
                end
            end 
            tmp = triu(tmp);
            tmp = tmp + tmp';
            tmp(1:m_list(i)*d+1:end) = 1;
            A{i,j} = tmp;
            
        else 
            
            tmp = zeros(d*m_list(i), d*m_list(j));
            for l1 = 1:m_list(i)
                for l2 = 1:m_list(j)
                    if rand(1) <= q
                        [u, ~, v] = svd(randn(d,d));
                        tmp((l1-1)*d+1:l1*d, (l2-1)*d+1:l2*d) = u*v';
                    else
                        tmp((l1-1)*d+1:l1*d, (l2-1)*d+1:l2*d) = zeros(d,d);
                    end
                end
            end
            A{i,j} = tmp;   
            A{j,i} = tmp';
            
        end
    end
    count = count + d*m_list(i);
end
A = cell2mat(A);
A = sparse(A);
    
end

