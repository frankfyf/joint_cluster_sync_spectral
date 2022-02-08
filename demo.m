%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo of the numerical experiment in the paper "A Spectral 
% Method for Joint Community Detection and Orthogonal Group 
% Synchronization." by Yifeng Fan, Yuehaw Khoo, and Zhizhen Zhao.
%
% Yifeng Fan (yifengf2@illinois.edu), Feb 2022

clear 
addpath(genpath('./'))
rng('default')

%%% Parameter setting
m_list = [100, 200, 300]; % The cluster sizes
d = 3; 
p = 0.8;
q = 0.2;
K = numel(m_list); 

%%% Generate the observation matrix A
[A, V, M, id_true] = gen_observation(m_list, p, q, d); 

%%% Apply the spectral method
[id, V_r] = sync_spectral(A, K, d);

%%% Evaluation
error_c = acc_measure(id, id_true) % Accuracy of clustering
error_r = sync_measure( V, V_r, id_true) % Synchronization error (log-scale)
    

