n=120;

A = rand_matrix(n,0.4);
B = rand_matrix(n,0.4);

% Normal BMM
C_ground = A*B > 0;

tic;
C = four_russians(A,B); toc;

spy(C-C_ground)



%% Test Case
% A = [1 1 0 0 0
%      0 0 1 1 1
%      1 0 0 1 0
%      1 0 0 1 1
%      1 0 1 0 1];

% B = [0 1 0 0 1
%      0 0 0 0 0
%      1 1 0 0 1
%      1 0 1 0 0
%      1 1 0 1 0];