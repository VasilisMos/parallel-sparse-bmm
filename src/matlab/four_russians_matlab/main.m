% function main(n=12)
n=18;

fnameA = "../../../datasets/test/A_test.mtx";
fnameB = "../../../datasets/test/B_test.mtx" ;
fnameC = "../../../datasets/test/C_result.mtx";

A = rand_matrix(n,0.1);
B = rand_matrix(n,0.1);

A = mmread(fnameA);
B = mmread(fnameB);
C = mmread(fnameC);

% Normal BMM
C_ground = A*B > 0;

tic;
% C = four_russians(A,B); 
toc;

dif = C-C_ground;

if(nnz(dif))
    spy(dif)
    pause; close all;
else
    disp('Results match, TEST PASSED');
end
