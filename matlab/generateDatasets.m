%% Function generateDatasets(n)
%  Generates 3 sparse matrices A,B,C, wuth dims = [N N] each and C = A * B,
%  with '*' describing Boolean Matrix Multiplication operation
%  Each matrix is stored in '../datasets' directory in the corresponding file

function generateDatasets(n=2.4e5)
%n = 1e4;
d = 2;
A = sprand( n, n, d/n) > 0;
B = sprand( n, n, d/n) > 0;
tic; C = (A*B) > 0; t = toc

fileID = fopen("../logs/times.csv",'a');
fprintf(fileID,"%d,%f,MATLAB,1\n",n,t);
fclose(fileID);

fprintf('nnz(C)=%d\n',nnz(C));

tag = 'test';

fprintf('Storing matrix A.. ');
tic; store_sparse_matrix(strcat('../datasets/',tag,'/A_', tag,'.mtx'),A); toc;
clear A;

fprintf('Storing matrix B.. ');
tic; store_sparse_matrix(strcat('../datasets/',tag,'/B_', tag, '.mtx'),B); toc;
clear B;

fprintf('Storing matrix C.. ');
tic; store_sparse_matrix(strcat('../datasets/',tag,'/C_', tag, '.mtx'),C); toc;
clear C;
end