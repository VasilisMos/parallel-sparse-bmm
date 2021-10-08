function F = generate_filter(n=1e6,d=1e-5)
%generate_filter - Generate a sparse Filter
%
% Syntax: F = generate_filter(n,d)
%
% Function to generate a nxn filter (Matrix) with sparsity d (0<d<1)
% To be used to perform masked Boolean Matrix Multiplication

    output = "../../datasets/test/filter.mtx";
    get_Filter = @(n,d) sprand(n,n,d);

    F = get_Filter(n,d)>1e-5;

    fprintf('Storing Filter Matrix.. ');
    tic; store_sparse_matrix(output,F); toc;

    d = 2;
    A = sprand( n, n, d/n) > 0;
    B = sprand( n, n, d/n) > 0;
    tic; C = F.*(A*B) > 0; t = toc

    fileID = fopen("../../logs/times.csv",'a');
    fprintf(fileID,"%d,%f,MATLAB,1\n",n,t);
    fclose(fileID);

    fprintf('nnz(C)=%d\n',nnz(C));

    tag = 'test';

    fprintf('Storing matrix A.. ');
    tic; store_sparse_matrix(strcat('../../datasets/',tag,'/A_', tag,'.mtx'),A); toc;
    clear A;

    fprintf('Storing matrix B.. ');
    tic; store_sparse_matrix(strcat('../../datasets/',tag,'/B_', tag, '.mtx'),B); toc;
    clear B;

    fprintf('Storing matrix C.. ');
    tic; store_sparse_matrix(strcat('../../datasets/',tag,'/C_', tag, '.mtx'),C); toc;
    clear C;
end

