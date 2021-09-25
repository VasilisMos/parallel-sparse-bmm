function store_sparse_matrix(filename,A)
   [i,j,k] = find(A);
   m_size = size(A);

   fileID = fopen(filename,'w');
   fprintf(fileID,'%%%%MatrixMarket matrix coordinate real general\n');
   fprintf(fileID,'%d %d %d\n',m_size(1),m_size(2),nnz(A));
   fprintf(fileID,'%d %d %d\n',[i j k]');
   fclose(fileID);
end