%% Octave Validator for BMM Multiplication
% Checks for differences between Octave/Matlab Baseline Results 
% and the current C++ implementation.
% If at least one mismatch found, prints some error info and the location
% of the mismatch using spy plot

filename_gold = "../datasets/test/C_test.mtx";
filenameC = "../datasets/test/C_result.mtx" ;

% Read Result and ground truth
C = mmread(filenameC);
C_gold = mmread(filename_gold);

%Differences between C (result from C++)
dif = C-C_gold;

if nnz(dif) == 0 % No error happened
  fprintf("Boolean Matrix Multiplication |M=%d|K=%d|N=%d|nnz(C)=%d|\n",size(C,1),size(C,2),size(C,2),nnz(C))
  disp("TEST PASSED, Matrices are equal!");
else % The two matrices disagree in at least on entry
  figure; spy(dif); title("Differnces between result and ground truth (Should be white)");
  xlabel("Dim2");
  ylabel("Dim1");

  fprintf("Res: (dim,nnz) = (%d,%d)\n",size(C,1),nnz(C));
  fprintf("Ground: (dim,nnz) = (%d,%d)\n",size(C_gold,1),nnz(C_gold));
end
