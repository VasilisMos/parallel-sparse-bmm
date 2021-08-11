clear all;

relative_folder = "../datasets/test/";
fnameA = relative_folder + "A_test.mtx";
fnameB = relative_folder + "B_test.mtx";
fnameC = relative_folder + "C_result.mtx";

differences = validate_results(fnameA,fnameB,fnameC);

function differences = validate_results(filenameA,filenameB,filenameC)
    A = mmread(filenameA);
    B = mmread(filenameB);
    C = mmread(filenameC);

    C_gold = double(logical(A*B));

    fprintf("Sparse Comparison:\n");
    fprintf("BMM Sparse Validation Test:");
    if isequal(C, C_gold);
        fprintf(" PASSED\n");
    else
        fprintf(" FAILED\n");
    end


    if max([size(A) size(B) size(C)]) > 5e3
        fprintf("Matrices size is too big for dense validation\n");
        differences = C-C_gold;
        return;
    end

    A = full(A);
    B = full(B);
    C = full(C);
    
    C_gold = double(logical(A*B));

    fprintf("Dense Comparison:");
    if isequal(C, C_gold);
        fprintf("Boolean Matrix Multiplication Validation Test PASSED\n");
    else
        fprintf("Boolean Matrix Multiplication Validation Test FAILED\n");
    end

    differences = C-C_gold;
end

%A = mmread('../datasets/A_small.mtx');
%A = full(A);
%C = mmread('../datasets/C_result.mtx');
%C = full(C);
%C_gold = double(logical(A*A));

%isequal(C_gold,C)

