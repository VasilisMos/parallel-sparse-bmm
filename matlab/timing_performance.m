%% timing_performance.m
% This function evaluates the time performance of 
% the MATLAB version of Sparse Boolean Matrix Multiplication
% (which performs normal sparse matrix multiplication followed by
% a value thresholding)
function timing_performance()
    d = 2;
    max_iters=13;
    t = [];
    sizes = 1e5:3e5:5e6;

    h = waitbar(0,"Matlab boolean matrix multiplication benchmark..");
    for n = sizes
        t1 = [];
        for iters = 1:max_iters
            A = sprand( n, n, d/n) > 0;
            B = sprand( n, n, d/n) > 0;
            tic; C = (A*B) > 0; temp = toc; t1 = [t1 temp];
        end

        t = [t median(t1)];
        waitbar(n/max(sizes));
    end
    
    figure;
    plot(sizes,t);
    title("Boolean Matrix Multiplication - MATLAB Implementation");
    xlabel("N"); ylabel("Time (sec)");
    close(h); pause; close all;

    fileID = fopen('MATLAB_times.txt','w');
    fprintf(fileID,'%s,%s\n','dim','times');
    fprintf(fileID,'%d,%.4f\n',[sizes; t]);
    fclose(fileID);
end
