n = (1:10) * 1e6;

for i=1:length(n)
    fprintf("Iteration %d/%d, n=%d\n",i,length(n),n(i))
    for k=1:5
        feval("generateDatasets",n(i));
    end
end



function mat_benchmark()
    d = 2;
    t = [];
    sizes = 1e3:5e3:5e5; % 1e3:5e3:5e6;
    waitbar(0,"Matlab boolean matrix multiplication benchmark..");
    for n = sizes
      A = sprand( n, n, d/n) > 0;
      B = sprand( n, n, d/n) > 0;
      tic; C = (A*B) > 0; t1 = toc; t = [t t1];
      waitbar(n/max(sizes));
    end
    
    figure;
    plot(sizes,t);
    title("Boolean Matrix Multiplication - MATLAB Implementation");
    xlabel("N"); ylabel("Time (sec)");
end