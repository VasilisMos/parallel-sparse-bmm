%Sample matlab Code to verify result
n=5e3;
d = 2;
A = sprand( n, n, d/n) > 0;
B = sprand( n, n, d/n) > 0;
tic; C = (A*B) > 0; toc

function logic_test()
    n = 5;

    disp('Integer')
    x = rand_bool(n);
    y = rand_bool(n);
    tic; res = x*y; toc
    
    disp('Logical')
    x_l = logical(x);
    y_l = logical(y);
    tic; res_l = x_l*y_l; toc

end

function x = rand_bool(n)
    x = randi(2,[n n])-1;
end