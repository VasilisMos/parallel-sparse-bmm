function M = rand_matrix(n,p)
    M = double(rand(n) > (1-p));
end