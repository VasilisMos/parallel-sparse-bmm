function A_pad = zero_padding(A,n,type)

    m = floor(log2(n));
    index_max = ceil(n/m) * m;

    if(index_max>n) %Needs Padding
        fprintf("Doing padding of %d\n",index_max-n);
        if(type==0)
            A_pad = [A zeros(n,index_max-n)];
        else
            A_pad = [A; zeros(index_max-n,n)];
        end
    else
        A_pad = A;
    end
end