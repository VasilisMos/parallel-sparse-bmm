function A_pad = zero_padding(A,n,type)

    index_max = ceil(n/log2(n)) * floor(log2(n));

    if(index_max>n) %Needs Padding
        if(type==0)
            A_pad = [A zeros(n,index_max-n)];
        else
            A_pad = [A; zeros(index_max-n,n)];
        end
    else
        A_pad = A;
    end
end