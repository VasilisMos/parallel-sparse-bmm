n=12;
m = floor(log2(n));

for i = 1:1
    Rs = zeros(2^m,n);
    bp=1;
    k=0;
    block_size = m;

    chunk = (i-1)*block_size:i*block_size-1; chunk = chunk+1; %Matlab/Octave Indexing

    for j = 1:2^m
        printf("j=%s, k=%d\n",dec2bin(j),k);
        printf("m-k=%d\n",m-k);
        % printf("j-2^k=%d\n",j-2^k);
        fprintf("\n");
        if(bp==1)
            bp=j+1;
            k=k+1;
        else
            bp=bp-1;
        end
    end
end