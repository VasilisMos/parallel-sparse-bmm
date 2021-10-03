%% Matlab Implementation of the pseudocode Algorithm given in
% https://louridas.github.io/rwa/assignments/four-russians/
% Input:
%  A: Dense nxn matrix
%  B: Dense nxn matrix
% Output:
%  C: Dense nxn matrix, result of boolean matrix multiplication of A and B

function C = four_russians(A_in,B_in)
    n = size(A_in,1);
    m = floor(log2(n));

    A = zero_padding(A_in,n,0);
    B = zero_padding(B_in,n,1);

    fprintf("Four Russians, (n,m)=(%d,%d)\n",n,m);

    C = zeros(n);

    for i = 1:ceil(n/m)
        Rs = zeros(2^m,n);
        bp=1;
        k=0;
        block_size = m;

        chunk = (i-1)*block_size:i*block_size-1; chunk = chunk+1; %Matlab/Octave Indexing

        Ai=A(:,chunk);
        Bi=B(chunk,:);

        num = @(v) (v)*( 2.^(size(v,2)-1:-1:0) )' +1;

        for j = 2:2^m+1
            Rs(j,:) = or(Rs(j-2^k,:),Bi(m+1 - (k+1),:));
            if(bp==1)
                bp=j+1;
                k=k+1;
            else
                bp=bp-1;
            end
        end

        Ci = zeros(n);
        for j=1:n 
            Ci(j,:) = Rs(num(Ai(j,:)),:);
        end
        C = or(C,Ci);
    end
end














% function C = four_russians(A_in,B_in)
%     n = size(A_in,1);
%     m = floor(log2(n));

%     A = zero_padding(A_in,n,0);
%     B = zero_padding(B_in,n,1);

%     fprintf("Four Russians, (n,m)=(%d,%d)\n",n,m);

%     C = zeros(n);

%     for i = 1:ceil(n/m)
%         Rs = zeros(2^m,n);
%         bp=1;
%         k=0;
%         block_size = m;

%         chunk = (i-1)*block_size:i*block_size-1; chunk = chunk+1 %Matlab/Octave Indexing

%         Ai=A(:,chunk);
%         Bi=B(chunk,:);

%         num = @(v) (v)*( 2.^(size(v,2)-1:-1:0) )' +1;

%         for j = 2:2^m+1
%             % fprintf("j=%d, k=%d\n",j,k);
%             Rs(j,:) = or(Rs(j-2^k,:),Bi(k+1,:));
%             if(bp==1)
%                 bp=j+1;
%                 k=k+1;
%             else
%                 bp=bp-1;
%             end
%         end

%         Ci = zeros(n);
%         for j=1:n 
%             Ci(j,:) = Rs(num(Ai(j,:)),:);
%         end
%         % Ci
%         C = or(C,Ci);
%     end
% end