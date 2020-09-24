function [ n ] = numberspin(N,state)
%Calculates the number operators a+a in Jordan Wigner form
number=zeros(2^N,2^N);
if strcmp(state,'up')
    for i=1:N/2
        number=number+jordanwigner(i,N,-1)*jordanwigner(i,N,1);
    end
else 
    for i=(N/2+1):N
        number=number+jordanwigner(i,N,-1)*jordanwigner(i,N,1);
    end
end
n=number;