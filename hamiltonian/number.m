function [ n ] = number(N)
%Calculates the number operators a+a in Jordan Wigner form
number=zeros(2^N,2^N);
for i=1:N
    number=number+jordanwigner(i,N,-1)*jordanwigner(i,N,1);
end
n=number;