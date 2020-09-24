function [ n ] = projection(N,K)
%The projection matrix filtering out states describing too many electrons
proj=eye(2^K);
for j = 0:K
    if j~=N
        proj=proj*(number(K)-j*eye(2^K))/(N-j);
    end
end
n=proj;
        