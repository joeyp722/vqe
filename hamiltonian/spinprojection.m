function [ n ] = spinprojection(N,K)
%The projection matrix filtering out states describing too many electrons
proj=eye(2^K);
for j = -K:K
    if abs(j)>N
        proj=proj*((numberspin(K,'up')-numberspin(K,'down'))-j*eye(2^K))/(N-j);
    end
end
n=proj;