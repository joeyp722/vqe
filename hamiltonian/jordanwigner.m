function [ a ] = jordanwigner(j,N,sign)
%defining the ospin matrices
oneop=[1 0; 0 1];
xop=[0 1; 1 0];
yop=[0 -1i; 1i 0];
zop=[1 0; 0 -1];
%calculating the JordanWigner Transform matrices by tensor products
aj=oneop;
if j>2
    for k = 1:(j-2)
    aj=kron(aj,oneop);
    end
end
if j==1
    aj=1/2*(xop+sign*1i*yop);
else
    aj=kron(aj,1/2*(xop+sign*1i*yop));
end
for k = (j+1):N
    aj=kron(aj,zop);
end
a=aj;

