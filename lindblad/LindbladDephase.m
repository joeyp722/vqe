function RhoOut = LindbladDephase(RhoIn,t,dephase)
%Kraus transformation for the dephasing channel
%of the qubit.
%input: RhoIn, t, dephase
%these are: input density matrix, time, decay rate
%ouput: RhoOut
%this is: output density matrix

Rho=RhoIn;

nqbits=log(length(Rho))/log(2);

Ed0=[1 0;0 exp(-t*dephase)];
Ed1=[0 0;0 sqrt(1-exp(-2*t*dephase))];

for i=1:nqbits
    M0=kron(eye(2^(i-1)),kron(Ed0,eye(2^(nqbits-i))));
    M1=kron(eye(2^(i-1)),kron(Ed1,eye(2^(nqbits-i))));
    Rho=M0*Rho*M0'+M1*Rho*M1';
end

RhoOut=Rho;