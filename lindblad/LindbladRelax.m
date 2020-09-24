function RhoOut = LindbladRelax(RhoIn,t,relax)
%Kraus transformation for the amplitude channel, that simulates relaxation
%of the qubit.
%input: RhoIn, t, amplitude
%these are: input density matrix, time, decay rate
%ouput: RhoOut
%this is: output density matrix

Rho=RhoIn;

nqbits=log(length(Rho))/log(2);

Ea0=[1 0;0 sqrt(exp(-t*relax))];
Ea1=[0 sqrt(1-exp(-t*relax));0 0];

for i=1:nqbits
    M0=kron(eye(2^(i-1)),kron(Ea0,eye(2^(nqbits-i))));
    M1=kron(eye(2^(i-1)),kron(Ea1,eye(2^(nqbits-i))));
    Rho=M0*Rho*M0'+M1*Rho*M1';
end

RhoOut=Rho;