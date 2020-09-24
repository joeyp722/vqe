function RhoOut = LindbladAmplitude(RhoIn,t,amplitude)
%Kraus transformation for the X-pauli channel, that simulates amplitude
%noise
%of the qubit.
%input: RhoIn, t, amplitude
%these are: input density matrix, time, decay rate
%ouput: RhoOut
%this is: output density matrix

Rho=RhoIn;

nqbits=log(length(Rho))/log(2);

I=[1 0;0 1];
X=[0 1;1 0];

for i=1:nqbits
    M0=kron(eye(2^(i-1)),kron(I,eye(2^(nqbits-i))));
    M1=kron(eye(2^(i-1)),kron(X,eye(2^(nqbits-i))));
    c0=0.5*(1+exp(-amplitude*t));
    c1=0.5*(1-exp(-amplitude*t));
    Rho=c0*M0*Rho*M0'+c1*M1*Rho*M1';
end

RhoOut=Rho;